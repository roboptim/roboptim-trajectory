// Copyright (C) 2009 by Thomas Moulard, AIST, CNRS, INRIA.
//
// This file is part of the roboptim.
//
// roboptim is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// roboptim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with roboptim.  If not, see <http://www.gnu.org/licenses/>.

#include <boost/format.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/spline.hh>
#include <roboptim/trajectory/trajectory-cost.hh>

#include "common.hh"

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;



int run_test ()
{
  Spline::vector_t params (4);

  // Initial and final position (one dimensional spline).
  params[0] = 0.,  params[1] = 0.33;
  params[2] = 0.66,  params[3] = 1.;

  Spline spline (std::make_pair (0., 1.), 1, params);
  discreteInterval_t window (0., 1., 0.001);

  std::cerr << spline << std::endl;

  std::cout << "set term wxt persist font ',5'" << std::endl;
  std::cout << "set multiplot layout 2,3" << std::endl;

  for (unsigned i = 0; i < 3; ++i)
    {
      std::string title;
      switch (i)
	{
	case 0:
	  title = "spline";
	  break;
	case 1:
	  title = "spline derivative";
	  break;
	case 2:
	  title = "spline 2nd derivative";
	  break;
	default:
	  assert (0);
	}

      std::cout << "plot '-' title '"<< title <<"' with line\n" << std::endl;
      for (double t = boost::get<0> (window); t < boost::get<1> (window);
	   t += boost::get<2> (window))
	{
	  Spline::vector_t grad = spline.derivative (t, i);
	  //std::cerr << jac.size1 () << "/" << jac.size2 () << std::endl;
	  std::cout << (boost::format ("%1% %2%\n") % t % grad (0)).str ();
	}
      std::cout << "e" << std::endl;
    }

  for (unsigned i = 0; i < 3; ++i)
    {
      std::string title;
      switch (i)
	{
	case 0:
	  title = "configuration variation";
	  break;
	case 1:
	  title = "derivative variation";
	  break;
	case 2:
	  title = "2nd derivative variation";
	  break;
	default:
	  assert (0);
	}

      std::cout << "plot '-' title '"<< title <<"' with line\n" << std::endl;
      for (double t = boost::get<0> (window); t < boost::get<1> (window);
	   t += boost::get<2> (window))
	{
	  Spline::matrix_t jac = spline.variationDerivWrtParam (t, i);
	  //std::cerr << jac.size1 () << "/" << jac.size2 () << std::endl;
	  std::cout << (boost::format ("%1% %2%\n") % t % jac (0, 0)).str ();
	}
      std::cout << "e" << std::endl;
    }

  std::cout << "unset multiplot" << std::endl;
  return 0;
}

GENERATE_TEST ()
