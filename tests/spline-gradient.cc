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

#include <roboptim/trajectory/sys.hh>

#include <boost/format.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <roboptim/core/finite-difference-gradient.hh>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/cubic-b-spline.hh>
#include <roboptim/trajectory/trajectory-cost.hh>

#include "shared-tests/common.hh"

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;

int run_test ()
{
  CubicBSpline::vector_t params (7);

  for (unsigned x = 0; x < 7; ++x)
    {
      std::cout
	<<
	(boost::format
	 ("set term wxt persist title 'Spline: base functions' %1% font ',5'")
	 % x) << std::endl
	<< "set grid xtics noytics linewidth 0.5" << std::endl
	<< "set xlabel 'time'" << std::endl
	<< "set ylabel 'spline value'" << std::endl
	<< "set multiplot layout 5,3" << std::endl;



      params.clear ();
      params[x] = 1.;

      // Build a cubic spline of dimension 1
      CubicBSpline spline (std::make_pair (0., 4.), 1, params);
      discreteInterval_t window (0., 4., 0.01);
      std::cerr << spline << std::endl;

      // Loop over the order of derivation
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

	  std::cout << "plot '-' title '"<< title <<"' with line\n"
		    << std::endl;
	  // Loop over the interval of definition
	  for (double t = boost::get<0> (window); t < boost::get<1> (window);
	       t += boost::get<2> (window))
	    {
	      CubicBSpline::vector_t grad = spline.derivative (t, i);
	      std::cout << (boost::format ("%1.2f %2.2f\n")
			    % normalize (t)
			    % normalize (grad (0))).str ();

	      try
	      {
		// Check gradient
		Function::vector_t x (1);
		x[0] = t;
		if ((t < boost::get<1> (window)-boost::get<2> (window))
		    && (t>boost::get<0> (window))) {
		  checkGradientAndThrow (spline, 0, x);
		}
	      }
	      catch (BadGradient& bg)
		{
		  std::cerr << bg << std::endl;
		  std::cerr << "t=" << t << ", spline=" << spline << std::endl;
		}
	    }
	  std::cout << "e" << std::endl;
	}

      // Loop over control points
      for (unsigned j = 0; j < 7; ++j)
	// Loop over order of derivation
	for (unsigned i = 0; i < 3; ++i)
	  {
	    std::string title;
	    switch (i)
	      {
	      case 0:
		title = "configuration";
		break;
	      case 1:
		title = "derivative";
		break;
	      case 2:
		title = "2nd derivative";
		break;
	      default:
		assert (0);
	      }
	    title.append ((boost::format (" (%1%)") % j).str ());

	    std::cout << "plot '-' title '"<< title <<"' with line\n" << std::endl;
	  // Loop over the interval of definition
	    for (double t = boost::get<0> (window); t < boost::get<1> (window);
		 t += boost::get<2> (window))
	      {
		CubicBSpline::matrix_t jac =
		  spline.variationDerivWrtParam (t, i);
		std::cout << (boost::format ("%1.2f %2.2f\n")
			      % normalize (t)
			      % normalize (jac (0, j))).str ();
	      }
	std::cout << "e" << std::endl;
	  }
      std::cout << "unset multiplot" << std::endl;
    }

  return 0;
}

GENERATE_TEST ()
