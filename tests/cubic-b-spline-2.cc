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

#undef NDEBUG

#include <iostream>
#include <fstream>

#include <boost/format.hpp>

#include "shared-tests/fixture.hh"

#include <roboptim/core/decorator/finite-difference-gradient.hh>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

#include <roboptim/trajectory/visualization/trajectory.hh>

#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/cubic-b-spline.hh>

using namespace roboptim;
using namespace roboptim::trajectory;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;

typedef CubicBSpline::size_type size_type;
typedef CubicBSpline::value_type value_type;
typedef CubicBSpline::interval_t interval_t;

typedef Eigen::MatrixXd matrix_t;

// Limited data => use full precision serialization
typedef matrix_t store_t;

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (trajectory_cubic_b_spline)
{
  using namespace roboptim::visualization::gnuplot;
  size_type nbRows = 0;

  matrix_t reference = readMatrix<store_t> ("cubic-b-spline-2.dat");

  if (reference.size () == 0) {
    throw std::runtime_error
      ("failed to open file \"cubic-b-spline-2.dat\"");
  }
  double tol = 1e-6;

  size_type nbKnots = 12;

  CubicBSpline::knots_t knots (nbKnots);
  knots << -3, -2, -1, 0, .2, 1., 1.2, 1.8, 1.9, 2.2, 3.2, 4.2;

  std::vector<CubicBSpline::vector_t> params;

  for (size_type i = 0; i < nbKnots-4; ++i)
  {
    params.push_back (CubicBSpline::vector_t (nbKnots - 4));
    std::size_t i_ = static_cast<std::size_t> (i);
    params[i_].setZero ();
    params[i_][i] = 1;
  }

  value_type start = knots[3];
  value_type end = knots[nbKnots-4];
  value_type eval_step = 0.002;
  size_type n_eval = static_cast<size_type> (std::floor ((end - start)/eval_step)+1);
  size_type total_n_eval = (nbKnots - 4) * n_eval;
  matrix_t spline_data (total_n_eval, 2);
  spline_data.setZero ();

  for (std::size_t i = 0; i < (std::size_t)nbKnots - 4; ++i)
    {
      CubicBSpline spline (1, knots, params [i], "Cubic B-spline");

      discreteInterval_t interval (spline.timeRange ().first,
				   spline.timeRange ().second, eval_step);

      std::cout
	<<
	(boost::format
	 ("set term wxt persist title 'Spline: base functions' %1% font ',5'")
	 % i) << std::endl
	<< "set grid xtics noytics linewidth 0.5" << std::endl
	<< "set xlabel 't'" << std::endl
	<< "set ylabel 'value'" << std::endl;
      std::string title ("spline");
      std::cout << "plot '-' title '"<< title <<"' with line\n";

      spline (0.95);

      // Loop over the interval of definition
      for (double t = boost::get<0> (interval); t < boost::get<1> (interval);
	   t += boost::get<2> (interval))
	{
	  CubicBSpline::vector_t value = spline (t);
	  std::cout << (boost::format ("%1.3f %2.8f\n")
			% normalize (t)
			% normalize (value (0))).str ();
	  BOOST_CHECK_SMALL (t - reference (nbRows, 0), tol);
	  BOOST_CHECK_SMALL (value (0) - reference (nbRows, 1), tol);
	  spline_data.row (nbRows) << t, value (0);
	  ++nbRows;
	}

      std::cout << "e" << std::endl;
    }

  // Generate the archive
  //writeMatrix<matrix_t> ("cubic-b-spline-2.dat", spline_data);

  // Compare data
  BOOST_CHECK (allclose (spline_data, reference, tol));

  std::cout << "unset multiplot" << std::endl;
}
BOOST_AUTO_TEST_SUITE_END ()
