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

#include <roboptim/core/finite-difference-gradient.hh>

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

typedef std::vector < std::vector < value_type > > matrix_t;

matrix_t getReference (const std::string& file)
{
  matrix_t result;
  size_type nbRows = 0;

  std::ifstream f; f.open (file.c_str ());
  std::string line;

  while (f.good ()) {
    std::vector <value_type> data;
    std::getline (f, line);
    std::stringstream ss (line);
    std::copy (std::istream_iterator< value_type > (ss),
	       std::istream_iterator< value_type > (),
	       std::back_inserter (data));
    ++nbRows;
    result.push_back (data);
  }
  return result;
}

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (trajectory_cubic_b_spline)
{
  using namespace roboptim::visualization::gnuplot;
  size_type nbRows = 0;
  matrix_t reference = getReference ("cubic-b-spline-2.stdout");
  if (reference.size () == 0) {
    throw std::runtime_error
      ("failed to open file \"cubic-b-spline-2.stdout\"");
  }
  double tol = 1e-6;

  size_type nbKnots = 12;
  double knotArray [] =
    { -3, -2, -1, 0, .2, 1., 1.2, 1.8, 1.9, 2.2, 3.2, 4.2};
  CubicBSpline::knots_t knots (nbKnots);
  for (std::size_t i=0; i<(std::size_t)nbKnots ; ++i) knots [i] = knotArray [i];
  std::vector <CubicBSpline::vector_t> params;

  for (std::size_t i=0; i < (std::size_t)(nbKnots-4); ++i) {
    params.push_back (CubicBSpline::vector_t (nbKnots - 4));
    params [i].setZero ();
    params [i][i] = 1;
  }

  std::ofstream f; f.open ("output");
  for (std::size_t i=0; i<(std::size_t)nbKnots-4; ++i) {
    CubicBSpline spline (1, knots, params [i], "Cubic B-spline");

    discreteInterval_t interval (spline.timeRange ().first,
				 spline.timeRange ().second, 0.002);

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
	 t += boost::get<2> (interval)) {
      CubicBSpline::vector_t value = spline (t);
      std::cout << (boost::format ("%1.3f %2.8f\n")
		    % normalize (t)
		    % normalize (value (0))).str ();
      BOOST_CHECK_SMALL (t - reference [nbRows][0], tol);
      BOOST_CHECK_SMALL (value (0) - reference [nbRows][1], tol);
      ++nbRows;
    }
    std::cout << "e" << std::endl;
  }
  std::cout << "unset multiplot" << std::endl;
}
BOOST_AUTO_TEST_SUITE_END ()
