// Copyright (C) 2014 by Thomas Moulard, CNRS-AIST JRL UMI/3218.
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

#include "shared-tests/common.hh"

#include <roboptim/trajectory/cubic-b-spline.hh>
#include <roboptim/trajectory/b-spline.hh>
#include <roboptim/trajectory/spline-length.hh>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

#include <cstdlib>
#include <limits>
#include <cassert>
#include <cmath>
#include <iostream>

using namespace roboptim;
using namespace std;

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (trajectory_bspline)
{
  CubicBSpline::vector_t params (22);

  params <<
    // Initial position.
    0., 0.,
    // Control point .
    .1, -2.,
    // Control point .
    .2, -2.,
    // Control point .
    .3, -2.0,
    // Control point .
    .4, -2.0,
    // Control point .
    .5, -2.0,
    // Control point .
    .6, -2.0,
    // Control point .
    .7, -2.0,
    // Control point .
    .8, -2.0,
    // Control point .
    .9, -2.,
    // Final position.
    1., .1;

  CubicBSpline::interval_t timeRange = CubicBSpline::makeInterval (0., 4.);

  CubicBSpline spline (timeRange, 2, params, "spline");
  SplineLength splineLength (spline);

  std::cout << "length: " << splineLength (params) << std::endl;
}

BOOST_AUTO_TEST_SUITE_END ()
