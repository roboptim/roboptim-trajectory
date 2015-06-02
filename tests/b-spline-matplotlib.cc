// Copyright (C) 2015 by Félix Darricau, AIST, CNRS, EPITA.
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

#include <iostream>
#include <fstream>

#include <boost/format.hpp>

#include "shared-tests/fixture.hh"

#include <roboptim/core/io.hh>

#include <roboptim/trajectory/visualization/b-spline-matplotlib.hh>

#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/b-spline.hh>

#include <roboptim/core/visualization/matplotlib-function.hh>

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::matplotlib;
using namespace roboptim::trajectory;
using namespace roboptim::trajectory::visualization::matplotlib;

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (trajectory_b_spline_matplotlib)
{
typedef BSpline<3> spline_t;
typedef spline_t::value_type value_type;

  boost::shared_ptr<boost::test_tools::output_test_stream>
    output = retrievePattern ("b-spline-matplotlib");

  spline_t::vector_t params (2*11);
  params << 0  , 1. ,
    0.5, 8. ,
    1. , 10.,
    1.5, 8. ,
    2. , 5. ,
    2.5, 10.,
    3  , 0. ,
    3.5, 2. ,
    4. , 5. ,
    4.5, 10.,
    5  , 0. ;

  spline_t::vector_t params2 (2*11);
  params2 << 0  , 1. ,
    0.5, 8. ,
    1. , 10.,
    1.5, 8. ,
    2. , 5. ,
    2.5, 1. ,
    3  , 0. ,
    3.5, 2. ,
    4. , 5. ,
    4.5, 10.,
    5  , 0. ;

  spline_t::vector_t params3 (2*11);
  params3 << 0  , 1. ,
    0.5, 8. ,
    1. , 10.,
    1.5, 8. ,
    2. , 5. ,
    2.5, 15.,
    3  , 0. ,
    3.5, 2. ,
    4. , 5. ,
    4.5, 10.,
    5  , 0. ;

  boost::shared_ptr<spline_t> spline = boost::make_shared<spline_t>
    (std::make_pair (0., 1.), 2, params, "BSpline1", true);

  boost::shared_ptr<spline_t> spline2 = boost::make_shared<spline_t>
    (std::make_pair (0., 1.), 2, params2, "BSpline2", true);

  boost::shared_ptr<spline_t> spline3 = boost::make_shared<spline_t>
    (std::make_pair (0., 1.), 2, params3, "BSpline3", true);

  Matplotlib matplotlib = Matplotlib::make_matplotlib (std::make_pair(3, 1));
  value_type step = 0.005;

  (*output)
    << (matplotlib
        << plot_spline (*spline, step)
        << title ("BSpline, with P4 = {2.5, 10}")
        << plot_spline (*spline2, step)
        << title ("BSpline, with P4 = {2.5, 1}")
        << plot_spline (*spline3, step)
        << title ("BSpline, with P4 = {2.5, 15}")
	);


  std::cout << output->str() << std::endl;
  BOOST_CHECK (output->match_pattern ());
}

BOOST_AUTO_TEST_SUITE_END ()
