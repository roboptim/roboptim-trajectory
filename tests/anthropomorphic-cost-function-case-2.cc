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

#include "shared-tests/common.hh"
#include "anthropomorphic-cost-function.hh"

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (trajectory_anthropomorphic_cost_function_case2)
{
  const double initialX = 0.;
  const double initialY = 0.;
  const double initialTheta = M_PI / 2.;

  const double finalX = 1.;
  const double finalY = 0.;
  const double finalTheta = 0.;

  int res = optimize (initialX,
		      initialY,
		      initialTheta,
		      finalX,
		      finalY,
		      finalTheta);

  BOOST_CHECK (res == 0);
}

BOOST_AUTO_TEST_SUITE_END ()
