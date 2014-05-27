// Copyright (C) 2013 by Alexander Werner, DLR.
// Copyright (C) 2014 by Benjamin Chrétien, CNRS-LIRMM.
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
#include <roboptim/trajectory/constrained-b-spline.hh>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

#include <cstdlib>
#include <limits>
#include <cmath>
#include <iostream>

using namespace roboptim;
using namespace std;

typedef Function::vector_t   vector_t;
typedef Function::argument_t argument_t;
typedef Function::interval_t interval_t;
typedef Function::value_type value_type;
typedef Function::size_type  size_type;

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

static double tol = 1e-6;

/// Simple test evaluating the spline throughout the interval.
/// Start and end values of the spline are fixed through constraints.
template <int N>
void check_evaluate (interval_t interval,
                     int dimension, const vector_t& params, int order)
{
  ConstrainedBSpline<N> spline (interval, dimension, params);

  // Initial parameters size
  size_type param_size = spline.parameters ().size ();

  // First point of the spline fixed at 0.
  spline.addFixedConstraint (interval.first, 0, 0);
  BOOST_CHECK (spline.parameters ().size () == param_size - 1);

  // Mid point of the spline fixed at 1.
  spline.addFixedConstraint (0.5*(interval.first+interval.second), 0, 1.);
  BOOST_CHECK (spline.parameters ().size () == param_size - 2);

  // Last point of the spline fixed at 0.
  spline.addFixedConstraint (interval.second, 0, 0);
  BOOST_CHECK (spline.parameters ().size () == param_size - 3);

  // Add couple constrained: 0.25 and 0.75 should have the same derivative.
  value_type range = interval.second - interval.first;
  spline.addCoupledConstraint (interval.first + 0.25 * range, 0,
                               interval.first + 0.75 * range, 0,
                               1, 1.);
  BOOST_CHECK (spline.parameters ().size () == param_size - 4);

  for (value_type t = interval.first; t < interval.second; t += 1e-3)
    {
      vector_t res (dimension);
      spline.derivative (res, t, order);
    }
}

template <int N>
void test_evaluate()
{
  interval_t interval (0., 1.);
  const int params_no = 10;
  BOOST_CHECK (params_no > N + 1);

  for (int dimension = 1; dimension < 4; dimension++)
    {
      argument_t params (dimension * params_no);
      params.setRandom ();

      for (int order = 0; order < N; order++)
        {
	  check_evaluate<N> (interval, dimension, params, order);
        }
    }
}


template <int N>
void test_plot (void)
{
  using namespace roboptim::visualization;
  using namespace roboptim::visualization::gnuplot;

  interval_t interval = std::make_pair (0., 1.);

  int min_params = N + 1;
  int params_c = N + 1 + N;
  BOOST_CHECK (params_c >= min_params);
  argument_t params (params_c);
  params.setRandom ();

  ConstrainedBSpline<N> spline (interval, 1, params);

  // First point of the spline fixed at 0.
  spline.addFixedConstraint (interval.first, 0, 0);

  // Mid point of the spline fixed at 1.
  spline.addFixedConstraint (0.5*(interval.first+interval.second), 0, 1.);

  // Last point of the spline fixed at 0.
  spline.addFixedConstraint (interval.second, 0, 0);

  // Add couple constrained: 0.25 and 0.75 should have the same derivative.
  value_type range = interval.second - interval.first;
  spline.addCoupledConstraint (interval.first + 0.25 * range, 0,
                               interval.first + 0.75 * range, 0,
                               1, 1.);
  value_type delta;
  delta = std::abs (0. - spline.derivative (interval.first, 0) (0));
  BOOST_CHECK_SMALL (delta, tol);

  delta = std::abs (0. - spline.derivative (interval.second, 0) (0));
  BOOST_CHECK_SMALL (delta, tol);

  delta = std::abs (spline.derivative (interval.first + 0.25 * range, 1) (0)
		    - spline.derivative (interval.first + 0.75 * range, 1) (0));
  BOOST_CHECK_SMALL (delta, tol);

  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();
  discreteInterval_t plot_interval (interval.first, interval.second, 0.01);

  // TODO: use RobOptim's gnuplot interface rather than stdout
  std::cout << "set terminal wxt persist" << std::endl;
  std::cout << "set multiplot layout 2,1" << std::endl;
  std::cout << "set xrange [" << interval.first << ":"
	    << interval.second << "]" << std::endl;
  std::cout << "set xtics 0.1" << std::endl;
  //		<< (gnuplot << plot (spline, plot_interval));
  std::cout << "plot '-' title 'B-Spline' with line" << std::endl;
  value_type delta_t = 1e-2;
  for (value_type t = interval.first; t < interval.second + delta_t;
       t += delta_t)
    {
      std::cout << t << " " << spline.derivative (t, 0)[0] << std::endl;
    }
  std::cout << "e" << std::endl;
  std::cout << "plot '-' title 'B-Spline 1-deriv' with line" << std::endl;
  for (value_type t = interval.first; t < interval.second + delta_t;
       t += delta_t)
    {
      std::cout << t << " " << spline.derivative (t, 1)[0] << std::endl;
    }
  std::cout << "e" << std::endl;
}

BOOST_AUTO_TEST_CASE (trajectory_constrained_bspline)
{
  srand (0);

  test_evaluate<4> ();
  test_plot<4> ();
}

BOOST_AUTO_TEST_SUITE_END ()
