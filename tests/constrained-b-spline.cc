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

#include <roboptim/trajectory/cubic-b-spline.hh>
#include <roboptim/trajectory/b-spline.hh>
#include <roboptim/trajectory/constrained-b-spline.hh>

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

typedef Function::vector_t vector_t;
typedef Function::argument_t argument_t;
typedef Function::value_type value_type;

/// Simple test evaluating the spline throughout the interval.
/// Start and end value of the spline are fixed through constraints.
template <int N>
void check_evaluate (std::pair<value_type, value_type> interval,
                     int dimension, vector_t const& params, int order)
{
  ConstrainedBSpline<N> spline (interval, dimension, params);

  spline.addFixedConstraint (interval.first, 0, 0);
  spline.addFixedConstraint (interval.second, 0, 0);

  for (value_type t = interval.first; t < interval.second; t += 1e-3)
    {
      vector_t res (dimension);
      spline.derivative (res, t, order);
    }

}

template <int N>
void test_evaluate()
{
  std::pair<value_type, value_type> interval (0., 1.);
  const int params_no = 10;
  assert (params_no > N + 1);
  for (int dimension = 1; dimension < 4; dimension++)
    {
      argument_t params (dimension * params_no);
      params.setRandom();
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

  std::pair<value_type, value_type> interval = std::make_pair (0., 1.);

  int min_params = N + 1;
  int params_c = N + 1 + N;
  assert (params_c >= min_params);
  argument_t params (params_c);
  params.setRandom();

  //BSpline<N> spline(interval,1,params,knots);
  ConstrainedBSpline<N> spline (interval, 1, params);

  spline.addFixedConstraint (interval.first, 0, 0);
  spline.addFixedConstraint (interval.second, 0, 0);
  spline.addFixedConstraint (interval.second, 0, 1, 1);

  /*std::cout << "spline.constraint_values_" << std::endl << spline.constraint_values_ << std::endl;
    std::cout << "spline.constraints_" << std::endl << spline.constraints_ << std::endl;
    std::cout << "spline.projector_offset_" << std::endl << spline.projector_offset_ << std::endl;
    std::cout << "spline.projector_" << std::endl << spline.projector_ << std::endl;*/

  value_type delta;
  delta = std::abs (0. - spline.derivative (interval.first, 0) (0));
  assert ( delta < 1e-5 );

  delta = std::abs (0. - spline.derivative (interval.second, 0) (0));
  assert ( delta < 1e-5 );

  delta = std::abs (1. - spline.derivative (interval.second, 1) (0));
  assert ( delta < 1e-5 );

  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();
  discreteInterval_t plot_interval (interval.first, interval.second, 0.01);

  std::cout << "set terminal wxt persist" << std::endl;
  std::cout << "set multiplot layout 2,1" << std::endl;
  //		<< (gnuplot << plot (spline, plot_interval));
  std::cout << "plot '-' title 'B-Spline' with line" << std::endl;
  for (value_type t = interval.first; t < interval.second; t += 1e-2)
    {
      std::cout << t << " " << spline.derivative (t, 0)[0] << std::endl;
    }
  std::cout << "e" << std::endl;
  std::cout << "plot '-' title 'B-Spline 1-deriv' with line" << std::endl;
  for (value_type t = interval.first; t < interval.second; t += 1e-2)
    {
      std::cout << t << " " << spline.derivative (t, 1)[0] << std::endl;
    }
  std::cout << "e" << std::endl;
}

int main()
{
  srand (time (NULL));
  test_evaluate<3>();
  test_plot<3>();

}
