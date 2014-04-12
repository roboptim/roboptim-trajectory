/*
 * bspline_test.cc
 *
 *  Created on: Mar 13, 2013
 *      Author: wern_al
 */

# define private public
# define protected public
# include <roboptim/trajectory/cubic-b-spline.hh>
# include <roboptim/trajectory/b-spline.hh>
# include <roboptim/trajectory/constraint-b-spline.hh>
# undef private
# undef protected

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

# include <cstdlib>
# include <limits>
# include <cassert>
# include <cmath>
# include <iostream>

using namespace roboptim;
using namespace std;

#undef NDEBUG

typedef Function::vector_t vector_t;

/**
 * simple test if evaluating the spline thoughout the interval
 * works of start and end value of the spline are fixed through
 * constraints.
 */
template <int N>
void check_evaluate (std::pair<double, double> interval,
                     int dimension, vector_t const& params, int order)
{
  ConstraintBSpline<N> spline (interval, dimension, params);

  spline.addfixedConstraint (interval.first, 0, 0);
  spline.addfixedConstraint (interval.second, 0, 0);

  for (double t = interval.first; t < interval.second; t += 1e-3)
    {
      Eigen::Matrix<double, Eigen::Dynamic, 1> res (dimension);
      spline.derivative (res, t, order);
    }

}

template <int N>
void test_evaluate()
{
  std::pair<double, double> interval (0., 1.);
  const int params_no = 10;
  assert (params_no > N + 1);
  for (int dimension = 1; dimension < 4; dimension++)
    {
      Eigen::Matrix<double, Eigen::Dynamic, 1> params (dimension * params_no);
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

  std::pair<double, double> interval = std::make_pair (0., 1.);

  int min_params = N + 1;
  int params_c = N + 1 + N;
  assert (params_c >= min_params);
  Eigen::Matrix<double, Eigen::Dynamic, 1> params (params_c);
  params.setRandom();

  //BSpline<N> spline(interval,1,params,knots);
  ConstraintBSpline<N> spline (interval, 1, params);

  spline.addfixedConstraint (interval.first, 0, 0);
  spline.addfixedConstraint (interval.second, 0, 0);
  spline.addfixedConstraint (interval.second, 0, 1, 1);

  /*std::cout << "spline.constraint_values_" << std::endl << spline.constraint_values_ << std::endl;
    std::cout << "spline.constraints_" << std::endl << spline.constraints_ << std::endl;
    std::cout << "spline.projector_offset_" << std::endl << spline.projector_offset_ << std::endl;
    std::cout << "spline.projector_" << std::endl << spline.projector_ << std::endl;*/

  double delta;
  delta = std::abs (0. - spline.derivative (interval.first, 0) (0));
  assert ( delta < 1e-5 );

  delta = std::abs (0. - spline.derivative (interval.second, 0) (0));
  assert ( delta < 1e-5 );

  delta = std::abs (1. - spline.derivative (interval.second, 1) (0));
  assert ( delta < 1e-5 );

  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();
  discreteInterval_t plot_interval (interval.first, interval.second, 0.01);

  std::cout << "set terminal wxt persist" << std::endl;
  std::cout << "set multiplot" << std::endl;
  //		<< (gnuplot << plot (spline, plot_interval));
  std::cout << "plot '-' title 'B-Spline' with line" << std::endl;
  for (double t = interval.first; t < interval.second; t += 1e-2)
    {
      std::cout << t << " " << spline.derivative (t, 0) << std::endl;
    }
  std::cout << "e" << std::endl;
  std::cout << "plot '-' title 'B-Spline 1-deriv' with line" << std::endl;
  for (double t = interval.first; t < interval.second; t += 1e-2)
    {
      std::cout << t << " " << spline.derivative (t, 1) << std::endl;
    }
  std::cout << "e" << std::endl;
}

int main()
{
  srand (time (NULL));
  test_evaluate<3>();
  test_plot<3>();

}
