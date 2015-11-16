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

#include "shared-tests/fixture.hh"

#include <roboptim/trajectory/cubic-b-spline.hh>
#include <roboptim/trajectory/b-spline.hh>

#include <roboptim/core/decorator/finite-difference-gradient.hh>
#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

#include <cstdlib>
#include <limits>
#include <cassert>
#include <cmath>
#include <iostream>

#include <boost/test/floating_point_comparison.hpp>

using namespace roboptim;
using namespace roboptim::trajectory;

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

typedef Function::vector_t vector_t;
typedef DifferentiableFunction::gradient_t gradient_t;
typedef Function::value_type value_type;
typedef Function::size_type size_type;

// TODO: track floating-point errors.
static value_type tol = 1e8 * std::numeric_limits<value_type>::epsilon ();

template <int N>
struct SplineDerivWrtParameters : public DifferentiableFunction
{
  typedef BSpline<N> spline_t;

  SplineDerivWrtParameters (const spline_t& spline, value_type t)
    : DifferentiableFunction
      (spline.parameters ().size (), spline.outputSize (),
       "spline differentiable w.r.t parameters"),
      spline_ (spline),
      t_ (t)
  {}

  virtual void
  impl_compute (result_ref result, const_argument_ref x)
    const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    bool cur_malloc_allowed = is_malloc_allowed ();
    set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    spline_t spline (spline_);
    spline.setParameters (x);
    result = spline (t_);

#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    set_is_malloc_allowed (cur_malloc_allowed);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION
  }

  virtual void
  impl_gradient (gradient_ref gradient,
		 const_argument_ref x,
		 size_type functionId = 0)
    const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    bool cur_malloc_allowed = is_malloc_allowed ();
    set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    spline_t spline (spline_);
    spline.setParameters (x);
    gradient = spline.variationConfigWrtParam (t_).row (functionId);

#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    set_is_malloc_allowed (cur_malloc_allowed);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION
  }

private:
  const spline_t& spline_;
  value_type t_;
};

template <int N>
struct spline_checks
{
  static void
  check_evaluate (const std::pair<value_type, value_type>&,
		  int, vector_t const&, int, bool);

  static void
  check_fd (const std::pair<value_type, value_type>&,
	    vector_t const&, bool);

  template <int ORDER>
  static void
  check_derivative (const std::pair<value_type, value_type>&,
		    int, vector_t const&);

  static void check_non_uniform
  (const std::pair<value_type, value_type>& interval, int dimension,
   vector_t const& params, vector_t const& knots, int order);
};


// Evaluate the spline and its derivatives, thus implicitly checking
// basisPolynomials_ and interval() against old CubicBSpline
// implementation.
template <>
void
spline_checks<3>::check_evaluate
(const std::pair<value_type, value_type>& interval, int dimension,
 vector_t const& params, int order, bool clamped)
{
  CubicBSpline old_spline (interval, dimension, params, "", clamped);
  BSpline<3> new_spline (interval, dimension, params, "", clamped);

  gradient_t old_res (dimension);
  gradient_t new_res (dimension);
  gradient_t delta (dimension);

  for (value_type t = interval.first; t < interval.second; t += 1e-3)
    {
      old_spline.derivative (old_res, t, order);
      new_spline.derivative (new_res, t, order);
      delta = old_res - new_res;

      value_type max_err = delta.lpNorm<Eigen::Infinity> ();
      BOOST_CHECK_SMALL (max_err, tol);
    }
}

template <>
void spline_checks<5>::check_evaluate
(const std::pair<value_type, value_type>& interval, int dimension,
 vector_t const& params, int order, bool clamped)
{
  BSpline<5> new_spline (interval, dimension, params, "", clamped);

  for (value_type t = interval.first; t < interval.second; t += 1e-3)
    {
      gradient_t new_res (dimension);
      new_spline.derivative (new_res, t, order);
      value_type min = params.minCoeff ();
      value_type max = params.maxCoeff ();
      if (order == 0)
        {
	  // test if spline is in between minimum and maximum
	  // parameter
	  BOOST_CHECK ((new_res.array () > min).all ()
		       && (new_res.array () < max).all ());
	  // FIXME: very generous assumption - especially for
	  // multi-dimensional splines
        }
    }
}

template <int N>
void
spline_checks<N>::check_fd
(const std::pair<value_type, value_type>& interval,
 vector_t const& params, bool clamped)
{
  BSpline<N> spline (interval, 1, params, "", clamped);

  // Check gradients with finite-differences
  for (value_type t = 0.; t < 1.; t += 1e-3)
    {
      try
        {
          vector_t x (1);
          x[0] = t;
          checkGradientAndThrow (spline, 0, x);

          SplineDerivWrtParameters<N> splineDerivWrtParams (spline, t);
          checkGradientAndThrow (splineDerivWrtParams, 0, spline.parameters ());
        }
      catch (BadGradient<EigenMatrixDense>& bg)
        {
          std::cerr << bg << std::endl;
          BOOST_CHECK(false);
        }
    }
}

template <int N>
template <int ORDER>
void
spline_checks<N>::check_derivative
(const std::pair<value_type, value_type>& interval, int dimension,
 vector_t const& params)
{
  BSpline<N> spline (interval, dimension, params);
  BSpline<N> deriv = spline.template derivative<ORDER> ();
  gradient_t delta;

  for (value_type t = interval.first; t < interval.second; t += 1e-3)
    {
      delta = deriv (t) - spline.derivative (t, ORDER);

      value_type max_err = delta.lpNorm<Eigen::Infinity> ();
      BOOST_CHECK_SMALL (max_err, tol);
    }
}

template <int N>
void
spline_checks<N>::check_non_uniform
(const std::pair<value_type, value_type>& interval, int dimension,
 vector_t const& params, vector_t const& knots, int order)
{
  BSpline<N> new_spline (interval, dimension, params, knots);

  for (value_type t = interval.first; t < interval.second; t += 1e-3)
    {
      gradient_t res (dimension);
      new_spline.derivative (res, t, order);
      value_type min = params.minCoeff ();
      value_type max = params.maxCoeff ();
      if (order == 0)
        {
	  // test if spline is in between minimum and maximum parameter
	  BOOST_CHECK (res.minCoeff() >= min);
	  BOOST_CHECK (res.maxCoeff() <= max);
	  // FIXME: very generous assumption - especially for
	  // multi-dimensional splines
        }
    }
}


template <int N>
void test_instantiate ()
{
  vector_t params (10);
  params.setZero ();

  typename BSpline<N>::interval_t interval = std::make_pair (0., 1.);
  {
    BSpline<N> spline (interval, 1, params);

    const vector_t& kv = spline.knotVector ();

    // Check intervals
    for (size_t i = 0; i <= 10; ++i)
      {
	value_type t = static_cast<value_type> (i) * 0.099;
	size_type k = static_cast<size_type> (spline.interval (t));
	BOOST_CHECK (kv[k] <= t);
	BOOST_CHECK (t <= kv[k+1]);
      }
  }

  vector_t knots (params.size () + N + 1);
  knots.setLinSpaced (0., 1.);
  {
    BSpline<N> spline (interval, 1, params, knots);
  }
}

template <int N>
void test_evaluate ()
{
  for (int derivative = 0; derivative <= N; derivative++)
    {
      for (int dimension = 1; dimension < 2; dimension++)
        {
	  std::pair<value_type, value_type> interval = std::make_pair (0., 1.);
	  int min_params = N + 1 + 10;
	  vector_t params (dimension * min_params);
	  params.setRandom ();
	  params *= 10.;

	  spline_checks<N>::check_evaluate
            (interval, dimension, params, derivative, false);

	  spline_checks<N>::check_evaluate
            (interval, dimension, params, derivative, true);
        }
    }
}

template <int N, int ORDER>
void test_derivative ()
{
  for (int dimension = 1; dimension < 2; dimension++)
    {
      typename BSpline<N>::interval_t interval = std::make_pair (0., 1.);
      int min_params = N + 1 + 10;
      vector_t params (dimension * min_params);
      params.setRandom ();
      params *= 10.;

      spline_checks<N>::template check_derivative<ORDER> (interval, dimension, params);
    }
}

template <int N>
void test_fd ()
{
  int min_params = N + 1 + 10;
  vector_t params (min_params);
  params.setRandom ();
  params *= 10.;

  typename BSpline<N>::interval_t interval = std::make_pair (0., 1.);

  spline_checks<N>::check_fd (interval, params, false);
  spline_checks<N>::check_fd (interval, params, true);
}

template <int N>
void test_non_uniform ()
{
  for (int derivative = 0; derivative <= N; derivative++)
    {
      for (int dimension = 1; dimension < 4; dimension++)
        {
	  std::pair<value_type, value_type> interval = std::make_pair (0., 1.);
	  int params_c = N + 1 + N;
	  vector_t params (dimension * params_c);
	  params.setRandom ();
	  params *= 10.;

	  vector_t knots (params_c + N + 1);
	  knots.head (N).setConstant (interval.first);
	  knots.segment (N, 2 + N).setLinSpaced (interval.first, interval.second);
	  knots.tail (N).setConstant (interval.second);

	  spline_checks<N>::check_non_uniform
            (interval, dimension, params, knots, derivative);
        }
    }
}


template <int N>
void test_plot ()
{
  using namespace roboptim::visualization;
  using namespace roboptim::visualization::gnuplot;

  // Fix clash with std::set
  using roboptim::visualization::gnuplot::set;

  std::stringstream filename;
  filename << "b-spline-" << N;
  boost::shared_ptr<boost::test_tools::output_test_stream>
    output = retrievePattern (filename.str ());

  std::pair<value_type, value_type> interval = std::make_pair (0., 5.);

  int params_c = N + 1 + N;
  vector_t params (params_c);

  for (size_type i = 0; i < params_c; ++i)
    {
      params (i) = std::pow (-1., i) * static_cast<value_type> (i);
    }

  std::stringstream name;
  name << "B-spline (order " << N << ")";
  BSpline<N> spline (interval, 1, params, name.str ());
  BSpline<N> deriv = spline.template derivative<1> ();
  name.str ("");

  name << "Clamped B-spline (order " << N << ")";
  BSpline<N> clamped_spline (interval, 1, params, name.str (), true);
  BSpline<N> clamped_deriv = clamped_spline.template derivative<1> ();

  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();
  discreteInterval_t plot_interval (interval.first, interval.second, 0.01);

  // Use RobOptim's gnuplot functions
  (*output) << (gnuplot
		<< set ("multiplot layout 2,2")

		<< comment (spline)

		<< comment ("Values:")
		<< comment (spline (0.))
		<< comment (spline (2.5))
		<< comment (spline (5.))

		<< comment ("1st derivative:")
		<< comment (spline.derivative (0., 1))
		<< comment (spline.derivative (2.5, 1))
		<< comment (spline.derivative (5., 1))

		<< comment ("2nd derivative:")
		<< comment (spline.derivative (0., 2))
		<< comment (spline.derivative (2.5, 2))
		<< comment (spline.derivative (5., 2))
		<< plot (spline, plot_interval)
		<< plot (deriv, plot_interval)

		<< comment (clamped_spline)

		<< comment ("Values:")
		<< comment (clamped_spline (0.))
		<< comment (clamped_spline (2.5))
		<< comment (clamped_spline (5.))

		<< comment ("1st derivative:")
		<< comment (clamped_spline.derivative (0., 1))
		<< comment (clamped_spline.derivative (2.5, 1))
		<< comment (clamped_spline.derivative (5., 1))

		<< comment ("2nd derivative:")
		<< comment (clamped_spline.derivative (0., 2))
		<< comment (clamped_spline.derivative (2.5, 2))
		<< comment (clamped_spline.derivative (5., 2))
		<< plot (clamped_spline, plot_interval)
		<< plot (clamped_deriv, plot_interval));

  std::cout << output->str () << std::endl;

  BOOST_CHECK (output->match_pattern ());
}


BOOST_AUTO_TEST_CASE (trajectory_bspline)
{
  test_instantiate<2> ();
  test_instantiate<3> ();
  test_instantiate<4> ();
  test_instantiate<5> ();

  srand (static_cast<unsigned int> (time (0)));
  test_plot<3> ();
  test_plot<5> ();

  test_evaluate<3> ();
  test_evaluate<5> ();

  test_derivative<3,1> ();
  test_derivative<3,2> ();
  test_derivative<5,1> ();
  test_derivative<5,3> ();

  test_fd<3> ();
  test_fd<5> ();

  test_non_uniform<3> ();
  test_non_uniform<4> ();
  test_non_uniform<5> ();
}

BOOST_AUTO_TEST_SUITE_END ()
