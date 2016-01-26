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

#include <roboptim/core/util.hh>

#include <roboptim/trajectory/polynomial.hh>
#include <roboptim/trajectory/polynomial-3.hh>

#include <limits>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/int.hpp>
#include <unsupported/Eigen/Polynomials>

using namespace roboptim;
using namespace roboptim::trajectory;
using namespace std;

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

static double tol = 1e-6;

template <int N>
void test_multiply ()
{
  typename roboptim::trajectory::Polynomial<N>::coefs_t params;
  params.setRandom ();
  double t0 = (double)rand () / RAND_MAX;
  double t1 = (double)rand () / RAND_MAX;

  roboptim::trajectory::Polynomial<N> p_0 (t0, params);
  p_0 (0.1);

  params.setZero ();
  params[0] = 0.3;
  roboptim::trajectory::Polynomial<N> p_1 (t1, params);
  p_1 (0.1);

  // Polynomial multiplication
  roboptim::trajectory::Polynomial<2*N> p_2 = p_0 * p_1;
  p_2 (0.1);
  // Polynomial multiplication + crop
  roboptim::trajectory::Polynomial<N> p_3 = p_0 * p_1;
  // Compare uncropped part
  BOOST_CHECK (allclose (p_2.coefs ().template head<N+1> (),
                         p_3.coefs ()));

  // Scalar multiplication
  roboptim::trajectory::Polynomial<N> q_0 = 2.0 * p_0;
  roboptim::trajectory::Polynomial<N> q_1 = p_0 * 2.0;

  BOOST_CHECK (allclose (q_0.coefs (), q_1.coefs ()));
}


template <int N>
struct poly_props
{
  static void check_derivative
    (const roboptim::trajectory::Polynomial<N>&, double, int, double)
  {
    // Not implemented for this N
    assert (0);
  }

  static void check_translate (const roboptim::trajectory::Polynomial<N>&, double,
			       const roboptim::trajectory::Polynomial<N>&)
  {
    // Not implemented for this N
    assert (0);
  }

  static void check_evaluate (double t0,
			      typename roboptim::trajectory::Polynomial<N>::coefs_t& coefs,
			      double t)
  {
    roboptim::trajectory::Polynomial<N> poly (t0, coefs);

    // Use Eigen's evaluator: translate polynomial from (t-t₀) to t
    roboptim::trajectory::Polynomial<N> eigen_poly = poly.translate (0.);
    double res1 = Eigen::poly_eval (eigen_poly.coefs (), t);
    double res2 = poly (t);

# if (defined ROBOPTIM_HAS_FENV_H && defined ENABLE_SIGFPE)
    // FPE may happen in Boost checks, in release
    roboptim::detail::DisableFPE d;
# endif //! (defined ROBOPTIM_HAS_FENV_H && defined ENABLE_SIGFPE)

    BOOST_CHECK_CLOSE (res1, res2, tol);
  }
};


// Code from Polynomial3
template <>
void poly_props<3>::check_derivative
(const roboptim::trajectory::Polynomial<3>& poly,
 double t, int order,
 double derivative)
{
  double dt = t - poly.t0 ();
  double dt2 = dt * dt;
  double a1 = poly.coefs () [1];
  double a2 = poly.coefs () [2];
  double a3 = poly.coefs () [3];

  double result;
  switch (order)
    {
    case 0:
      result = poly (t);
      break;
    case 1:
      result = a1 + 2 * a2 * dt + 3 * a3 * dt2;
      break;
    case 2:
      result = 2 * a2 + 6 * a3 * dt;
      break;
    case 3:
      result = 6 * a3;
      break;
    default:
      result = 0;
      break;
    }

  BOOST_CHECK_CLOSE (result, derivative, tol);
}

template <>
void poly_props<5>::check_derivative
(const roboptim::trajectory::Polynomial<5>& poly,
 double t, int order, double derivative)
{
  double dt = t - poly.t0 ();
  double dt2 = dt * dt;
  double dt3 = dt * dt2;
  double dt4 = dt * dt3;
  double a1 = poly.coefs () [1];
  double a2 = poly.coefs () [2];
  double a3 = poly.coefs () [3];
  double a4 = poly.coefs () [4];
  double a5 = poly.coefs () [5];

  double result;
  switch (order)
    {
    case 0:
      result = poly (t);
      break;
    case 1:
      result = a1 + 2 * a2 * dt + 3 * a3 * dt2 + 4 * a4 * dt3 + 5 * a5 * dt4;
      break;
    case 2:
      result = 2 * a2 + 6 * a3 * dt + 12 * a4 * dt2 + 20 * a5 * dt3;
      break;
    case 3:
      result = 6 * a3 + 24 * a4 * dt + 60 * a5 * dt2;
      break;
    case 4:
      result = 24 * a4 + 120 * a5 * dt;
      break;
    case 5:
      result = 120 * a5;
      break;
    default:
      result = 0;
      break;
    }

  BOOST_CHECK_CLOSE (result, derivative, tol);
}


template <>
void poly_props<3>::check_translate
(const roboptim::trajectory::Polynomial<3>& poly,
 double t1,
 const roboptim::trajectory::Polynomial<3>& ref_poly)
{
  double dt = t1 - poly.t0 ();
  double dt2 = dt * dt;
  double dt3 = dt * dt2;
  double a0 = poly.coefs () [0];
  double a1 = poly.coefs () [1];
  double a2 = poly.coefs () [2];
  double a3 = poly.coefs () [3];

  roboptim::trajectory::Polynomial<3>::coefs_t temp;
  temp [0] = a0 + a1 * dt + a2 * dt2 + a3 * dt3;
  temp [1] = a1 + 2 * dt * a2 + 3 * dt2 * a3;
  temp [2] = a2 + 3 * dt * a3;
  temp [3] = a3;
  roboptim::trajectory::Polynomial<3> p_new (t1, temp);

# if (defined ROBOPTIM_HAS_FENV_H && defined ENABLE_SIGFPE)
  // FPE may happen in Boost checks, in release
  roboptim::detail::DisableFPE d;
# endif //! (defined ROBOPTIM_HAS_FENV_H && defined ENABLE_SIGFPE)

  BOOST_CHECK_CLOSE (p_new.t0 (),
                     ref_poly.t0 (),
                     tol);

  for (int idx = 0; idx < 3 + 1; idx++)
    {
      BOOST_CHECK_CLOSE (p_new.coefs ()[idx],
                         ref_poly.coefs ()[idx],
                         tol);
    }
}


template <>
void poly_props<5>::check_translate
(const roboptim::trajectory::Polynomial<5>& poly,
 double t1,
 const roboptim::trajectory::Polynomial<5>& ref_poly)
{
  double dt = t1 - poly.t0 ();
  double dt2 = dt * dt;
  double dt3 = dt * dt2;
  double dt4 = dt * dt3;
  double dt5 = dt * dt4;
  double a0 = poly.coefs ()[0];
  double a1 = poly.coefs ()[1];
  double a2 = poly.coefs ()[2];
  double a3 = poly.coefs ()[3];
  double a4 = poly.coefs ()[4];
  double a5 = poly.coefs ()[5];


  roboptim::trajectory::Polynomial<5>::coefs_t temp;
  temp [0] = a0 +     a1 * dt +      a2 * dt2 +      a3 * dt3 +     a4 * dt4 + a5 * dt5;
  temp [1] = a1 + 2 * a2 * dt +  3 * a3 * dt2 +  4 * a4 * dt3 + 5 * a5 * dt4;
  temp [2] = a2 + 3 * a3 * dt +  6 * a4 * dt2 + 10 * a5 * dt3;
  temp [3] = a3 + 4 * a4 * dt + 10 * a5 * dt2;
  temp [4] = a4 + 5 * a5 * dt;
  temp [5] = a5;
  roboptim::trajectory::Polynomial<5> p_new (t1, temp);

# if (defined ROBOPTIM_HAS_FENV_H && defined ENABLE_SIGFPE)
  // FPE may happen in Boost checks, in release
  roboptim::detail::DisableFPE d;
# endif //! (defined ROBOPTIM_HAS_FENV_H && defined ENABLE_SIGFPE)

  BOOST_CHECK_CLOSE (p_new.t0 (), ref_poly.t0 (), tol);

  for (int idx = 0; idx < 5 + 1; idx++)
    {
      BOOST_CHECK_CLOSE (p_new.coefs ()[idx],
                         ref_poly.coefs ()[idx],
                         tol);
    }
}

template <int N>
struct test_derivative_loop
{
  test_derivative_loop (const roboptim::trajectory::Polynomial<N>& p,
                        Function::value_type t)
    : p_ (p),
      t_ (t)
  {
  }

  template <typename U>
  void operator () (U)
  {
    const int order = U::value;
    double derivative = p_.derivative (t_, order);
    double val = p_.template derivative<order> () (t_);

    {
# if (defined ROBOPTIM_HAS_FENV_H && defined ENABLE_SIGFPE)
      // FPE may happen in Boost checks, in release
      roboptim::detail::DisableFPE d;
# endif //! (defined ROBOPTIM_HAS_FENV_H && defined ENABLE_SIGFPE)

      BOOST_CHECK_CLOSE (derivative, val, tol);

      if (order > N)
        BOOST_CHECK_SMALL (derivative, tol);
    }

    poly_props<N>::check_derivative (p_, t_, order, derivative);
  }

  const roboptim::trajectory::Polynomial<N>& p_;
  Function::value_type t_;
};

template <int N>
void test_derivative ()
{
  typename roboptim::trajectory::Polynomial<N>::coefs_t params;
  params.setRandom ();
  double t0 = (double)rand () / RAND_MAX;
  roboptim::trajectory::Polynomial<N> p_1 (t0, params);

  double t = (double)rand () / RAND_MAX;

  boost::mpl::for_each<boost::mpl::range_c<int,0,N+1> >
    (test_derivative_loop<N> (p_1, t));
}

template <int N>
void test_translate ()
{
  typename roboptim::trajectory::Polynomial<N>::coefs_t params;
  params.setRandom ();
  double t0 = (double)rand () / RAND_MAX;
  double t1 = (double)rand () / RAND_MAX;
  roboptim::trajectory::Polynomial<N> p_1 (t0, params);
  roboptim::trajectory::Polynomial<N> p_2 = p_1.translate (t1);

  {
# if (defined ROBOPTIM_HAS_FENV_H && defined ENABLE_SIGFPE)
    // FPE may happen in Boost checks, in release
    roboptim::detail::DisableFPE d;
# endif //! (defined ROBOPTIM_HAS_FENV_H && defined ENABLE_SIGFPE)
    BOOST_CHECK_CLOSE (p_1 (0.), p_2 (0.), tol);
  }

  poly_props<N>::check_translate (p_1, t1, p_2);

  // Test translate in place.
  roboptim::trajectory::Polynomial<N> p_3 = p_1;
  p_3.translateInPlace (t1);
  BOOST_CHECK_CLOSE (p_2.t0 (), p_3.t0 (), tol);
  BOOST_CHECK (allclose (p_2.coefs (), p_3.coefs ()));
}

template <int N>
void test_copy ()
{
  typename roboptim::trajectory::Polynomial<N>::coefs_t params;
  params.setRandom ();

  // Same degree
  roboptim::trajectory::Polynomial<N> p_0 (3., params);
  roboptim::trajectory::Polynomial<N> p_1 = p_0;

  BOOST_CHECK_CLOSE (p_0.t0 (), p_1.t0 (), tol);
  BOOST_CHECK (allclose (p_0.coefs (), p_1.coefs ()));

  // Different degrees
  roboptim::trajectory::Polynomial<N+2> p_2 = p_0;
  BOOST_CHECK_CLOSE (p_0.t0 (), p_2.t0 (), tol);
  BOOST_CHECK (allclose (p_0.coefs (),
                         p_2.coefs ().template head<N+1> ()));

  roboptim::trajectory::Polynomial<N-1> p_3 = p_0;
  BOOST_CHECK_CLOSE (p_0.t0 (), p_3.t0 (), tol);
  BOOST_CHECK (allclose (p_3.coefs (),
                         p_0.coefs ().template head<N> ()));
}

template <int N>
void test_evaluate ()
{
  typename roboptim::trajectory::Polynomial<N>::coefs_t params;
  params.setRandom ();
  double t0 = (double)rand () / RAND_MAX;
  double t = (double)rand () / RAND_MAX;
  poly_props<N>::check_evaluate (t0, params, t);
}

template <int N>
void test_roots ()
{
  typedef typename roboptim::trajectory::Polynomial<N>::value_type value_type;

  typename roboptim::trajectory::Polynomial<N>::coefs_t params;
  params.setRandom ();
  double t0 = (double)rand () / RAND_MAX;
  roboptim::trajectory::Polynomial<N> p (t0, params);

  std::vector<value_type> roots = p.realRoots ();
  for (typename std::vector<value_type>::const_iterator
	 iter = roots.begin ();
       iter != roots.end ();
       ++iter)
    BOOST_CHECK_SMALL (p (*iter), tol);

  // Test other cases: null leading coefficient or null polynomial
  roots.clear ();
  p.coefs ()[N] = 0.;
  roots = p.realRoots ();
  for (typename std::vector<value_type>::const_iterator
	 iter = roots.begin ();
       iter != roots.end ();
       ++iter)
    BOOST_CHECK_SMALL (p (*iter), tol);

  roots.clear ();
  p.coefs ().setZero ();
  BOOST_CHECK_THROW (roots = p.realRoots (), std::runtime_error);

  // Test constant polynomial
  p.coefs ()[0] = 12.;
  BOOST_CHECK_THROW (roots = p.realRoots (), std::runtime_error);

  // Test linear polynomial
  p.coefs ()[1] = 3.;
  BOOST_CHECK_NO_THROW (roots = p.realRoots ());
  BOOST_CHECK (roots.size () == 1);
  BOOST_CHECK_CLOSE (roots[0], t0-4, 1e-12);
}

template <int N>
void test_min ()
{
  typedef typename roboptim::trajectory::Polynomial<N>::value_type value_type;
  typedef typename roboptim::trajectory::Polynomial<N>::min_t      min_t;
  typedef typename roboptim::trajectory::Polynomial<N>::interval_t interval_t;
  typedef typename roboptim::trajectory::Polynomial<N>::coefs_t    coefs_t;

  // Start with some dummy examples
  coefs_t params;
  interval_t interval;
  value_type t0 = 0.;
  min_t res_min;

  // Linear: 1+x on [-1,1]
  params.setZero ();
  params[0] = 1.;
  params[1] = 1.;
  interval.first  = -1.;
  interval.second =  1.;
  roboptim::trajectory::Polynomial<N> linear (t0, params);
  res_min = linear.min (interval);
  BOOST_CHECK_CLOSE (res_min.first, -1, tol);
  BOOST_CHECK_SMALL (res_min.second, tol);

  // Quadratic: 1+x² on [-1,1]
  params.setZero ();
  params[0] = 1.;
  params[2] = 1.;
  interval.first  = -1.;
  interval.second =  1.;
  roboptim::trajectory::Polynomial<N> quadratic (t0, params);
  res_min = quadratic.min (interval);
  BOOST_CHECK_SMALL (res_min.first, tol);

  // Cubic: 1+x³ on [0,1]
  params.setZero ();
  params[0] = 1.;
  params[3] = 1.;
  interval.first  = 0.;
  interval.second = 1.;
  roboptim::trajectory::Polynomial<N> cubic (t0, params);
  res_min = cubic.min (interval);
  BOOST_CHECK_SMALL (res_min.first, tol);

  // Test with a random polynomial
  params.setRandom ();
  t0 = (value_type)rand () / RAND_MAX;
  roboptim::trajectory::Polynomial<N> p (t0, params);
  roboptim::trajectory::Polynomial<N-1> dp = p.template derivative<1> ();

  interval.first  = -10. * std::abs ((double)rand () / RAND_MAX);
  interval.second =  10. * std::abs ((double)rand () / RAND_MAX);

  res_min = p.min (interval);
  // TODO: use a better check (not just dP(t_min) = 0)
  if (res_min.first != interval.first && res_min.first != interval.second)
    BOOST_CHECK_SMALL (dp (res_min.first), tol);

  // Test other cases: null leading coefficient or null polynomial
  p.coefs ()[N] = 0.;
  dp = p.template derivative<1> ();
  res_min = p.min (interval);
  if (res_min.first != interval.first && res_min.first != interval.second)
    BOOST_CHECK_SMALL (dp (res_min.first), tol);

  p.coefs ().setZero ();
  BOOST_CHECK_THROW (res_min = p.min (interval, false), std::runtime_error);
}

template <int N>
void test_max ()
{
  typedef typename roboptim::trajectory::Polynomial<N>::value_type value_type;
  typedef typename roboptim::trajectory::Polynomial<N>::max_t      max_t;
  typedef typename roboptim::trajectory::Polynomial<N>::interval_t interval_t;
  typedef typename roboptim::trajectory::Polynomial<N>::coefs_t    coefs_t;

  // Start with some dummy examples
  coefs_t params;
  interval_t interval;
  value_type t0 = 0.;
  max_t res_max;

  // Linear: 1+x on [-1,1]
  params.setZero ();
  params[0] = 1.;
  params[1] = 1.;
  interval.first  = -1.;
  interval.second =  1.;
  roboptim::trajectory::Polynomial<N> linear (t0, params);
  res_max = linear.max (interval);
  BOOST_CHECK_CLOSE (res_max.first, 1, tol);
  BOOST_CHECK_SMALL (res_max.second, 2 + tol);

  // Quadratic: 1+x² on [-1,1]
  params.setZero ();
  params[0] = 1.;
  params[2] = 1.;
  interval.first  = -1.;
  interval.second =  1.;
  roboptim::trajectory::Polynomial<N> quadratic (t0, params);
  res_max = quadratic.max (interval);
  BOOST_CHECK_CLOSE (res_max.second, 2, tol);

  // Cubic: 1+x³ on [0,1]
  params.setZero ();
  params[0] = 1.;
  params[3] = 1.;
  interval.first  = 0.;
  interval.second = 1.;
  roboptim::trajectory::Polynomial<N> cubic (t0, params);
  res_max = cubic.max (interval);
  BOOST_CHECK_CLOSE (res_max.first, 1, tol);
  BOOST_CHECK_CLOSE (res_max.second, 2, tol);

  // Test with a random polynomial
  params.setRandom ();
  t0 = (value_type)rand () / RAND_MAX;
  roboptim::trajectory::Polynomial<N> p (t0, params);
  roboptim::trajectory::Polynomial<N-1> dp = p.template derivative<1> ();

  interval.first  = -10. * std::abs ((double)rand () / RAND_MAX);
  interval.second =  10. * std::abs ((double)rand () / RAND_MAX);

  res_max = p.max (interval);
  // TODO: use a better check (not just dP(t_max) = 0)
  if (res_max.first != interval.first && res_max.first != interval.second)
    BOOST_CHECK_SMALL (dp (res_max.first), tol);

  // Test other cases: null leading coefficient or null polynomial
  p.coefs ()[N] = 0.;
  dp = p.template derivative<1> ();
  res_max = p.max (interval);
  if (res_max.first != interval.first && res_max.first != interval.second)
    BOOST_CHECK_SMALL (dp (res_max.first), tol);

  p.coefs ().setZero ();
  BOOST_CHECK_THROW (res_max = p.max (interval, false), std::runtime_error);
}

template <int N>
void test_misc ()
{
  typename roboptim::trajectory::Polynomial<N>::coefs_t params;
  params.setZero ();
  roboptim::trajectory::Polynomial<N> p (1., params);

  BOOST_CHECK (p.isNull ());
  BOOST_CHECK (p.isConstant ());
  BOOST_CHECK (p.isLinear ());
  BOOST_CHECK (p.trueOrder () == 0);

  p.coefs ()[0] = 2.;
  BOOST_CHECK (!p.isNull ());
  BOOST_CHECK (p.isConstant ());
  BOOST_CHECK (p.isLinear ());
  BOOST_CHECK (p.trueOrder () == 0);

  p.coefs ()[1] = 2.;
  BOOST_CHECK (!p.isNull ());
  BOOST_CHECK (!p.isConstant ());
  BOOST_CHECK (p.isLinear ());
  BOOST_CHECK (p.trueOrder () == 1);

  p.coefs ()[2] = 2.;
  BOOST_CHECK (!p.isNull ());
  BOOST_CHECK (!p.isConstant ());
  BOOST_CHECK (!p.isLinear ());
  BOOST_CHECK (p.trueOrder () == 2);

  typedef typename roboptim::trajectory::Polynomial<N>::polynomialFunction_t
    polynomialFunction_t;
  polynomialFunction_t func = p.asFunction ();
  typename polynomialFunction_t::argument_t x (func.inputSize ());
  x[0] = 0.4;
  BOOST_CHECK (func (x)[0] == p (x[0]));
}

void test_print ()
{
  boost::shared_ptr<boost::test_tools::output_test_stream>
    output = retrievePattern ("polynomial");

  typedef roboptim::trajectory::Polynomial<3> poly3_t;

  poly3_t::coefs_t params;
  params.setZero ();

  poly3_t null (1., params);
  (*output) << null << std::endl;

  params[0] = -42;
  poly3_t negative_constant (2., params);
  (*output) << negative_constant << std::endl;

  params[0] = 42;
  poly3_t constant (2., params);
  (*output) << constant << std::endl;

  params[1] = 1;
  poly3_t linear1 (0., params);
  (*output) << linear1 << std::endl;

  poly3_t linear2 (3., params);
  (*output) << linear2 << std::endl;

  params[0] = -42;
  poly3_t linear3 (0., params);
  (*output) << linear3 << std::endl;

  params[1] = -1;
  poly3_t linear4 (3., params);
  (*output) << linear4 << std::endl;

  params[0] = 42;
  poly3_t linear5 (3., params);
  (*output) << linear5 << std::endl;

  std::cout << output->str () << std::endl;

  BOOST_CHECK (output->match_pattern ());
}

BOOST_AUTO_TEST_CASE (trajectory_polynomial)
{
  srand (static_cast<unsigned int> (time (NULL)));

  test_copy<3> ();
  test_copy<5> ();

  test_evaluate<3> ();
  test_evaluate<5> ();

  test_derivative<3> ();
  test_derivative<5> ();

  test_translate<3> ();
  test_translate<5> ();

  test_multiply<3> ();
  test_multiply<5> ();

  test_roots<3> ();
  test_roots<5> ();

  test_min<3> ();
  test_min<5> ();

  test_max<3> ();
  test_max<5> ();

  test_misc<3> ();
  test_misc<5> ();

  test_print ();
}

BOOST_AUTO_TEST_SUITE_END ()
