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

#include <roboptim/trajectory/polynomial.hh>
#include <roboptim/trajectory/polynomial-3.hh>

#include <limits>
#include <unsupported/Eigen/Polynomials>

using namespace roboptim;
using namespace std;

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

static double tol = 1e-6;

template <int N>
void test_multiply ()
{
  typename Polynomial<N>::coefs_t params;
  params.setRandom ();
  double t0 = (double)rand () / RAND_MAX;
  double t1 = (double)rand () / RAND_MAX;

  Polynomial<N> p_0 (t0, params);
  p_0 (0.1);

  params.setZero ();
  params[0] = 0.3;
  Polynomial<N> p_1 (t1, params);
  p_1 (0.1);

  // Polynomial multiplication
  Polynomial<2*N> p_2 = p_0 * p_1;
  p_2 (0.1);
  // Polynomial multiplication + crop
  Polynomial<N> p_3 = p_0 * p_1;
  // Compare uncropped part
  BOOST_CHECK (allclose (p_2.coefs ().template head<N+1> (),
                         p_3.coefs ()));

  // Scalar multiplication
  Polynomial<N> q_0 = 2.0 * p_0;
  Polynomial<N> q_1 = p_0 * 2.0;

  BOOST_CHECK (allclose (q_0.coefs (), q_1.coefs ()));
}


template <int N>
struct poly_props
{
  static void check_derivative (const Polynomial<N>&, double, int, double)
  {
    // Not implemented for this N
    assert (0);
  }

  static void check_translate (const Polynomial<N>&, double,
			       const Polynomial<N>&)
  {
    // Not implemented for this N
    assert (0);
  }

  static void check_evaluate (double t0,
			      typename Polynomial<N>::coefs_t& coefs,
			      double t)
  {
    Polynomial<N> poly (t0, coefs);

    // Use Eigen's evaluator: translate polynomial from (t-t₀) to t
    Polynomial<N> eigen_poly = poly.translate (0.);
    double res = Eigen::poly_eval (eigen_poly.coefs (), t);

    BOOST_CHECK_CLOSE (poly (t), res, tol);
  }
};


// Code from Polynomial3
template <>
void poly_props<3>::check_derivative (const Polynomial<3>& poly,
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
void poly_props<5>::check_derivative (const Polynomial<5>& poly,
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
void poly_props<3>::check_translate (const Polynomial<3>& poly,
                                     double t1,
                                     const Polynomial<3>& ref_poly)
{
  double dt = t1 - poly.t0 ();
  double dt2 = dt * dt;
  double dt3 = dt * dt2;
  double a0 = poly.coefs () [0];
  double a1 = poly.coefs () [1];
  double a2 = poly.coefs () [2];
  double a3 = poly.coefs () [3];

  Polynomial<3>::coefs_t temp;
  temp [0] = a0 + a1 * dt + a2 * dt2 + a3 * dt3;
  temp [1] = a1 + 2 * dt * a2 + 3 * dt2 * a3;
  temp [2] = a2 + 3 * dt * a3;
  temp [3] = a3;
  Polynomial<3> p_new (t1, temp);

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
void poly_props<5>::check_translate (const Polynomial<5>& poly,
                                     double t1,
                                     const Polynomial<5>& ref_poly)
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


  Polynomial<5>::coefs_t temp;
  temp [0] = a0 +     a1 * dt +      a2 * dt2 +      a3 * dt3 +     a4 * dt4 + a5 * dt5;
  temp [1] = a1 + 2 * a2 * dt +  3 * a3 * dt2 +  4 * a4 * dt3 + 5 * a5 * dt4;
  temp [2] = a2 + 3 * a3 * dt +  6 * a4 * dt2 + 10 * a5 * dt3;
  temp [3] = a3 + 4 * a4 * dt + 10 * a5 * dt2;
  temp [4] = a4 + 5 * a5 * dt;
  temp [5] = a5;
  Polynomial<5> p_new (t1, temp);

  BOOST_CHECK_CLOSE (p_new.t0 (), ref_poly.t0 (), tol);

  for (int idx = 0; idx < 5 + 1; idx++)
    {
      BOOST_CHECK_CLOSE (p_new.coefs ()[idx],
                         ref_poly.coefs ()[idx],
                         tol);
    }
}


template <int N>
void test_derivative ()
{
  typename Polynomial<N>::coefs_t params;
  params.setRandom ();
  double t0 = (double)rand () / RAND_MAX;
  Polynomial<N> p_1 (t0, params);

  double t = (double)rand () / RAND_MAX;

  for (int order = 0; order < N + 2; order++)
    {
      double derivative = p_1.derivative (t, order);

      if (order > N)
        BOOST_CHECK_SMALL (derivative, tol);

      poly_props<N>::check_derivative (p_1, t, order, derivative);
    }
}

template <int N>
void test_translate ()
{
  typename Polynomial<N>::coefs_t params;
  params.setRandom ();
  double t0 = (double)rand () / RAND_MAX;
  double t1 = (double)rand () / RAND_MAX;
  Polynomial<N> p_1 (t0, params);
  Polynomial<N> p_2 = p_1.translate (t1);

  BOOST_CHECK_CLOSE (p_1 (0.), p_2 (0.), tol);

  poly_props<N>::check_translate (p_1, t1, p_2);

  // Test translate in place.
  Polynomial<N> p_3 = p_1;
  p_3.translateInPlace (t1);
  BOOST_CHECK_CLOSE (p_2.t0 (), p_3.t0 (), tol);
  BOOST_CHECK (allclose (p_2.coefs (), p_3.coefs ()));
}

template <int N>
void test_copy ()
{
  typename Polynomial<N>::coefs_t params;
  params.setRandom ();

  // Same degree
  Polynomial<N> p_0 (3., params);
  Polynomial<N> p_1 = p_0;

  BOOST_CHECK_CLOSE (p_0.t0 (), p_1.t0 (), tol);
  BOOST_CHECK (allclose (p_0.coefs (), p_1.coefs ()));

  // Different degrees
  Polynomial<N+2> p_2 = p_0;
  BOOST_CHECK_CLOSE (p_0.t0 (), p_2.t0 (), tol);
  BOOST_CHECK (allclose (p_0.coefs (),
                         p_2.coefs ().template head<N+1> ()));

  Polynomial<N-1> p_3 = p_0;
  BOOST_CHECK_CLOSE (p_0.t0 (), p_3.t0 (), tol);
  BOOST_CHECK (allclose (p_3.coefs (),
                         p_0.coefs ().template head<N> ()));
}

template <int N>
void test_evaluate ()
{
  typename Polynomial<N>::coefs_t params;
  params.setRandom ();
  double t0 = (double)rand () / RAND_MAX;
  double t = (double)rand () / RAND_MAX;
  poly_props<N>::check_evaluate (t0, params, t);
}

template <int N>
void test_roots ()
{
  typedef typename Polynomial<N>::value_type value_type;

  typename Polynomial<N>::coefs_t params;
  params.setRandom ();
  double t0 = (double)rand () / RAND_MAX;
  Polynomial<N> p (t0, params);

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
}

template <int N>
void test_misc ()
{
  typename Polynomial<N>::coefs_t params;
  params.setZero ();
  Polynomial<N> p (1., params);

  BOOST_CHECK (p.isConstant ());
  p.coefs ()[1] = 2.;
  BOOST_CHECK (!p.isConstant ());
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

  test_misc<3> ();
  test_misc<5> ();
}

BOOST_AUTO_TEST_SUITE_END ()
