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

#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/polynomial-3.hh>

using namespace roboptim;

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (trajectory_polynomial3)
{
  double tol = 1e-6;

  // Polynomial centered on 0
  Polynomial3 p (0., 1., 2., 3., 4.);
  BOOST_CHECK_CLOSE (p (0.), 1., tol);
  BOOST_CHECK_CLOSE (p[0], 1., tol);
  BOOST_CHECK_CLOSE (p[1], 2., tol);
  BOOST_CHECK_CLOSE (p[2], 3., tol);
  BOOST_CHECK_CLOSE (p[3], 4., tol);

  // Polynomial centered on 1
  Polynomial3 q (1., 1., 2., 3., 4.);
  BOOST_CHECK_CLOSE (q (1.), 1., tol);

  // Translated polynomial
  Polynomial3 r = q.translate (0.);
  BOOST_CHECK_CLOSE (q (3.), r (3.), tol);

  // Real roots
  std::vector<double> roots;
  p = Polynomial3 (0., 1., 1., 1., 1.);
  roots = p.realRoots ();
  BOOST_CHECK_EQUAL (roots.size (), 1);
  BOOST_CHECK_CLOSE (roots[0], -1., tol);

  p = Polynomial3 (0., 1., 2., 2., 1.);
  roots = p.realRoots ();
  BOOST_CHECK_EQUAL (roots.size (), 1);
  BOOST_CHECK_CLOSE (roots[0], -1., tol);

  p = Polynomial3 (0., -6., 11., -6., 1.);
  roots = p.realRoots ();
  BOOST_CHECK_EQUAL (roots.size (), 3);
  BOOST_CHECK_CLOSE (roots[0], 1., tol);
  BOOST_CHECK_CLOSE (roots[1], 2., tol);
  BOOST_CHECK_CLOSE (roots[2], 3., tol);

  p = Polynomial3 (1., -6., 11., -6., 1.);
  roots = p.realRoots ();
  BOOST_CHECK_EQUAL (roots.size (), 3);
  BOOST_CHECK_CLOSE (roots[0], 2., tol);
  BOOST_CHECK_CLOSE (roots[1], 3., tol);
  BOOST_CHECK_CLOSE (roots[2], 4., tol);

  p = Polynomial3 (1., -6., 11., -6., 0.);
  roots = p.realRoots ();

  p = Polynomial3 (1., 0., 0., 0., 0.);
  BOOST_CHECK_THROW (roots = p.realRoots (), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END ()
