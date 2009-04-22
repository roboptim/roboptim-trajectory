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

/**
 * \brief Class Spline implentation.
 */

#include <roboptim-trajectory/spline.hh>

#include <spline/bspline.h>

namespace roboptim
{
  Spline::Spline (size_type m, const vector_t& p) throw ()
    : Trajectory<4> (m, p)
  {
  }

  Spline::~Spline () throw ()
  {
  }

  Spline::vector_t
  Spline::derivative (double x, size_type order) const throw ()
  {
    vector_t res (42);
    return res;
  }

  Spline::jacobian_t
  Spline::variationConfigWrtParam (double t) const throw ()
  {
    jacobian_t jac (42, 42);
    return jac;
  }

  Spline::jacobian_t
  Spline::variationDerivWrtParam (double t, size_type order)
    const throw ()
  {
    jacobian_t jac (42, 42);
    return jac;
  }

  Spline::value_type
  Spline::singularPointAtRank (size_type rank) const
  {
    return 42.;
  }

  Spline::vector_t
  Spline::derivBeforeSingularPoint (size_type rank, size_type order) const
  {
    vector_t res (42);
    return res;
  }

  Spline::vector_t
  Spline::derivAfterSingularPoint (size_type rank, size_type order) const
  {
    vector_t res (42);
    return res;
  }
} // end of namespace roboptim.
