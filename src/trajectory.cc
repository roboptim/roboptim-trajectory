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
 * \brief Class Trajectory implementation.
 */

#include "roboptim-trajectory/trajectory.hh"

namespace roboptim
{
  Trajectory::Trajectory (size_type derivo, size_type m) throw ()
    : TwiceDerivableFunction (1, m),
      derivabilityOrder_ (derivo),
      parameters_ (),
      singularPoints_ ()
  {
  }


  Trajectory::~Trajectory () throw ()
  {
  }


  Trajectory::vector_t&
  Trajectory::parameters () throw ()
  {
    return parameters_;
  }


  const Trajectory::vector_t&
  Trajectory::parameters () const throw ()
  {
    return parameters_;
  }


  Trajectory::size_type
  Trajectory::derivabilityOrder () const throw ()
  {
    return derivabilityOrder_;
  }


  Trajectory::vector_t
  Trajectory::state (double t, size_type order) const throw ()
  {
    size_type dimension = m;
    vector_t result ((order + 1) * dimension);

    for (size_type o = 0; o <= order; ++o)
      {
	vector_t df = derivative (t, o);
	for (size_type i = 0; i < dimension; ++i)
	  result (i + dimension * o) = df (i);
    }
    return result;
  }

  Trajectory::jacobian_t
  Trajectory::variationStateWrtParam (double t, size_type order) const throw ()
  {
    size_type dimension = m;
    size_type parameterSize = parameters ().size ();
    jacobian_t result ((dimension + 1) * order, parameterSize);

    for (size_type o = 0; o <= order; ++o)
      {
	jacobian_t jacobian = variationDerivWrtParam (t, order);
	for (size_type i = 0; i < dimension; ++i)
	  for (size_type j = 0; j < parameterSize; ++j)
	    result (i + dimension * order, j) = jacobian (i, j);
      }
    return result;
  }

  Trajectory::size_type
  Trajectory::singularPoints () const throw ()
  {
    return singularPoints_;
  }

} // end of namespace roboptim.
