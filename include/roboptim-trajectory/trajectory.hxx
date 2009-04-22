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

#ifndef ROBOPTIM_TRAJECTORY_TRAJECTORY_HXX
# define ROBOPTIM_TRAJECTORY_TRAJECTORY_HXX

namespace roboptim
{
  template <unsigned dorder>
  Trajectory<dorder>::Trajectory (size_type m, const vector_t& p) throw ()
    : parent_t (m),
      parameters_ (p),
      singularPoints_ ()
  {
  }


  template <unsigned dorder>
  Trajectory<dorder>::~Trajectory () throw ()
  {
  }


  template <unsigned dorder>
  typename Trajectory<dorder>::vector_t&
  Trajectory<dorder>::parameters () throw ()
  {
    return parameters_;
  }


  template <unsigned dorder>
  const typename Trajectory<dorder>::vector_t&
  Trajectory<dorder>::parameters () const throw ()
  {
    return parameters_;
  }


  template <unsigned dorder>
  typename Trajectory<dorder>::vector_t
  Trajectory<dorder>::state (double t, size_type order) const throw ()
  {
    size_type dimension = this->m;
    vector_t result ((order + 1) * dimension);

    for (size_type o = 0; o <= order; ++o)
      {
	vector_t df = derivative (t, o);
	for (size_type i = 0; i < dimension; ++i)
	  result (i + dimension * o) = df (i);
    }
    return result;
  }


  template <unsigned dorder>
  typename Trajectory<dorder>::jacobian_t
  Trajectory<dorder>::variationStateWrtParam (double t, size_type order)
    const throw ()
  {
    size_type dimension = this->m;
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


  template <unsigned dorder>
  typename Trajectory<dorder>::size_type
  Trajectory<dorder>::singularPoints () const throw ()
  {
    return singularPoints_;
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HXX
