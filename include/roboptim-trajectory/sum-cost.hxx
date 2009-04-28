// Copyright (C) 2009 by Florent Lamiraux, Thomas Moulard, AIST, CNRS, INRIA.
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

#ifndef ROBOPTIM_TRAJECTORY_SUM_COST_HXX
# define ROBOPTIM_TRAJECTORY_SUM_COST_HXX

namespace roboptim
{
  template <typename T>
  SumCost<T>::SumCost (const trajectory_t& traj,
		       stateCost sc,
		       const vector_t& pts) throw ()
    : TrajectoryCost (traj),
      stateCost_ (sc),
      points_ (pts)
  {
  }

  template <typename T>
  SumCost<T>::~SumCost () throw ()
  {
  }

  template <typename T>
  typename SumCost<T>::vector_t
  SumCost<T>::operator () (const vector_t& x) const throw ()
  {
    double result = 0.;

    typedef vector_t::const_iterator citer_t;
    for (citer_t it = points_.begin (); it != points_.end (); ++it)
      result += stateCost_ (trajectory_ (*it))[0];

    vector_t res (m);
    res[0] = result;
    return res;
  }

  template <typename T>
  typename SumCost<T>::gradient_t
  SumCost<T>::gradient (const vector_t& x, int) const throw ()
  {
    gradient_t result (n, m);

    typedef vector_t::const_iterator citer_t;
    for (citer_t it = points_.begin (); it != points_.end (); ++it)
      result += stateCost_.gradient (trajectory_ (*it)) *
	trajectory_.variationDerivWrtParam (*it, 1);
    return result;
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_SUM_COST_HXX
