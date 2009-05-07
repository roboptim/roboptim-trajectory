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

#ifndef ROBOPTIM_TRAJECTORY_TRAJECTORY_SUM_COST_HXX
# define ROBOPTIM_TRAJECTORY_TRAJECTORY_SUM_COST_HXX
# include <boost/tuple/tuple.hpp>

namespace roboptim
{
  template <typename T>
  TrajectorySumCost<T>::TrajectorySumCost (const trajectory_t& traj,
					   const stateCost_t& sc,
					   const vector_t& pts) throw ()
    : parent_t (traj),
      stateCost_ (sc),
      points_ (pts)
  {
  }

  template <typename T>
  TrajectorySumCost<T>::TrajectorySumCost
  (const trajectory_t& traj,
   const stateCost_t& sc,
   const discreteInterval_t& interval) throw ()
    : parent_t (traj),
      stateCost_ (sc),
      points_ ()
  {
    using namespace boost;
    for (double i = get<0> (interval); i <= get<1> (interval);
	 i += get<2> (interval))
      points_.push_back (i);
  }


  template <typename T>
  TrajectorySumCost<T>::~TrajectorySumCost () throw ()
  {
  }

  template <typename T>
  typename TrajectorySumCost<T>::vector_t
  TrajectorySumCost<T>::operator () (const vector_t& x) const throw ()
  {
    double result = 0.;

    typedef typename std::vector<double>::const_iterator citer_t;
    for (citer_t it = points_.begin (); it != points_.end (); ++it)
      result += stateCost_ (this->trajectory_ (*it))[0];

    vector_t res (this->m);
    res[0] = result;
    return res;
  }

  template <typename T>
  typename TrajectorySumCost<T>::gradient_t
  TrajectorySumCost<T>::gradient (const vector_t& x, int) const throw ()
  {
    gradient_t result (this->n);
    result.clear ();

    typedef typename std::vector<double>::const_iterator citer_t;
    for (citer_t it = points_.begin (); it != points_.end (); ++it)
      {
	assert (this->trajectory_.timeRange ().first <= *it
		&& *it <= this->trajectory_.timeRange ().second);

	result += prod (stateCost_.gradient (this->trajectory_ (*it), 0),
			this->trajectory_.variationConfigWrtParam (*it));
      }
    return result;
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_SUM_COST_HXX
