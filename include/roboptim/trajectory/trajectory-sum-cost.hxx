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
  TrajectorySumCost<T>::
  TrajectorySumCost (const trajectory_t& gamma,
		     boost::shared_ptr<DerivableFunction> cost,
		     const discreteStableTimePointInterval_t& interval,
		     size_type order) throw ()
    : DerivableFunction (gamma.parameters ().size (),
			 cost->outputSize (),
			 (boost::format ("sum cost using function ``%1%''")
			  % cost->getName ()).str ()),
      trajectory_ (gamma),
      function_ (cost),
      interval_ (interval),
      order_ (order)
  {
    assert (function_->inputSize () == trajectory_.outputSize () * (order + 1));
  }

  template <typename T>
  TrajectorySumCost<T>::~TrajectorySumCost() throw ()
  {
  }

  template <typename T>
  typename TrajectorySumCost<T>::size_type
  TrajectorySumCost<T>::order () const throw ()
  {
    return order_;
  }

  template <typename T>
  void
  TrajectorySumCost<T>::impl_compute (result_t& res,
				      const argument_t& p) const throw ()
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
      Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    static trajectory_t updatedTrajectory = trajectory_;
    updatedTrajectory.setParameters (p);

    // Loop over sample points.
    vector_t cost (1);
    cost.setZero ();

    value_type min =
      this->getLowerBound (interval_).getTime (updatedTrajectory.timeRange ());
    value_type max =
      this->getUpperBound (interval_).getTime (updatedTrajectory.timeRange ());
    value_type step =
      this->getStep (interval_).getTime (updatedTrajectory.timeRange ());

    for (value_type t = min; t < max; t += step)
      {
	(*function_) (cost, updatedTrajectory.state(t, this->order_));
	res += cost;
      }
  }

  template <typename T>
  void
  TrajectorySumCost<T>::impl_gradient (gradient_t& grad, const argument_t& p,
				       size_type i) const throw ()
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
      Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    using namespace boost::numeric::ublas;
    static trajectory_t updatedTrajectory = trajectory_;
    updatedTrajectory.setParameters (p);
    grad.setZero ();

    // Loop over sample points.
    gradient_t gr (grad.size ());

    value_type min =
      this->getLowerBound (interval_).getTime (updatedTrajectory.timeRange ());
    value_type max =
      this->getUpperBound (interval_).getTime (updatedTrajectory.timeRange ());
    value_type step =
      this->getStep (interval_).getTime (updatedTrajectory.timeRange ());

    for (value_type t = min; t < max; t += step)
      {
	gr = function_->gradient
	  (updatedTrajectory.state (t, this->order_), i).adjoint () *
	  updatedTrajectory.variationStateWrtParam (t, this->order_);
	grad += gr;
      }
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_SUM_COST_HXX
