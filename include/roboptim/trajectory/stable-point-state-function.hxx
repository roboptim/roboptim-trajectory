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

#ifndef ROBOPTIM_TRAJECTORY_STABLE_POINT_STATE_COST_HXX
# define ROBOPTIM_TRAJECTORY_STABLE_POINT_STATE_COST_HXX
# include <boost/format.hpp>

namespace roboptim
{
namespace trajectory
{
  template <typename T>
  StablePointStateFunction<T>::StablePointStateFunction
  (const trajectory_t& trajectory,
   boost::shared_ptr<DerivableFunction> function,
   const StableTimePoint tpt,
   size_type order)
    : DerivableFunction (trajectory.parameters ().size (),
			 function->outputSize (),
			 (boost::format
			  ("stable point state cost using function ``%1%''")
			  % function->getName ()).str ()),
      trajectory_ (trajectory),
      function_ (function),
      tpt_ (tpt),
      order_ (order)
  {
    assert (function_->inputSize () == trajectory_.outputSize () * (order + 1));
  }

  template <typename T>
  StablePointStateFunction<T>::~StablePointStateFunction()
  {
  }

  template <typename T>
  typename StablePointStateFunction<T>::size_type
  StablePointStateFunction<T>::order () const
  {
    return order_;
  }

  template <typename T>
  void
  StablePointStateFunction<T>::impl_compute (result_ref res,
				  const_argument_ref p) const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
      Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    static boost::shared_ptr<trajectory_t> updatedTrajectory =
      boost::shared_ptr<trajectory_t> (trajectory_.clone ());
    updatedTrajectory->setParameters (p);
    (*function_) (res, updatedTrajectory->state (tpt_, this->order_));
  }

  template <typename T>
  void
  StablePointStateFunction<T>::impl_gradient (gradient_ref grad,
				   const_argument_ref p,
				   size_type i) const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
      Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    assert (i == 0);

    static boost::shared_ptr<trajectory_t> updatedTrajectory =
      boost::shared_ptr<trajectory_t> (trajectory_.clone ());
    updatedTrajectory->setParameters (p);

    // Compute derivatives w.r.t parameters.
    // Derivative w.r.t p_0 is wrong here.
    const value_type t = tpt_.getTime (updatedTrajectory->timeRange ());
    grad =
      function_->gradient
      (updatedTrajectory->state (t, this->order_), i)
      * updatedTrajectory->variationStateWrtParam (t, this->order_);

    // Compute derivatives w.r.t p_0.
    const gradient_t df_dstate =
      function_->gradient (updatedTrajectory->state (tpt_, this->order_), i);
    const vector_t dgamma_dt =
      updatedTrajectory->getFixedTimeTrajectory ().state (tpt_, this->order_);
    grad[0] = df_dstate.dot (dgamma_dt);
  }

} // end of namespace trajectory.
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_STABLE_POINT_STATE_COST_HXX
