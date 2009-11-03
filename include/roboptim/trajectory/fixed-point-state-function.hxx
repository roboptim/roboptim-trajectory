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

#ifndef ROBOPTIM_TRAJECTORY_FIXED_POINT_STATE_COST_HXX
# define ROBOPTIM_TRAJECTORY_FIXED_POINT_STATE_COST_HXX
# include <boost/format.hpp>

namespace roboptim
{
  template <typename T>
  FixedPointStateFunction<T>::FixedPointStateFunction
  (const trajectory_t& trajectory,
   boost::shared_ptr<DerivableFunction> function,
   const StableTimePoint tpt,
   size_type order) throw ()
    : DerivableFunction (trajectory.parameters ().size (),
			 function->outputSize (),
			 (boost::format
			  ("fixed point state cost using function ``%1%''")
			  % function->getName ()).str ()),
      trajectory_ (trajectory),
      function_ (function),
      tpt_ (tpt),
      order_ (order)
  {
    assert (function_->inputSize () == trajectory_.outputSize () * (order + 1));
  }

  template <typename T>
  FixedPointStateFunction<T>::~FixedPointStateFunction() throw ()
  {
  }

  template <typename T>
  typename FixedPointStateFunction<T>::size_type
  FixedPointStateFunction<T>::order () const throw ()
  {
    return order_;
  }

  template <typename T>
  void
  FixedPointStateFunction<T>::impl_compute (result_t& res,
				  const argument_t& p) const throw ()
  {
    static boost::shared_ptr<T> updatedTrajectory =
      boost::shared_ptr<trajectory_t>
      (trajectory_->getFixedPointTrajectory ().clone ());
    updatedTrajectory->setParameters (removeScaleFromParameters (p));
    (*function_) (res, updatedTrajectory->state
		  (tpt_.getTime (updatedTrajectory->timeRange ()),
		   this->order_));
  }

  template <typename T>
  void
  FixedPointStateFunction<T>::impl_gradient (gradient_t& grad,
				   const argument_t& p,
				   size_type i) const throw ()
  {
    using namespace boost::numeric::ublas;

    static boost::shared_ptr<T> updatedTrajectory =
      boost::shared_ptr<trajectory_t>
      (trajectory_->getFixedPointTrajectory ().clone ());
    updatedTrajectory->setParameters (removeScaleFromParameters (p));

    const value_type t = tpt_.getTime (updatedTrajectory->timeRange ());

    grad[0] = 0.;
    subrange (grad, 1, grad.size ()) =
      prod (function_->gradient
	    (updatedTrajectory->state (t, this->order_), i),
	    updatedTrajectory->variationStateWrtParam (t, this->order_));
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_FIXED_POINT_STATE_COST_HXX

