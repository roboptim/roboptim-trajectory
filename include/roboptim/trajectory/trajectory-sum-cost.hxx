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
		     const discreteInterval_t& interval,
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
    static trajectory_t updatedTrajectory = trajectory_;
    updatedTrajectory.setParameters (p);

    // Loop over sample points.
    // TODO replace by stable time points
    vector_t cost(1);
    cost.clear();
    for (value_type t=interval_.get<0>(); t < interval_.get<1>(); 
	 t += interval_.get<2>()) {
      (*function_) (cost, updatedTrajectory.state(t, this->order_));
      res += cost;
    }
  }

  template <typename T>
  void
  TrajectorySumCost<T>::impl_gradient (gradient_t& grad, const argument_t& p, 
				       size_type i) const throw ()
  {
    using namespace boost::numeric::ublas;
    static trajectory_t updatedTrajectory = trajectory_;
    updatedTrajectory.setParameters (p);
    grad.clear();
    // Loop over sample points.
    // TODO replace by stable time points
    gradient_t gr(grad.size());
    for (value_type t=interval_.get<0>(); t < interval_.get<1>(); 
	 t += interval_.get<2>()) {
      gr = prod (function_->gradient (updatedTrajectory.state (t, this->order_),
				      i),
		 updatedTrajectory.variationStateWrtParam (t, this->order_));
      grad += gr;
    }
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_SUM_COST_HXX
