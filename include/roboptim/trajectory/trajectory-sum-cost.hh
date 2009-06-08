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

#ifndef ROBOPTIM_TRAJECTORY_TRAJECTORY_SUM_COST_HH
# define ROBOPTIM_TRAJECTORY_TRAJECTORY_SUM_COST_HH
# include <vector>
# include <roboptim/core/derivable-function.hh>

# include <roboptim/trajectory/fwd.hh>
# include <roboptim/trajectory/state-cost.hh>
# include <roboptim/trajectory/trajectory-cost.hh>

namespace roboptim
{
  /// \addtogroup roboptim_meta_function
  /// @{

  /// \brief Define trajectory cost as the sum of state costs.
  ///
  /// Define a generic cost on a trajectory as the sum of costs
  /// on several states.
  ///
  /// \tparam T trajectory type
  template <typename T>
  class TrajectorySumCost : public TrajectoryCost<T>
  {
  public:
    /// \brief Parent type.
    typedef TrajectoryCost<T> parent_t;
    /// \brief Trajectory type.
    typedef T trajectory_t;
    /// \brief State cost type.
    typedef StateCost<T> stateCost_t;

    /// \brief Import vector type.
    typedef typename parent_t::vector_t vector_t;
    /// \brief Import gradient type.
    typedef typename parent_t::gradient_t gradient_t;
    /// \brief Import discrete interval type.
    typedef typename parent_t::discreteInterval_t discreteInterval_t;

    /// \brief Instantiate from a trajectory, a state cost and a list
    /// of discrete points.
    ///
    /// \param traj trajectory on which the cost is computed
    /// \param statecost state cost object
    /// \param vector point list
    TrajectorySumCost (const trajectory_t& traj,
		       const stateCost_t& statecost,
		       const vector_t& vector) throw ();

    /// \brief Instantiate from a trajectory, a state cost and a list
    /// of discrete points.
    ///
    /// \param traj trajectory on which the cost is computed
    /// \param statecost state cost object
    /// \param interval discrete interval used to generate point list
    TrajectorySumCost (const trajectory_t& traj,
		       const stateCost_t& statecost,
		       const discreteInterval_t& interval) throw ();

    virtual ~TrajectorySumCost () throw ();

    virtual vector_t operator () (const vector_t&) const throw ();
    virtual gradient_t gradient (const vector_t&, int) const throw ();

  protected:
    /// \brief State cost object.
    const stateCost_t& stateCost_;
    /// \brief Point list.
    std::vector<double> points_;
  };

  /// @}

} // end of namespace roboptim.

# include <roboptim/trajectory/trajectory-sum-cost.hxx>
#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_SUM_COST_HH
