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

#ifndef ROBOPTIM_TRAJECTORY_TRAJECTORY_COST_HH
# define ROBOPTIM_TRAJECTORY_TRAJECTORY_COST_HH
# include <roboptim/trajectory/sys.hh>

# include <roboptim/trajectory/fwd.hh>
# include <roboptim/core/derivable-function.hh>

namespace roboptim
{
  /// \addtogroup roboptim_meta_function
  /// @{

  /// \brief Meta-function for trajectory cost.
  ///
  /// This is a \f$\mathbb{R}^n \rightarrow \mathbb{R}\f$
  /// meta-function for costs on trajectories.
  ///
  /// It takes a trajectory as its input.
  ///
  /// \tparam T trajectory type
  template <typename T>
  class TrajectoryCost : public DerivableFunction
  {
  public:
    /// \brief Trajectory type.
    typedef T trajectory_t;

    /// \brief Concret class should call this constructor.
    ///
    /// \param traj trajectory on which the cost is computed
    /// \param name function's name
    TrajectoryCost (const trajectory_t& traj,
		    std::string name = std::string ());
    virtual ~TrajectoryCost ();
  protected:
    /// \brief Input trajectory.
    const trajectory_t& trajectory_;
  };

  /// @}

} // end of namespace roboptim.

# include <roboptim/trajectory/trajectory-cost.hxx>
#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_COST_HH
