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

#ifndef ROBOPTIM_TRAJECTORY_TRAJECTORY_COST_HXX
# define ROBOPTIM_TRAJECTORY_TRAJECTORY_COST_HXX

namespace roboptim
{
  template <typename T>
  TrajectoryCost<T>::TrajectoryCost (const trajectory_t& traj,
				     std::string name)
    : DerivableFunction (traj.parameters ().size (), 1, name),
      trajectory_ (traj)
  {
  }

  template <typename T>
  TrajectoryCost<T>::~TrajectoryCost () throw()
  {
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_COST_HXX
