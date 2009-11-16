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

#ifndef ROBOPTIM_TRAJECTORY_FWD_HH
# define ROBOPTIM_TRAJECTORY_FWD_HH
# include <roboptim/trajectory/sys.hh>

namespace roboptim
{
  template <typename T>
  class FreeTimeTrajectory;

  template <unsigned dorder>
  class Trajectory;

  template <typename T>
  class StateFunction;

  template <typename T>
  class SumCost;

  template <typename T>
  class TrajectoryCost;

  class Spline;
  class SplineLength;
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_FWD_HH
