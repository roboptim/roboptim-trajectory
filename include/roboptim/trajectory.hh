// Copyright (C) 2010 by Thomas Moulard, AIST, CNRS, INRIA.
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

#ifndef ROBOPTIM_TRAJECTORY_HH
# define ROBOPTIM_TRAJECTORY_HH
# include <roboptim/trajectory/sys.hh>

// Generic headers.
# include <roboptim/trajectory/fwd.hh>


// Main headers.
# include <roboptim/trajectory/b-spline.hh>
# include <roboptim/trajectory/constrained-b-spline.hh>
# include <roboptim/trajectory/cubic-b-spline.hh>
# include <roboptim/trajectory/free-time-trajectory.hh>
# include <roboptim/trajectory/freeze.hh>
# include <roboptim/trajectory/frontal-speed.hh>
# include <roboptim/trajectory/limit-omega.hh>
# include <roboptim/trajectory/limit-speed.hh>
# include <roboptim/trajectory/orthogonal-speed.hh>
# include <roboptim/trajectory/polynomial.hh>
# include <roboptim/trajectory/polynomial-3.hh>
# include <roboptim/trajectory/spline-length.hh>
# include <roboptim/trajectory/stable-point-state-function.hh>
# include <roboptim/trajectory/stable-time-point.hh>
# include <roboptim/trajectory/state-function.hh>
# include <roboptim/trajectory/sys.hh>
# include <roboptim/trajectory/trajectory-cost.hh>
# include <roboptim/trajectory/trajectory-sum-cost.hh>
# include <roboptim/trajectory/trajectory.hh>
# include <roboptim/trajectory/vector-interpolation.hh>


// Visualization.
# include <roboptim/trajectory/visualization/limit-speed.hh>
# include <roboptim/trajectory/visualization/speed.hh>
# include <roboptim/trajectory/visualization/trajectory.hh>

#endif //! ROBOPTIM_TRAJECTORY_HH
