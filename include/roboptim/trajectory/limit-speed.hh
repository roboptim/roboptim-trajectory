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

#ifndef ROBOPTIM_TRAJECTORY_LIMIT_SPEED_HH
# define ROBOPTIM_TRAJECTORY_LIMIT_SPEED_HH
# include <roboptim/core/derivable-function.hh>
# include <roboptim/trajectory/fwd.hh>

namespace roboptim
{
  struct LimitSpeed : public DerivableFunction
  {
    LimitSpeed (double t, const GenericTrajectory& spline) throw ();
    ~LimitSpeed () throw ();
  protected:
    void impl_compute (result_t& res, const argument_t& p) const throw ();
    void impl_gradient (gradient_t& grad, const argument_t& p, size_type i)
      const throw ();
  private:
    double t_;
    const GenericTrajectory& spline_;
  };
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_LIMIT_SPEED_HH
