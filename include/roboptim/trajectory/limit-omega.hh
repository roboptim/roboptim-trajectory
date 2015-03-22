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

#ifndef ROBOPTIM_TRAJECTORY_LIMIT_OMEGA_HH
# define ROBOPTIM_TRAJECTORY_LIMIT_OMEGA_HH
# include <roboptim/trajectory/sys.hh>

# include <boost/shared_ptr.hpp>

# include <roboptim/core/derivable-function.hh>
# include <roboptim/trajectory/fwd.hh>
# include <roboptim/trajectory/stable-time-point.hh>

namespace roboptim
{
namespace trajectory
{
  template <typename T>
  class LimitOmega : public DerivableFunction
  {
  public:
    LimitOmega (StableTimePoint timePoint,
		const T& trajectory);
    ~LimitOmega ();

    template <typename F, typename CLIST>
    static void addToProblem (const T& trajectory,
			      Problem<F, CLIST>& problem,
			      typename Function::interval_t vRange,
			      unsigned nConstraints);

  protected:
    void impl_compute (result_ref res, const_argument_ref p) const;
    void impl_gradient (gradient_ref grad, const_argument_ref p, size_type i)
      const;

  private:
    StableTimePoint timePoint_;
    const T& trajectory_;
  };

} // end of namespace trajectory.
} // end of namespace roboptim.

# include <roboptim/trajectory/limit-omega.hxx>
#endif //! ROBOPTIM_TRAJECTORY_LIMIT_OMEGA_HH
