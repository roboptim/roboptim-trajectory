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
# include <boost/shared_ptr.hpp>

# include <roboptim/core/derivable-function.hh>
# include <roboptim/trajectory/fwd.hh>
# include <roboptim/trajectory/stable-time-point.hh>

namespace roboptim
{
  struct LimitSpeed : public DerivableFunction
  {
    LimitSpeed (StableTimePoint timePoint, const GenericTrajectory& spline) throw ();
    ~LimitSpeed () throw ();

    template <typename F, typename CLIST>
    static void addToProblem (const GenericTrajectory&,
			      Problem<F, CLIST>&,
			      typename Function::interval_t);

  protected:
    void impl_compute (result_t& res, const argument_t& p) const throw ();
    void impl_gradient (gradient_t& grad, const argument_t& p, size_type i)
      const throw ();
  private:
    StableTimePoint timePoint_;
    const GenericTrajectory& trajectory_;
  };

  template <typename F, typename CLIST>
  void
  LimitSpeed::addToProblem (const GenericTrajectory& trajectory,
			    Problem<F, CLIST>& problem,
			    typename Function::interval_t vRange)
  {
    using namespace boost;
    for (double i = 0; i < 1.; i += 0.1)
      {
	shared_ptr<LimitSpeed> speed (new LimitSpeed (i * tMax, trajectory));
	problem.addConstraint
	  (static_pointer_cast<DerivableFunction> (speed),
	   vRange);
      }
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_LIMIT_SPEED_HH
