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

#ifndef ROBOPTIM_TRAJECTORY_LIMIT_OMEGA_HXX
# define ROBOPTIM_TRAJECTORY_LIMIT_OMEGA_HXX
# include <roboptim/trajectory/sys.hh>

# include <boost/format.hpp>
# include <boost/scoped_ptr.hpp>

# include <roboptim/core/finite-difference-gradient.hh>

# include <roboptim/trajectory/trajectory.hh>

namespace roboptim
{
  template <typename T>
  LimitOmega<T>::LimitOmega (StableTimePoint timePoint,
			     const T& trajectory)
    : DerivableFunction (trajectory.parameters ().size (), 1,
			 (boost::format ("omega limit (%1%)")
			  % timePoint.getAlpha ()).str ()),
      timePoint_ (timePoint),
      trajectory_ (trajectory)
  {}

  template <typename T>
  LimitOmega<T>::~LimitOmega ()
  {}

  template <typename T>
  template <typename F, typename CLIST>
  void
  LimitOmega<T>::addToProblem (const T& trajectory,
			       Problem<F, CLIST>& problem,
			       typename Function::interval_t vRange,
			       unsigned nConstraints)
  {
    using namespace boost;
    if (nConstraints == 0)
      return;

    const value_type delta = 1. / nConstraints;

    for (double i = delta; i < 1. - delta; i += delta)
      {
	shared_ptr<LimitOmega> speed
	  (new LimitOmega (i * tMax, trajectory));
	problem.addConstraint
	  (static_pointer_cast<DerivableFunction> (speed),
	   vRange);
      }
  }

  template <typename T>
  void
  LimitOmega<T>::impl_compute (result_t& res, const argument_t& p) const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
      Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    static T updatedTrajectory  = trajectory_;
    updatedTrajectory.setParameters (p);
    res[0] = updatedTrajectory.derivative (timePoint_, 1)[2];
  }

  template <typename T>
  void
  LimitOmega<T>::impl_gradient
  (gradient_t& grad, const argument_t& p, size_type i)
    const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
      Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    GenericFiniteDifferenceGradient<EigenMatrixDense> fd (*this);
    fd.gradient (grad, p, i);
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_LIMIT_OMEGA_HXX
