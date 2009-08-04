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

#ifndef ROBOPTIM_TRAJECTORY_FRONTAL_SPEED_HXX
# define ROBOPTIM_TRAJECTORY_FRONTAL_SPEED_HXX
# include <cmath>
# include <boost/format.hpp>
# include <boost/numeric/ublas/vector.hpp>
# include <boost/scoped_ptr.hpp>

# include <roboptim/core/finite-difference-gradient.hh>

# include <roboptim/trajectory/trajectory.hh>

namespace roboptim
{
  namespace detail
  {
    std::string getFrontalSpeedName (const StableTimePoint& timePoint);

    std::string getFrontalSpeedName (const StableTimePoint& timePoint)
    {
      using boost::format;
      return (format ("frontal speed (%1%)") % timePoint.getAlpha ()).str ();
    }
  }


  template <typename T>
  FrontalSpeed<T>::FrontalSpeed (StableTimePoint timePoint,
			     const T& trajectory) throw ()
    : DerivableFunction (trajectory.parameters ().size (), 1,
			 detail::getFrontalSpeedName (timePoint)),
      timePoint_ (timePoint),
      trajectory_ (trajectory)
  {}

  template <typename T>
  FrontalSpeed<T>::~FrontalSpeed () throw ()
  {}

  template <typename T>
  const T&
  FrontalSpeed<T>::trajectory () const throw ()
  {
    return trajectory_;
  }

  template <typename T>
  void
  FrontalSpeed<T>::impl_compute (result_t& res, const argument_t& p) const throw ()
  {
    using namespace boost::numeric::ublas;
    res.clear ();

    boost::scoped_ptr<T> updatedTrajectory (trajectory_.clone ());
    updatedTrajectory->setParameters (p);

    result_t position = (*updatedTrajectory) (timePoint_);
    // const value_type& x = position[0];
    // const value_type& y = position[1];
    const value_type& theta = position[2];

    gradient_t speed = updatedTrajectory->derivative (timePoint_, 1);
    const value_type& xdot = speed[0];
    const value_type& ydot = speed[1];
    // const value_type& thetadot = speed[2];

    res[0] = std::sin (theta) * ydot - std::cos (theta) * xdot;
  }

  template <typename T>
  void
  FrontalSpeed<T>::impl_gradient (gradient_t& grad, const argument_t& p, size_type i)
    const throw ()
  {
    FiniteDifferenceGradient fdfunction (*this);
    fdfunction.gradient (grad, p, 0);
  }

  template <typename T>
  template <typename F, typename CLIST>
  void
  FrontalSpeed<T>::addToProblem (const T& trajectory,
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
	shared_ptr<FrontalSpeed> speed (new FrontalSpeed (i * tMax, trajectory));
	problem.addConstraint
	  (static_pointer_cast<DerivableFunction> (speed),
	   vRange);
      }
  }
} // end of namespace roboptim.


#endif //! ROBOPTIM_TRAJECTORY_FRONTAL_SPEED_HXX
