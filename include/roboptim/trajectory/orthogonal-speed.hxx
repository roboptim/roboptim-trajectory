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

#ifndef ROBOPTIM_TRAJECTORY_ORTHOGONAL_SPEED_HXX
# define ROBOPTIM_TRAJECTORY_ORTHOGONAL_SPEED_HXX
# include <cmath>
# include <boost/format.hpp>
# include <boost/numeric/ublas/vector.hpp>
# include <boost/scoped_ptr.hpp>

# include <roboptim/core/finite-difference-gradient.hh>

# include <roboptim/trajectory/trajectory.hh>

namespace roboptim
{
  template <typename T>
  OrthogonalSpeed<T>::OrthogonalSpeed (const T& trajectory) throw ()
    : DerivableFunction (1, 1, "orthogonal speed"),
      trajectory_ (trajectory)
  {}

  template <typename T>
  OrthogonalSpeed<T>::~OrthogonalSpeed () throw ()
  {}


  template <typename T>
  void
  OrthogonalSpeed<T>::impl_compute (result_t& res, const argument_t& t) const throw ()
  {
    using namespace boost::numeric::ublas;
    res.clear ();

    result_t position = trajectory_ (t[0]);
    const value_type& theta = position[2];

    gradient_t speed = trajectory_.derivative (t[0], 1);
    const value_type& xdot = speed[0];
    const value_type& ydot = speed[1];

    res[0] = std::cos (theta) * ydot - std::sin (theta) * xdot;
  }

  template <typename T>
  void
  OrthogonalSpeed<T>::impl_gradient (gradient_t& grad, const argument_t& t, size_type i)
    const throw ()
  {
    grad.clear ();

    result_t position = trajectory_ (t[0]);
    const value_type& theta = position[2];

    gradient_t speed = trajectory_.derivative (t[0], 1);
    const value_type& dx = speed[0];
    const value_type& dy = speed[1];
    const value_type& dtheta = speed[2];

    gradient_t acceleration = trajectory_.derivative (t[0], 2);
    const value_type& ddx = acceleration[0];
    const value_type& ddy = acceleration[1];

    grad[0] =
      std::cos (theta) * (ddy - dtheta * dx)
      - std::sin (theta) * (ddx + dtheta * dy);
  }


  template <typename T>
  LimitOrthogonalSpeed<T>::LimitOrthogonalSpeed (StableTimePoint timePoint,
					   const T& trajectory) throw ()
    : DerivableFunction (trajectory.parameters ().size (), 1,
			 (boost::format ("frontal speed limit (%1%)")
			  % timePoint.getAlpha ()).str ()),
      timePoint_ (timePoint),
      trajectory_ (trajectory)
  {}

  template <typename T>
  LimitOrthogonalSpeed<T>::~LimitOrthogonalSpeed () throw ()
  {}

  template <typename T>
  void
  LimitOrthogonalSpeed<T>::impl_compute (result_t& res, const argument_t& p) const throw ()
  {
    using namespace boost::numeric::ublas;
    res.clear ();

    boost::scoped_ptr<T> updatedTrajectory (trajectory_.clone ());
    updatedTrajectory->setParameters (p);

    OrthogonalSpeed<T> frontalSpeed (*updatedTrajectory);
    vector_t t (1);
    t[0] = this->timePoint_.getTime (updatedTrajectory->timeRange ());
    res = frontalSpeed (t);
  }

  template <typename T>
  void
  LimitOrthogonalSpeed<T>::impl_gradient (gradient_t& grad, const argument_t& p, size_type i)
    const throw ()
  {
    assert (i == 0);
    grad.clear ();

    //FIXME: compute gradient analytically.
    FiniteDifferenceGradient fdfunction (*this);
    fdfunction.gradient (grad, p, 0);
  }

  template <typename T>
  template <typename F, typename CLIST>
  void
  LimitOrthogonalSpeed<T>::addToProblem (const T& trajectory,
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
	shared_ptr<LimitOrthogonalSpeed> speed
	  (new LimitOrthogonalSpeed (i * tMax, trajectory));
	problem.addConstraint
	  (static_pointer_cast<DerivableFunction> (speed),
	   vRange);
      }
  }

} // end of namespace roboptim.


#endif //! ROBOPTIM_TRAJECTORY_ORTHOGONAL_SPEED_HXX
