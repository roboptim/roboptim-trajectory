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

    res[0] = std::cos (theta) * ydot + std::sin (theta) * xdot;
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

    grad[0] = dtheta * (std::cos (theta) * (ddy + dx)
			+ std::sin (theta) * (ddx - dy));
  }
} // end of namespace roboptim.


#endif //! ROBOPTIM_TRAJECTORY_ORTHOGONAL_SPEED_HXX
