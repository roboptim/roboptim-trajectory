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

#include <roboptim/trajectory/sys.hh>

#include <cmath>
#include <roboptim/trajectory/orthogonal-speed.hh>

namespace roboptim
{
  OrthogonalSpeed::OrthogonalSpeed ()
    : DerivableFunction (2 * 3, 1, "orthogonal speed")
  {}

  OrthogonalSpeed::~OrthogonalSpeed () throw()
  {}

  void
  OrthogonalSpeed::impl_compute (result_t& res, const argument_t& x) 
    const throw()
  {
    const value_type& theta = x[2];

    const value_type& xdot = x[3 + 0];
    const value_type& ydot = x[3 + 1];
    res[0] = std::cos (theta) * ydot - std::sin (theta) * xdot;
  }


  void
  OrthogonalSpeed::impl_gradient (gradient_t& grad, const argument_t& arg, size_type)
    const throw()
  {
    //const value_type& x = arg[0];
    //const value_type& y = arg[1];
    const value_type& theta = arg[2];

    const value_type& dx = arg[3 + 0];
    const value_type& dy = arg[3 + 1];
    //const value_type& dtheta = arg[3 + 2];

    grad.setZero();
    // 0, 1, 5 components are null.
    grad[2] = -1. * (std::cos (theta) * dx + std::sin (theta) * dy);
    grad[3] = -std::sin (theta);
    grad[4] = std::cos (theta);
  }
} // end of namespace roboptim.
