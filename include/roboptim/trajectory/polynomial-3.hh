// Copyright (C) 2012 by Florent Lamiraux, CNRS.
// Copyright (C) 2014 by Benjamin Chrétien, CNRS-LIRMM.
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

#ifndef ROBOPTIM_TRAJECTORY_POLYNOMIAL_3_HH
# define ROBOPTIM_TRAJECTORY_POLYNOMIAL_3_HH

# include <roboptim/trajectory/polynomial.hh>

namespace roboptim
{
  namespace trajectory
  {

    /// \brief Polynomial of degree at most 3.
    /// \f[
    /// P (t) = \sum_{i=0}^{3} a_i (t-t_0)^i
    /// \f]
    typedef Polynomial<3> Polynomial3;

    /// \brief Monomial3
    /// \f[
    /// M (t) = t-t_0
    /// \f]
    typedef Monomial<3> Monomial3;

  } // end of namespace trajectory.
} // end of namespace roboptim.

#endif // ROBOPTIM_TRAJECTORY_POLYNOMIAL_3_HH
