// Copyright (C) 2012 by Florent Lamiraux, CNRS.
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

# include <roboptim/core/function.hh>

# include <ostream>

namespace roboptim {
  /// \brief Polynomial of degree at most 3.

  /// \f[
  /// P (t) = \sum_{i=0}{3} a_i (t-t_0)^i
  /// \f]
  class Polynomial3
  {
  public:
    typedef Function::size_type size_type;

    /// \brief Default constructor: create a null polynomial.
    Polynomial3 ();
    Polynomial3 (double t0, double a0, double a1, double a2, double a3);
    Polynomial3 operator* (const Polynomial3& poly) const;
    Polynomial3 operator+ (const Polynomial3& poly) const;
    Polynomial3 operator- (const Polynomial3& poly) const;
    Polynomial3 translate (const double &t1) const;

    double operator () (const double& t) const;
    double derivative (const double& t, size_type order = 1) const;

    /// \brief Print Polynomial3.
    std::ostream& print (std::ostream& o) const throw ();

    double coefs_ [4];
    double t0_;
  }; // class Polynomial3

  /// \brief Scalar multiplication of a 3rd degree polynomial.
  Polynomial3 operator* (const double& lambda, const Polynomial3& poly);

  /// \brief Monomial
  /// \f[
  /// M (t) = t-t_0
  /// \f]
  struct Monomial : public Polynomial3
  {
    Monomial (double t0) : Polynomial3 (t0, 0., 1., 0., 0.) {}
  }; // class Monomial

  std::ostream& operator << (std::ostream& o, const Polynomial3& p);

} // namespace roboptim
#endif // ROBOPTIM_TRAJECTORY_POLYNOMIAL_3_HH
