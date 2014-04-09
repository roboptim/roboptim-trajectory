// Copyright (C) 2012 by Florent Lamiraux, CNRS.
// Copyright (C) 2013 by Alexander Werner, DLR.
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

#ifndef ROBOPTIM_TRAJECTORY_POLYNOMIAL_HH
# define ROBOPTIM_TRAJECTORY_POLYNOMIAL_HH
# include <roboptim/core/function.hh>

namespace roboptim
{
  /// Polynomial of degree at most N
  ///
  /// \f[
  /// P (t) = \sum_{i=0}{3} a_i (t-t_0)^i
  /// \f]
  template<int N>
  class Polynomial
  {
  public:
    typedef Function::size_type size_type;
    typedef Function::vector_t vector_t;

    Polynomial (double t0, const  vector_t& coefs);

    /**
     * these operations are degree specific, and are only implemented
     * for N={3,5} FIXME
     */
    Polynomial<N> translate (const double &t1) const;

    double derivative (const double& t, size_type order = 1) const;
    Polynomial<N> operator* (const Polynomial<N>& poly) const;
    Polynomial<N> operator+ (const Polynomial<N>& poly) const;
    Polynomial<N> operator- (const Polynomial<N>& poly) const;
    double operator () (const double& t) const;

    double coefs_[N+1];
    // could not use fixed size Eigen::Matrix here since this would
    // create a direct reference to Eigen
    double t0_;
  protected:
    double
    impl_derivative
    (const double &t, size_type order, size_type start_coef=0) const;

    static const int order_=N;

    enum special_polynomials
      {
	all_zero_coefficients=0,
	monomial_coefficients=1
      };

    /**
     * \brief special constructor for Monomial<N> and some operators
     * \param key one of all_zero_coefficients, monomial_coefficients
     */
    Polynomial (double t0, special_polynomials key);
  }; // class Polynomial


  /// Monomial

  /// \f[
  /// M (t) = t-t_0
  /// \f]
  template<int N>
  struct Monomial : public Polynomial<N>
  {
  public:
    Monomial (double t0)
      : Polynomial<N> (t0,Polynomial<N>::monomial_coefficients)
    {}
  }; // class Monomial
} // end of namespace roboptim.

# include <roboptim/trajectory/polynomial.hxx>
#endif // ROBOPTIM_TRAJECTORY_POLYNOMIAL_HH
