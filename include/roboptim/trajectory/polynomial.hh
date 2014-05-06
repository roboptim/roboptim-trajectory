// Copyright (C) 2012 by Florent Lamiraux, CNRS.
// Copyright (C) 2013 by Alexander Werner, DLR.
// Copyright (C) 2014 by Benjamin Chrétien, LIRMM-CNRS.
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
  /// \brief Polynomial of degree at most N.
  ///
  /// \f[
  /// P (t) = \sum_{i=0}^{N} a_i (t-t_0)^i
  /// \f]
  template <int N>
  class Polynomial
  {
  public:
    typedef Function::size_type size_type;
    typedef Function::value_type value_type;
    typedef Function::vector_t vector_t;
    typedef Eigen::Matrix<value_type, N+1, 1> coefs_t;

    /// \brief Default constructor: return a null polynomial.
    Polynomial ();

    /// \brief Construct of a polynomial from its center and its
    /// coefficients.
    ///
    /// \param t0 polynomial of (t-t₀).
    /// \param coefs polynomial coefficients.
    Polynomial (value_type t0, const vector_t& coefs);

    /// \brief Variadic constructor.
    /// Note: this is a legacy constructor used to keep Polynomial3's API.
    /// \param t0 polynomial of (t-t₀).
    /// \param ... variadic arguments containing [a₀,a₁,...,a_N]
    Polynomial (value_type t0, ...);

    /// \brief Copy constructor of polynomials of different orders.
    ///
    /// \warning If N < M, coefficients of higher degrees will be discarded.
    /// This is only valid if you are considering |t-t₀|«1.
    ///
    /// \tparam M degree of the copied polynomial.
    /// \param p copied polynomial.
    template <int M>
    Polynomial (const Polynomial<M>& p);

    /// \brief Return a new polynomial translated from (t-t₀) to (t-t₁).
    /// \param t1 new center.
    /// \return polynomial translated to (t-t₁)
    Polynomial<N> translate (value_type t1) const;

    /// \brief Translate the polynomial (in place) from (t-t₀) to (t-t₁).
    /// \param t1 new center.
    void translateInPlace (value_type t1);

    /// \brief Compute the derivative polynomial of a given order.
    /// \tparam K derivative order.
    /// \return K-th derivative of the polynomial.
    template <int K>
    Polynomial<N-K> derivative () const;

    /// \brief Evaluate the derivative of a given order.
    ///
    /// \param t time of the evaluation.
    /// \param order order of the derivative.
    ///
    /// \return derivative of a given order evaluated at t.
    value_type derivative (value_type t, size_type order = 1) const;

    /// \brief Multiply polynomials of different orders.
    /// \tparam M order of the polynomial.
    /// \param poly polynomial to multiply.
    /// \return P₀P₁
    template <int M>
    Polynomial<N+M> operator* (const Polynomial<M>& poly) const;

    /// \brief Addition of polynomials
    /// \param poly polynomial to add.
    /// \return P₀ + P₁
    Polynomial<N> operator+ (const Polynomial<N>& poly) const;

    /// \brief Subtraction of polynomials.
    /// \param poly polynomial to subtract.
    /// \return P₀ - P₁
    Polynomial<N> operator- (const Polynomial<N>& poly) const;

    /// \brief Scalar multiplication of a polynomial.
    /// \param lambda scalar.
    /// \return λP
    Polynomial<N> operator* (value_type lambda) const;

    /// \brief Evaluate the polynomial.
    /// \param t point of evaluation.
    /// \return P(t)
    value_type operator () (value_type t) const;

    /// \brief Const getter to coefs.
    ///
    /// \return const reference to polynomial coefficients.
    const coefs_t& coefs () const;

    /// \brief Getter to coefs.
    ///
    /// \return reference to polynomial coefficients.
    coefs_t& coefs ();

    /// \brief Const getter to t0.
    ///
    /// \return t0.
    value_type t0 () const;

    /// \brief Reference to t0.
    ///
    /// \return reference to t0.
    value_type& t0 ();

    /// \brief Get the i-th polynomial coefficient.
    /// \param i number of the coefficient to get.
    /// \return i-th polynomial coefficient αi, with: α_i (t-t₀)^i
    value_type operator [] (int i) const;

    /// \brief Return the real roots of the polynomial.
    /// \warning This function relies on Eigen's experimental polynomial
    /// solver.
    /// \return vector of the real roots of the polynomial.
    /// \throw std::runtime_error invalid polynomial (e.g. null).
    std::vector<value_type> realRoots () const;

    /// \brief Return whether the polynomial is constant.
    /// \return true if the polynomial is constant, false otherwise.
    bool isConstant () const;

    /// \brief Return the order of such a polynomial.
    /// \return order of such a polynomial.
    static int order ()
    {
      return order_;
    }

  private:

    /// \brief vector of polynomial coefficients (ordered from lowest to
    /// highest). The vector is of size N+1.
    coefs_t coefs_;

    /// \brief Point on which the polynomial is centered, i.e. a polynomial
    /// of (t-t₀).
    value_type t0_;

  protected:
    value_type
    impl_derivative
    (value_type t, size_type order, size_type start_coef = 0) const;

    /// \brief Implementation of the compile-time derivative.
    /// \tparam K order of the derivation.
    /// \return K-th derivative.
    template <int K>
    Polynomial<N-K> impl_derivative () const;

    coefs_t
    impl_translate (value_type t1) const;

    /// \brief order of the polynomial.
    static const int order_ = N;

    /// \brief Enum for special polynomials.
    enum special_polynomials
      {
        /// \brief Null polynomial: 0
        all_zero_coefficients = 0,
        /// \brief Monomial polynomial: (t-t₀)
        monomial_coefficients = 1
      };

    /// \brief Special constructor for Monomial<N> and some operators.
    ///
    /// \param t0
    /// \param key one of all_zero_coefficients, monomial_coefficients.
    Polynomial (value_type t0, special_polynomials key);
  }; // class Polynomial

  template <int N>
  Polynomial<N> operator* (typename Polynomial<N>::value_type lambda,
                           const Polynomial<N>& poly);

  /// \brief Monomial
  /// \f[
  /// M (t) = t-t_0
  /// \f]
  template <int N>
  struct Monomial : public Polynomial<N>
  {
    typedef Polynomial<N> parent_t;
    typedef typename parent_t::value_type value_type;

  public:
    /// \brief Constructor of a monomial: (t-t₀)
    ///
    /// \param t0 "center" of the monomial.
    Monomial (value_type t0)
      : parent_t (t0, parent_t::monomial_coefficients)
    {}
  }; // class Monomial
} // end of namespace roboptim.

# include <roboptim/trajectory/polynomial.hxx>
#endif // ROBOPTIM_TRAJECTORY_POLYNOMIAL_HH
