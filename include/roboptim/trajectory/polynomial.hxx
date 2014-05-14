// Copyright (C) 2012 by Florent Lamiraux, CNRS.
// Copyright (C) 2013 by Alexander Werner, DLR.
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

#ifndef ROBOPTIM_TRAJECTORY_POLYNOMIAL_HXX
# define ROBOPTIM_TRAJECTORY_POLYNOMIAL_HXX

# include <algorithm> // transform
# include <limits> // numeric_limits
# include <cstdarg> // va_list, va_start, va_arg, va_end

// Polynomial solver
# include <unsupported/Eigen/Polynomials>

namespace roboptim
{
  template <int N>
  Polynomial<N>::Polynomial ()
    : t0_ (0)
  {
    coefs_.setZero ();
  }

  template <int N>
  Polynomial<N>::Polynomial (value_type t0, const vector_t& coefs)
    : t0_ (t0)
  {
    assert (coefs.rows () == N + 1);
    coefs_ = coefs;
  }

  template <int N>
  Polynomial<N>::Polynomial (value_type t0, ...)
    :  t0_ (t0)
  {
    va_list vl;
    // Note: void va_start (va_list ap, paramN);
    // ---> "paramN: Name of the last named parameter in the function
    // definition. The arguments extracted by subsequent calls to
    // va_arg are those after paramN."
    // See: http://www.cplusplus.com/reference/cstdarg/va_start/
    va_start (vl, t0);

    for (size_type i = 0; i < N + 1; ++i)
      // We expect coefs to be of size N+1.
      coefs_[i] = va_arg (vl, value_type);
    va_end (vl);
  }

  template <int N>
  template <int M>
  Polynomial<N>::Polynomial (const Polynomial<M>& p)
    : t0_ (p.t0 ())
  {
    if (M > N)
      {
	// Warning: this discards coefficients of higher degree.
	coefs_ = p.coefs ().template head<N+1> ();
      }
    else // M <= N
      {
	coefs_.setZero ();
	coefs_.template head<M+1> () = p.coefs ();
      }
  }

  template <int N>
  std::ostream& operator<< (std::ostream& stream, const Polynomial<N>& p)
  {
    typedef typename Polynomial<N>::value_type value_type;
    typedef typename Polynomial<N>::coefs_t coefs_t;
    const coefs_t& coefs = p.coefs ();

    stream << "f = ";
    bool printed_before = false;
    for (int idx = N; idx > 0; idx--)
      {
	if (std::abs (coefs[idx])
	    > std::numeric_limits<value_type>::epsilon ())
	  {
	    if (printed_before)
	      stream << " + ";
	    stream << coefs[idx] << " * x**" << idx;
	    printed_before = true;
	  }
      }
    if (std::abs (coefs[0])
        > std::numeric_limits<value_type>::epsilon ())
      {
	if (printed_before)
	  stream << " + ";
	stream << coefs[0];
      }
    else if (!printed_before)
      {
	stream << coefs[0];
      }
    stream << "\t" << "( x = t - " << p.t0 () << " )";
    return stream;
  }

  template <int N>
  template <int K>
  Polynomial<N-K>
  Polynomial<N>::impl_derivative () const
  {
    assert ((K >= 0 && K <= N) && "Wrong derivation order.");
    if (K == 0) return Polynomial<N> (*this);

    typename Polynomial<N-K>::coefs_t coefs;
    coefs.setZero ();
    value_type alpha = 1.;

    for (size_type i = 0; i < coefs.size (); ++i)
      {
	alpha = 1.;

	for (size_type j = i+K; j > i; --j)
	  alpha *= static_cast<value_type> (j);

	coefs[i] = alpha * coefs_[i+K];
      }

    return Polynomial<N-K> (t0_, coefs);
  }

  // start_coef is only used in translate
  template <int N>
  typename Polynomial<N>::value_type
  Polynomial<N>::impl_derivative
  (value_type t, size_type order, size_type start_coef) const
  {
    if (order > N)
      return 0.;
    else if (order == 0 && start_coef == 0)
      return (*this) (t);

    value_type dt = 1.;
    value_type result = 0.;

    size_type idx = order;
    if (start_coef > idx)
      {
	for (int diff = 0; diff < start_coef - idx ; diff++)
	  dt = dt * (t - t0_);
	idx = start_coef;
      }

    // loop over all polynomial coefficients which are in the
    // derivative
    for (; idx <= order_; idx++)
      {
	// determine factor of this coefficient for the requested
	// derivative
	value_type deriv_coef = 1.;
	for (size_type foo = 0, exponent = idx;
	     foo < order; foo++, exponent--)
	  deriv_coef *= static_cast<value_type> (exponent);

	result += deriv_coef * coefs_[idx] * dt;

	dt = dt * (t - t0_);
      }
    return result;
  }

  template <int N>
  typename Polynomial<N>::coefs_t
  Polynomial<N>::impl_translate (typename Polynomial<N>::value_type t1) const
  {
    coefs_t coefs;
    for (int order = 0; order <= order_; order++)
      {
	coefs[order] = coefs_[order];
	int factorial = 1;
	for (int idx = order; idx > 0; idx--)
	  {
	    // overflow check
	    assert ((factorial <= std::numeric_limits<int>::max ()/idx)
		    && "Overflow detected in factorial computation.");
	    factorial *= idx; // calculate order!
	  }

	const value_type remaining_derivs = impl_derivative (t1,
							     order,
							     order + 1);
	coefs[order] += remaining_derivs / factorial;
      }
    return coefs;
  }

  template <int N>
  Polynomial<N> Polynomial<N>::translate (value_type t1) const
  {
    return Polynomial<N> (t1, impl_translate (t1));
  }

  template <int N>
  void Polynomial<N>::translateInPlace (value_type t1)
  {
    coefs_ = impl_translate (t1);
    t0_ = t1;
  }

  template <int N>
  template <int K>
  Polynomial<N-K> Polynomial<N>::derivative () const
  {
    return impl_derivative<K> ();
  }

  template <int N>
  typename Polynomial<N>::value_type
  Polynomial<N>::derivative (value_type t, size_type order) const
  {
    return impl_derivative (t, order);
  }

  template <int N>
  template <int M>
  Polynomial<N+M> Polynomial<N>::operator* (const Polynomial<M>& poly) const
  {
    Polynomial<M> other = poly.translate (t0_);
    typename Polynomial<N+M>::coefs_t temp;
    temp.setZero ();

    // unrolling should be possible here since the loops are only
    // dependent on constant values
    for (size_type idx_a = 0; idx_a < N + 1; idx_a++)
      {
	for (size_type idx_b = 0; idx_b < M + 1; idx_b++)
	  {
	    temp[idx_a + idx_b] += coefs_[idx_a] * other.coefs ()[idx_b];
	  }
      }
    return Polynomial<N+M> (t0_, temp);
  }

  template <int N>
  Polynomial<N> Polynomial<N>::operator+ (const Polynomial<N>& poly) const
  {
    Polynomial<N> other = poly.translate (t0_);
    return Polynomial<N> (t0_, coefs_ + other.coefs_);
  }

  template <int N>
  Polynomial<N> Polynomial<N>::operator- (const Polynomial<N>& poly) const
  {
    Polynomial<N> other = poly.translate (t0_);
    return Polynomial<N> (t0_, coefs_ - other.coefs_);
  }

  template <int N>
  Polynomial<N> Polynomial<N>::operator* (value_type lambda) const
  {
    return Polynomial<N> (t0_, lambda * coefs_);
  }

  template <int N>
  Polynomial<N> operator* (typename Polynomial<N>::value_type lambda,
                           const Polynomial<N>& poly)
  {
    return poly * lambda;
  }

  template <int N>
  typename Polynomial<N>::value_type
  Polynomial<N>::operator () (value_type t) const
  {
    value_type dt = 1.;
    value_type result = 0.;
    for (int idx = 0; idx < N + 1; idx++)
      {
	result += coefs_[idx] * dt;
	dt = dt * ( t - t0_ );
      }
    return result;
  }

  template <int N>
  typename Polynomial<N>::value_type
  Polynomial<N>::t0 () const
  {
    return t0_;
  }

  template <int N>
  typename Polynomial<N>::value_type&
  Polynomial<N>::t0 ()
  {
    return t0_;
  }

  template <int N>
  const typename Polynomial<N>::coefs_t&
  Polynomial<N>::coefs () const
  {
    return coefs_;
  }

  template <int N>
  typename Polynomial<N>::coefs_t&
  Polynomial<N>::coefs ()
  {
    return coefs_;
  }

  template <int N>
  typename Polynomial<N>::value_type
  Polynomial<N>::operator [] (int i) const
  {
    assert (i >= 0 && i < N+1);
    return coefs_[i];
  }

  template <int N>
  std::vector<typename Polynomial<N>::value_type>
  Polynomial<N>::realRoots () const
  {
    const value_type eps = 1e-6;
    std::vector<value_type> roots;

    // Eigen expects a polynomial in the form: Σ α_i t^i
    // Thus, we compute the roots of the "translated" polynomial:
    // P(u) = Σ α_i u^i
    // Note: the leading coefficient has to be different than 0
    vector_t::Index n = coefs_.size ();

    // Usual case: leading coefficient ≠ 0 (i.e. known polynomial size)
    if (std::abs (coefs_[n-1]) > eps)
      {
	Eigen::PolynomialSolver<value_type, N> solver (coefs_);
	solver.realRoots (roots);
      }
    else // leading coefficient = 0, i.e. dynamic size
      {
	while (n > 0 && std::abs (coefs_[n-1]) < eps)
	  n--;

	if (n == 0)
	  throw std::runtime_error ("solver cannot process null polynomials.");

	Eigen::PolynomialSolver<value_type, Eigen::Dynamic> solver (coefs_.head (n));
	solver.realRoots (roots);
      }

    // Then we shift the roots to get the roots of the actual polynomial:
    // u = t - t₀
    // P(u) = 0 => P(t = u + t₀) = 0
    std::transform (roots.begin (), roots.end (), roots.begin (),
                    std::bind2nd (std::plus<value_type> (), t0_));

    return roots;
  }

  template <int N>
  bool Polynomial<N>::isConstant (value_type epsilon) const
  {
    // Check whether all coefs after the first one are null
    return coefs_.template tail<N> ().template lpNorm<Eigen::Infinity> ()
      < epsilon;
  }

  template <int N>
  Polynomial<N>::Polynomial (value_type t0, special_polynomials key)
    : t0_ (t0)
  {
    switch (key)
      {
      case all_zero_coefficients:
        coefs_.setZero ();
        break;
      case monomial_coefficients:
        coefs_.setZero ();
        coefs_[1] = 1.;
        break;
      default:
        assert (0);
        break;
      }
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_POLYNOMIAL_HXX
