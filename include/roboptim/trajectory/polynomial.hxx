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

#ifndef ROBOPTIM_TRAJECTORY_POLYNOMIAL_HXX
# define ROBOPTIM_TRAJECTORY_POLYNOMIAL_HXX

namespace roboptim
{
  template<int N>
  Polynomial<N>
  operator* (const double& lambda, const Polynomial<N>& poly)
  {
    Eigen::Matrix<double,N+1,1> temp;
    for(int idx=0;idx<N+1;idx++)
      temp[idx] = lambda * poly.coefs_ [idx];
    return Polynomial<N> (poly.t0_,temp);
  }

  template<int N>
  std::ostream& operator<< (std::ostream& stream, const Polynomial<N>& p)
  {
    stream << "f = ";
    bool printed_before=false;
    for (int idx = N; idx > 0; idx--)
      {
	if(std::abs (p.coefs_[idx]) > std::numeric_limits<double>::epsilon ())
	  {
	    if (printed_before)
	      stream << " + ";
	    stream << p.coefs_[idx] << " * x**" << idx;
	    printed_before = true;
	  }
      }
    if (std::abs (p.coefs_[0]) > std::numeric_limits<double>::epsilon ())
      {
	if (printed_before)
	  stream << " + ";
	stream << p.coefs_[0];
      }
    else if (!printed_before)
      {
	stream << p.coefs_[0];
      }
    stream << "\t" << "( x = t - " << p.t0_ << " )";
    return stream;
  }


  /**
   * start_coef is only used in translate
   */
  template <int N>
  double
  Polynomial<N>::impl_derivative
  (const double &t, size_type order, size_type start_coef) const
  {
    if(order == 0 && start_coef == 0)
      return (*this) (t);

    double dt = 1.;
    double result = 0.;

    size_type idx = order;
    if (start_coef > idx)
      {
	for (int diff = 0; diff< start_coef - idx ; diff++)
	  dt = dt * (t - t0_);
	idx = start_coef;
      }

    // loop over all polynomial coefficients which are in the
    // derivative
    for (; idx <= order_; idx++)
      {
	// determine factor of this coefficient for the requested
	// derivative
	double deriv_coef = 1.;
	for (size_type foo = 0, exponent = idx;
	     foo < order; foo++, exponent--)
	  deriv_coef *= static_cast<double> (exponent);

	result += deriv_coef * coefs_[idx] * dt;

	dt = dt * ( t - t0_ );
      }
    return result;
  }

  template <int N>
  Polynomial<N> Polynomial<N>::translate (const double &t1) const
  {
    Polynomial<N> result (t1, all_zero_coefficients);
    for(int order = 0; order <= order_; order++)
      {
	result.coefs_[order] = coefs_[order];
	int factorial = 1;
	for(int idx = order; idx > 0; idx--)
	  factorial *= idx; // calculate order!
	const double remaining_derivs = impl_derivative (t1, order, order + 1);
	result.coefs_[order] += remaining_derivs / factorial;
      }
    return result;
  }

  template <int N>
  double Polynomial<N>::derivative(const double &t, size_type order) const
  {
    return impl_derivative(t, order);
  }

  template <int N>
  Polynomial<N> Polynomial<N>::operator* (const Polynomial<N>& poly) const
  {
    Polynomial<order_> other = poly.translate (t0_);
    Eigen::Matrix<double, order_+1, 1> temp;
    temp.setZero ();

    /* unrolling should be possible here since the loops are only
     * dependant on constants and the if/else goes away in release mode
     */
    for (size_type idx_a = 0; idx_a < order_ + 1; idx_a++)
      {
	for (size_type idx_b = 0; idx_b < order_ + 1; idx_b++)
	  {
	    const size_type combined_idx = idx_a + idx_b;
	    if (combined_idx < order_ + 1)
	      temp[combined_idx] += coefs_[idx_a] * other.coefs_[idx_b];
	    else
	      // Those combinations whould yield a Polynomial of higher order
	      assert (coefs_[idx_a] * other.coefs_[idx_b] == 0);
	  }
      }
    return Polynomial<order_> (t0_, temp);
  }

  /**
   * implementations of general Polynomial<N> operators follow
   */
  template <int N>
  Polynomial<N> Polynomial<N>::operator+ (const Polynomial<N>& poly) const
  {
    Polynomial<N> other = poly.translate (t0_);
    Eigen::Matrix<double, order_ + 1, 1> temp;
    for (int idx = 0; idx  <N + 1; idx++)
      temp[idx] = coefs_ [idx] + other.coefs_ [idx];
    return Polynomial<N> (t0_,temp);
  }

  template <int N>
  Polynomial<N> Polynomial<N>::operator- (const Polynomial<N>& poly) const
  {
    Polynomial<N> other = poly.translate (t0_);
    Eigen::Matrix<double,order_+1,1> temp;
    for(int idx=0;idx<N+1;idx++)
      temp[idx] = coefs_ [idx] + other.coefs_ [idx];
    return Polynomial<N> (t0_,temp);
  }

  template <int N>
  double Polynomial<N>::operator () (const double& t) const
  {
    double dt = 1.;
    double result = 0.;
    for(int idx=0;idx<N+1;idx++)
      {
	result += coefs_[idx] * dt;
	dt = dt * ( t - t0_ );
      }
    return result;
  }

  template <int N>
  Polynomial<N>::Polynomial (double t0, const  vector_t& coefs)
    : t0_(t0)
  {
    assert (coefs.rows () == N + 1);
    std::copy (coefs.data (), coefs.data () + coefs.rows (), coefs_);
  }

  template <int N>
  Polynomial<N>::Polynomial (double t0, special_polynomials key)
    : t0_ (t0)
  {
    switch(key)
      {
      case all_zero_coefficients:
	std::fill_n(coefs_,order_+1,0);
	break;
      case monomial_coefficients:
	std::fill_n(coefs_,order_+1,0);
	coefs_[1]=1.;
	break;
      default:
	assert(0);
	break;
      }
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_POLYNOMIAL_HXX
