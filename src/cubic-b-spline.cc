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

#include <exception>

#include <boost/numeric/conversion/converter.hpp>
#include <boost/format.hpp>

#include <roboptim/core/indent.hh>
#include <roboptim/core/alloc.hh>

#include <roboptim/trajectory/sys.hh>
#include <roboptim/trajectory/cubic-b-spline.hh>

namespace roboptim
{
namespace trajectory
{

  //FIXME: defined_lc_in has to be true (false untested).
  CubicBSpline::CubicBSpline (interval_t tr, size_type outputSize,
			      const vector_t& p,
			      std::string name, bool clamped)
    : Trajectory<3> (tr, outputSize, p, name),
      nbp_ (p.size () / outputSize), uniform_ (true)
  {
    // Parameter size should be a multiple of spline dimension.
    assert (parameters_.size () % outputSize == 0);

    // Number of control points should be at least 4.
    assert (nbp_ >= 4);

    // Fill vector of regularly spaced knots.
    size_type m = nbp_ + 4;

    value_type delta_t = (tr.second - tr.first) / (static_cast<value_type> (m) - 7.);

    knots_.resize (m);
    size_type kv_iter = 0;

    // Clamped B-spline.
    if (clamped)
      {
	// The first 4 knots should be equal to tr.first.
	// The 4th one will be added in the main loop.
	for (size_type i = 0; i < 3; i++) {
	  knots_[kv_iter++] = tr.first;
	}

	// Note: we do not use an accumulator to get improved numerical precision
	for (size_type i = 0; i < nbp_ - 2; i++) {
	  knots_[kv_iter++] = tr.first + static_cast<value_type> (i) * delta_t;
	}

	// The last 4 knots should be equal to tr.second.
	// The 1st one was added in the main loop.
	for (size_type i = 0; i < 3; i++) {
	  knots_[kv_iter++] = tr.second;
	}
      }
    else // Default case.
      {
	// Note: we do not use an accumulator to get improved numerical precision
	for (size_type i = 0; i < m; i++) {
	  knots_[kv_iter++] = tr.first + static_cast<value_type> (i-3) * delta_t;
	}
      }

    assert (knots_.size () == m);

    // interval lower bound should be rigorously equal to knot 3.
    assert (knots_ [3] == tr.first);

    setParameters (p);
    computeBasisPolynomials ();
  }

  CubicBSpline::CubicBSpline (size_type outputSize, const knots_t& knots,
			      const vector_t& p, std::string name)
    // Note: cast to unsigned to avoid strict-overflow warning
    : Trajectory<3> (std::make_pair (knots[3],
                                     knots[static_cast<unsigned> (knots.size ()-4)]),
		     outputSize, p, name), nbp_ (p.size () / outputSize),
      knots_ (knots), uniform_ (false)
  {
    // Parameter size should be a multiple of spline dimension.
    assert (parameters_.size () % outputSize == 0);

    // Number of control points should be at least 4.
    assert (nbp_ >= 4);

    // Fill vector of regularly spaced knots.
    assert (knots_.size () == nbp_ + 4);
    setParameters (p);
    computeBasisPolynomials ();
  }

  CubicBSpline::CubicBSpline (const CubicBSpline& spline)
    : Trajectory<3> (spline.timeRange (), spline.outputSize (),
                     spline.parameters (), spline.getName ()),
      nbp_ (spline.parameters ().size () / spline.outputSize ()),
      knots_ (spline.knots_),
      basisPolynomials_ (spline.basisPolynomials_),
      uniform_ (spline.uniform_)
  {
    // Parameter size should be a multiple of spline dimension.
    assert (parameters_.size () % outputSize () == 0);

    // Number of control points should be at least 4.
    assert (nbp_ >= 5);

    setParameters (spline.parameters ());
  }


  void CubicBSpline::computeBasisPolynomials ()
  {
    basisPolynomials_.clear();
    for (size_type j = 0; j < nbp_; j++)
      {
	basisPolynomials_.push_back
	  (basisPolynomials_t ());
	// t_j
	value_type t0 = knots_[j];
	// t_{j+1}
	value_type t1 = knots_[j + 1];
	// t_{j+2}
	value_type t2 = knots_[j + 2];
	// t_{j+3}
	value_type t3 = knots_[j + 3];
	// t_{j+4}
	value_type t4 = knots_[j + 4];

	Polynomial3 B0 =
	  1./((t3-t0)*(t2-t0)*(t1-t0))
	  * Monomial3 (t0) * Monomial3 (t0) * Monomial3 (t0);

	Polynomial3 B1 =
	  -1./((t3-t0)*(t2-t1)*(t2-t0))
	  * Monomial3 (t0) * Monomial3 (t0) * Monomial3 (t2)
	  -1./((t3-t0)*(t3-t1)*(t2-t1))
	  * Monomial3 (t0) * Monomial3 (t3) * Monomial3 (t1)
	  -1./((t4-t1)*(t3-t1)*(t2-t1))
	  * Monomial3 (t4) * Monomial3 (t1) * Monomial3 (t1);

	Polynomial3 B2 =
	  1./((t3-t0)*(t3-t1)*(t3-t2))
	  * Monomial3 (t0) * Monomial3 (t3) * Monomial3 (t3)
	  + 1./((t4-t1)*(t3-t1)*(t3-t2))
	  * Monomial3 (t4) * Monomial3 (t1) * Monomial3 (t3)
	  + 1./((t4-t1)*(t4-t2)*(t3-t2))
	  * Monomial3 (t4) * Monomial3 (t4) * Monomial3 (t2);

	Polynomial3 B3 =
	  -1./((t4-t1)*(t4-t2)*(t4-t3))
	  * Monomial3 (t4) * Monomial3 (t4) * Monomial3 (t4);

	basisPolynomials_.back ().push_back (B0);
	basisPolynomials_.back ().push_back (B1);
	basisPolynomials_.back ().push_back (B2);
	basisPolynomials_.back ().push_back (B3);
      }
  }

  CubicBSpline::~CubicBSpline ()
  {
  }

  void
  CubicBSpline::setParameters (const vector_t& p)
  {
    assert (p.size () == parameters_.size ());
    parameters_ = p;
  }

  void
  CubicBSpline::impl_compute (result_ref derivative, double t) const
  {
    t = detail::fixTime (t, *this);
    assert (timeRange ().first <= t && t <= timeRange ().second);
    this->derivative (derivative, t, 0);
  }

  CubicBSpline::value_type
  CubicBSpline::Dt () const
  {
    assert (uniform_);
    return length () / ((value_type)nbp_ - 3.);
  }

  CubicBSpline::size_type
  CubicBSpline::interval (value_type t) const
  {
    t = detail::fixTime (t, *this);

    // Use unsigned to avoid strict-overflow warning
    unsigned i = 0;
    unsigned imin = 3;
    unsigned imax = static_cast<unsigned> (knots_.size () - 5);

    typedef boost::numeric::converter<unsigned, double> Double2Unsigned;

    // In the uniform case, we can access the interval directly
    if (uniform_)
      {
	size_type m = nbp_ + 4;
	double delta_t = (timeRange ().second - timeRange ().first) / (static_cast<double> (m) - 7.);

	i = imin + Double2Unsigned::convert (std::floor ((t - timeRange ().first)/delta_t));
      }
    else // Default case: use dichotomy
      {
	bool found = false;

	if (t == knots_[imax+1]) //We are exactly at the end of the spline
	{
		return imax;
	}
	while (!found)
	  {
	    i = Double2Unsigned::convert
	      (std::floor (.5 * static_cast<double> (imin + imax) + .5));
	    if (t < knots_ [i])
	      {
		imax = i - 1;
	      }
	    else if (t >= knots_ [i + 1])
	      {
		imin = i + 1;
	      }
	    else
	      {
		found = true;
	      }
	    assert (imin <= imax);
	  }
      }

    if (i > nbp_ - 1)
      i = static_cast<unsigned> (nbp_) - 1;
    if (i < 3)
      i = 3;

    assert (knots_ [i] <= t);
    assert (t <= knots_ [i+1]);

    return i;
  }

  CubicBSpline::vector_t
  CubicBSpline::basisFunctions (value_type t, size_type order) const
  {
    assert (uniform_);
    t = detail::fixTime (t, *this);

    const double Dt = length () / ((value_type)nbp_ - 3.);
    const double t3 = getLowerBound (timeRange ());

    const size_type i = interval (t);

    // Nonzero basis functions are b_{i-3,3}(t), b_{i-2,3}(t), b_{i-1,3}(t),
    // b_{i,3}(t).
    const value_type t_i = t3 + ((value_type)i - 3) * Dt;

    const value_type tau_i = t - t_i;
    const value_type tau_i_2 = tau_i * tau_i;
    const value_type tau_i_3 = tau_i * tau_i_2;
    assert (tau_i <= Dt + tolerance ());

    const value_type Dt_2 = Dt * Dt;
    const value_type Dt_3 = Dt * Dt_2;

    // Evaluate basis functions, up to division by 6 Dt_3 that will
    // be performed at the end.
    vector_t b_i (4);
    switch (order)
      {
      case 0:
	{
	  b_i[3] = (  -tau_i_3 + 3*Dt*tau_i_2 - 3*Dt_2*tau_i +   Dt_3);
	  b_i[2] = ( 3*tau_i_3 - 6*Dt*tau_i_2                + 4*Dt_3);
	  b_i[1] = (-3*tau_i_3 + 3*Dt*tau_i_2 + 3*Dt_2*tau_i +   Dt_3);
	  b_i[0] = (   tau_i_3                                       );
	  b_i /= 6.;
	}
	break;
      case 1:
	{
	  b_i[3] = (-3*tau_i_2 +  6*Dt*tau_i - 3*Dt_2);
	  b_i[2] = ( 9*tau_i_2 - 12*Dt*tau_i         );
	  b_i[1] = (-9*tau_i_2 +  6*Dt*tau_i + 3*Dt_2);
	  b_i[0] = ( 3*tau_i_2                       );
	  b_i /= 6.;
	}
	break;
      case 2:
	{
	  b_i[3] = (  -tau_i +  Dt);
	  b_i[2] = ( 3*tau_i - 2*Dt);
	  b_i[1] = (-3*tau_i +   Dt);
	  b_i[0] = (   tau_i       );
	  // division by six not needed here.
	}
	break;
      default:
	assert (0);
      }
    return b_i;
  }

  void
  CubicBSpline::impl_derivative (derivative_ref derivative, value_type t,
				 size_type order)
    const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    bool cur_malloc_allowed = is_malloc_allowed ();
    set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    t = detail::fixTime (t, *this);
    const size_type k = interval (t);
    const std::size_t k_ = static_cast<std::size_t> (k);
    const size_type n = outputSize ();

    const_vector_ref P_k_3 = parameters().segment((k - 3) * n,n);
    const_vector_ref P_k_2 = parameters().segment((k - 2) * n,n);
    const_vector_ref P_k_1 = parameters().segment((k - 1) * n,n);
    const_vector_ref P_k   = parameters().segment((k - 0) * n,n);

    const Polynomial3& B_k_3_k = basisPolynomials_[k_ - 3][3];
    const Polynomial3& B_k_2_k = basisPolynomials_[k_ - 2][2];
    const Polynomial3& B_k_1_k = basisPolynomials_[k_ - 1][1];
    const Polynomial3& B_k_k = basisPolynomials_[k_][0];

    derivative = B_k_3_k.derivative(t, order) * P_k_3 +
      B_k_2_k.derivative(t, order) * P_k_2 +
      B_k_1_k.derivative(t, order) * P_k_1 +
      B_k_k.derivative(t, order) * P_k;

#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    set_is_malloc_allowed (cur_malloc_allowed);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION
  }


  void
  CubicBSpline::translateBasisPolynomials (double t1)
  {
    for (basisPolynomialsVector_t::iterator
           it = basisPolynomials_.begin ();
	 it != basisPolynomials_.end ();
	 ++it)
      for (basisPolynomials_t::iterator
	     p = it->begin ();
	   p != it->end ();
	   ++p)
	p->translateInPlace (t1);
  }

  void
  CubicBSpline::toPolynomials (basisPolynomials_t& res)
    const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    bool cur_malloc_allowed = is_malloc_allowed ();
    set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    const std::size_t nbp = static_cast<std::size_t> (nbp_);
    const size_type n = outputSize ();
    const size_type offset = n - 1;

    // Cubic B-spline ---> nbp_-3 intervals (aka segments or bays)
    if (res.size () != nbp - 3)
      res.resize (nbp - 3);

    if (n >= 3)
    {
      throw std::runtime_error
        ("Invalid use of toPolynomials: dimension must be 1 or 2");
    }

    for (size_type k = 3; k < nbp_; ++k)
      {
	const std::size_t k_ = static_cast<std::size_t> (k);

	const value_type& P_k_3 = parameters()(n * (k - 3) + offset);
	const value_type& P_k_2 = parameters()(n * (k - 2) + offset);
	const value_type& P_k_1 = parameters()(n * (k - 1) + offset);
	const value_type& P_k   = parameters()(n * (k - 0) + offset);

	const Polynomial3& B_k_3_k = basisPolynomials_[k_ - 3][3];
	const Polynomial3& B_k_2_k = basisPolynomials_[k_ - 2][2];
	const Polynomial3& B_k_1_k = basisPolynomials_[k_ - 1][1];
	const Polynomial3& B_k_k   = basisPolynomials_[k_][0];

	res[k_- 3] =
	  P_k_3 * B_k_3_k
	  + P_k_2 * B_k_2_k
	  + P_k_1 * B_k_1_k
	  + P_k   * B_k_k;
      }

#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    set_is_malloc_allowed (cur_malloc_allowed);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION
  }


  CubicBSpline CubicBSpline::operator+ (const CubicBSpline& s) const
  {
    if (parameters ().size () != s.parameters ().size ()
        || timeRange () != s.timeRange ()
        || outputSize () != s.outputSize ())
      throw std::runtime_error ("mismatch in cubic B-spline dimensions.");

    // Since the splines have the same dimension, the sum of the splines
    // is simply the spline with the sum of their control points.
    CubicBSpline sum (timeRange (), outputSize (),
                      parameters () + s.parameters (),
                      this->getName () + " + " + s.getName ());

    return sum;
  }


  void CubicBSpline::operator+= (const CubicBSpline& s)
  {
    if (parameters ().size () != s.parameters ().size ()
        || timeRange () != s.timeRange ()
        || outputSize () != s.outputSize ())
      throw std::runtime_error ("mismatch in cubic B-spline dimensions.");

    // Since the splines have the same dimension, the sum of the splines
    // is simply the spline with the sum of their control points.
    setParameters (parameters () + s.parameters ());
  }


  void
  CubicBSpline::impl_derivative (derivative_ref derivative,
				 StableTimePoint stp,
				 size_type order) const
  {
    this->impl_derivative (derivative,
			   stp.getTime (this->timeRange ()),
			   order);
  }

  CubicBSpline::jacobian_t
  CubicBSpline::variationConfigWrtParam (double t) const
  {
    return variationDerivWrtParam (t, 0);
  }

  CubicBSpline::jacobian_t
  CubicBSpline::variationDerivWrtParam (double t, size_type order)
    const
  {
    t = detail::fixTime (t, *this);
    const size_type k = interval (t);
    const std::size_t k_ = static_cast<std::size_t> (k);
    const size_type n = outputSize ();

    const Polynomial3& B_k_3_k = basisPolynomials_[k_ - 3][3];
    const Polynomial3& B_k_2_k = basisPolynomials_[k_ - 2][2];
    const Polynomial3& B_k_1_k = basisPolynomials_[k_ - 1][1];
    const Polynomial3& B_k_k = basisPolynomials_[k_][0];

    jacobian_t jac (n, nbp_ * n);
    jac.setZero ();
    matrix_t In (n,n);
    In.setIdentity ();

    jac.middleCols((k - 3) * n, n) =
      B_k_3_k.derivative(t, order) * In;
    jac.middleCols((k - 2) * n, n) =
      B_k_2_k.derivative(t, order) * In;
    jac.middleCols((k - 1) * n, n) =
      B_k_1_k.derivative(t, order) * In;
    jac.middleCols((k + 0) * n, n) =
      B_k_k.derivative(t, order) * In;

    return jac;
  }

  CubicBSpline::jacobian_t
  CubicBSpline::variationConfigWrtParam (StableTimePoint stp)
    const
  {
    return this->variationConfigWrtParam (stp.getTime (this->timeRange ()));
  }


  CubicBSpline::jacobian_t
  CubicBSpline::variationDerivWrtParam (StableTimePoint stp, size_type order)
    const
  {
    return this->variationDerivWrtParam
      (stp.getTime (this->timeRange ()), order);
  }


  CubicBSpline::value_type
  CubicBSpline::singularPointAtRank (size_type rank) const
  {
    return (value_type)rank * length () / ((value_type)nbp_- 3);
  }

  CubicBSpline::vector_t
  CubicBSpline::derivBeforeSingularPoint (size_type rank, size_type order) const
  {
    return derivative (singularPointAtRank (rank), order);
  }

  CubicBSpline::vector_t
  CubicBSpline::derivAfterSingularPoint (size_type rank, size_type order) const
  {
    return derivative (singularPointAtRank (rank), order);
  }

  std::ostream&
  CubicBSpline::print (std::ostream& o) const
  {
    using roboptim::operator <<;

    o << "Cubic B-spline:" << incindent
      << iendl << "Name: " << getName ()
      << iendl << "Number of parameters per spline function: " << nbp_
      << iendl << "Length: " << length ()
      << iendl << "Parameters: " << parameters ()
      << iendl << "Knot vector: " << knots_
      << decindent;
    return o;
  }

  int CubicBSpline::dimension () const
  {
    return 3;
  }
} // end of namespace trajectory.
} // end of namespace roboptim.
