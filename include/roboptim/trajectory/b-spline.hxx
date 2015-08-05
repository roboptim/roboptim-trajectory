// Copyright (C) 2010 by Thomas Moulard, AIST, CNRS, INRIA.
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

#ifndef ROBOPTIM_TRAJECTORY_B_SPLINE_HXX
# define ROBOPTIM_TRAJECTORY_B_SPLINE_HXX
# include <boost/shared_ptr.hpp>
# include <boost/numeric/conversion/converter.hpp>
# include <roboptim/core/numeric-linear-function.hh>

namespace roboptim
{
namespace trajectory
{
  //FIXME: defined_lc_in has to be true (false untested).
  template <int N>
  BSpline<N>::BSpline (const interval_t& tr, size_type outputSize,
                       const vector_t& p,
                       std::string name, bool clamped)
    : Trajectory<N> (tr, outputSize, p, name),
    nbp_ (p.size () / outputSize),
    knots_ (),
    basisPolynomials_ (),
    uniform_ (true)
  {
    //Parameter size should be a multiple of spline dimension
    assert (this->parameters_.size () % outputSize == 0);
    // number of control points should be at least order + 1 so that
    // there is at least one segment in spline
    assert (nbp_ >= order_ + 1);

    // Initialize the knot vector
    initializeKnots (clamped);

    setParameters (p);
    computeBasisPolynomials ();
  }

  template <int N>
  BSpline<N>::BSpline (const interval_t& tr, size_type outputSize,
                       const vector_t& p, const_vector_ref knots,
                       std::string name)
    : Trajectory<N> (tr, outputSize, p, name),
    nbp_ (p.size () / outputSize), knots_ (knots), uniform_ (false)
  {
    //Parameter size should be a multiple of spline dimension
    assert (this->parameters_.size () % outputSize == 0);
    // number of control points should be at least 6.
    assert (nbp_ >= order_ + 1);
    // calculated number of control points must match the recieved
    // control point set
    assert (nbp_ + order_  + 1 == knots.size());
    // control points must be a monotonically increasing series
    for (size_type idx = 0; idx < knots.size () - 1; idx++)
      assert (knots[idx] <= knots[idx + 1]);

    // not more then the order_ of knot points with the same value are allowed

    setParameters (p);
    computeBasisPolynomials ();
  }

  template <int N>
  BSpline<N>::BSpline (const BSpline<N>& spline)
    : Trajectory<N> (spline.timeRange (), spline.outputSize (),
		     spline.parameters_),
    nbp_ (spline.parameters_.size () / spline.outputSize ()),
    knots_ (spline.knots_),
    basisPolynomials_ (spline.basisPolynomials_),
    uniform_ (spline.uniform_)
  {
    // Parameter size should be a multiple of spline dimension
    assert (this->parameters_.size () % this->outputSize () == 0);
    // number of control points should be at least 6.
    assert (nbp_ >= order_ + 1);

    setParameters (spline.parameters ());
  }

  template <int N>
  BSpline<N>::BSpline (ConstructionMode mode,
                       const interval_t& tr, size_type outputSize,
                       const vector_t& p,
                       std::string name, bool clamped)
    : Trajectory<N> (tr, outputSize, p, name),
    nbp_ (p.size () / outputSize),
    knots_ (),
    basisPolynomials_ (),
    uniform_ (true)
  {
    // Parameter size should be a multiple of spline dimension
    assert (this->parameters_.size () % this->outputSize () == 0);
    // number of control points should be at least 6.
    assert (nbp_ >= order_ + 1);

    setParameters (p);

    if (mode == NORMAL)
      {
	initializeKnots (clamped);
	computeBasisPolynomials ();
      }
  }

  template <int N>
  template <int M>
  BSpline<N> BSpline<N>::derivative () const
  {
    assert (M > 0);

    // Copy interval, dimension, name, number of control points, knots
    // and compute the basis polynomials by simple derivation
    std::string name = (boost::format ("%s (order %i derivative)")
			% this->getName () % M).str ();

    // TODO: deduce basis polynomials instead
    BSpline<N> ds (UNINITIALIZED,
		   this->timeRange (), this->outputSize (),
		   this->parameters_, name);
    ds.knots_ = this->knots_;
    ds.basisPolynomials_ = this->deriveBasisPolynomials<M> ();

    return ds;
  }

  template <int N>
  Trajectory<N>* BSpline<N>::resize (interval_t timeRange) const
  {
    return new BSpline<N> (timeRange,
                           this->outputSize (),
                           this->parameters (),
                           this->knots_,
                           this->getName());
  }


  template <int N>
  void BSpline<N>::initializeKnots (bool clamped)
  {
    // Fill vector of regularly spaced knots
    size_type m = nbp_ + order_ + 1;
    knots_.resize (m);

    const interval_t& tr = this->timeRange ();

    value_type delta_t =
      (tr.second - tr.first)
      / static_cast<value_type> (nbp_ - order_);

    if (clamped)
      {
	// The first order_+1 knots should be equal to tr.first.
	// The last one will be added in the main loop.
	for (size_type i = 0; i < order_; i++) {
	  knots_ (i) = tr.first;
	}

	// Note: we do not use an accumulator to get improved numerical precision
	for (size_type i = 0; i < nbp_ - order_ + 1; i++) {
	  knots_ (order_ + i) = tr.first + static_cast<double> (i) * delta_t;
	}

	// The last order_+1 knots should be equal to tr.second.
	// The 1st one was added in the main loop.
	for (size_type i = 0; i < order_; i++) {
	  knots_ (m-1-i) = tr.second;
	}
      }
    else
      {
	for (size_type i = 0; i < m; i++)
	  {
	    knots_ (i) = tr.first + static_cast<value_type> (i-order_) * delta_t;
	  }
      }
  }


  template <int N>
  template <int M>
  typename BSpline<N>::basisPolynomialsVector_t
  BSpline<N>::deriveBasisPolynomials () const
  {
    basisPolynomialsVector_t result (static_cast<size_t> (nbp_));
    assert (result.size () == basisPolynomials_.size ());

    for (size_t i = 0; i < basisPolynomials_.size (); ++i)
      {
	const basisPolynomials_t& polynomials = basisPolynomials_[i];
	for (typename basisPolynomials_t::const_iterator
	       p  = polynomials.begin ();
	     p != polynomials.end (); ++p)
	  {
	    // Conversion to proper polynomial degree done automatically
	    result[i].push_back (p->template derivative<M> ());
	  }
      }

    return result;
  }


  template <int N>
  typename BSpline<N>::cox_map
  BSpline<N>::cox_de_boor (size_type j, size_type n) const
  {
    std::stringstream label;
    label << "[j=" << j << "/n=" << n << "]";

    if (n == 0) //end of recursion
      {
	cox_map map;
	Eigen::Matrix < value_type, N + 1, 1 > temp_params;
	temp_params.setZero();
	temp_params[0] = 1.;
	std::pair<cox_map_itr_t, bool> ptr =
	  map.insert (std::make_pair (j, polynomial_t (0., temp_params)));

	LOG4CXX_DEBUG (this->logger, label.str() << "return(end)    : "
		       << ptr.first->second);
	return map;
      }
    else
      {
	// t_{j}
	const value_type t0 = knots_[j + 0];
	// t_{j+1}
	const value_type t1 = knots_[j + 1];
	// t_{j+n}
	const value_type tn = knots_[j + n];
	// t_{j+n+1}
	const value_type tn1 = knots_[j + n + 1];

	polynomial_t p_1_rat = 1. / (tn - t0) * monomial_t (t0);
	// http://wolftype.com/ucsb/spatial/bspline.html, uniform b-splines
	if (std::isinf (p_1_rat.coefs ()[1]))
	  {
	    //FIXME: this is probably not how it should work.
	    p_1_rat = monomial_t (t0);
	  }

	LOG4CXX_DEBUG (this->logger,
		       label.str() << "p_1_rat        : "
		       << p_1_rat);

	cox_map p_1_cox = cox_de_boor (j, n - 1);
	cox_map p_1;

	for (cox_map_itr_t itr = p_1_cox.begin (); itr != p_1_cox.end(); itr++)
	  {
	    LOG4CXX_DEBUG (this->logger,
			   label.str()
			   << "p_1_cox        : " << itr->first
			   << " : " << itr->second);

	    polynomial_t p_prod = itr->second * p_1_rat;

	    LOG4CXX_DEBUG (this->logger,
			   label.str() << "p_prod (1)     : "
			   << p_prod);
	    p_1.insert (std::make_pair (itr->first, p_prod));
	  }

	polynomial_t p_2_rat = 1. / (t1 - tn1) * monomial_t (tn1);

	// http://wolftype.com/ucsb/spatial/bspline.html, uniform b-splines
	if (std::isinf (p_2_rat.coefs ()[1]))
	  {
	    //FIXME: this is probably not how it should work.
	    p_2_rat = Monomial<N> (tn1);
	  }

	LOG4CXX_DEBUG (this->logger,
		       label.str() << "p_2_rat        : "
		       << p_2_rat);

	cox_map p_2_cox = cox_de_boor (j + 1, n - 1);
	cox_map p_2;
	for (cox_map_itr_t itr = p_2_cox.begin(); itr != p_2_cox.end(); itr++)
	  {
	    LOG4CXX_DEBUG (this->logger,
			   label.str()
			   << "p_2_cox        : " << itr->first
			   << " : " << itr->second);

	    polynomial_t p_prod = itr->second * p_2_rat;

	    LOG4CXX_DEBUG (this->logger,
			   label.str()
			   << "p_prod (2)     : " << p_prod);

	    p_2.insert (std::make_pair (itr->first, p_prod));
	  }

	cox_map p;
	p.insert (p_1.begin(), p_1.end());
	for (cox_map_itr_t itr = p_2.begin(); itr != p_2.end(); itr++)
	  {
	    cox_map_itr_t existing_element = p.find (itr->first);
	    if (existing_element != p.end())
	      existing_element->second = existing_element->second + itr->second;
	    else
	      p.insert (std::make_pair (itr->first, itr->second));
	  }

        LOG4CXX_DEBUG (this->logger,
                       label.str() << "result of recursion branch");

	for (cox_map_itr_t itr = p.begin(); itr != p.end(); itr++)
	  {
	    LOG4CXX_DEBUG (this->logger,
			   label.str() << "B_" << j
			   << "_" << itr->first << "  : "
			   << itr->second);
	  }
	return p;
      }
  }

  template <int N>
  void BSpline<N>::toPolynomials (basisPolynomials_t& res) const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    bool cur_malloc_allowed = is_malloc_allowed ();
    set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    const std::size_t nbp = static_cast<std::size_t> (nbp_);
    const size_type n = this->outputSize ();
    const size_type offset = n - 1;

    // B-spline ---> nbp_-order_ intervals (aka segments or bays)
    if (res.size () != nbp - order_)
      res.resize (nbp - order_);

    if (n >= order_)
    {
      throw std::runtime_error
        ("Invalid use of toPolynomials: dimension must be 1 or 2");
    }

    for (size_type k = order_; k < nbp_; ++k)
      {
	const std::size_t k_ = static_cast<std::size_t> (k);

	res[k_ - order_] = polynomial_t ();
  for (unsigned long i = 0; i <= static_cast<unsigned long> (order_); ++i)
  {
    res[k_ - order_] +=
      basisPolynomials_[k_-i][i]*this->parameters()(n*(k-static_cast<long>(i))+offset);
  }
      }

#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    set_is_malloc_allowed (cur_malloc_allowed);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION
  }


  template <int N>
  void BSpline<N>::computeBasisPolynomials ()
  {
    basisPolynomials_.clear();
    for (size_type j = 0; j < nbp_; j++)
      {
	//calculate basis polynomials for each interval of the knot vector
	basisPolynomials_.push_back (basisPolynomials_t ());
	cox_map map = cox_de_boor (j, order_);

	for (cox_map_itr_t itr = map.begin(); itr != map.end(); itr++)
	  {
	    LOG4CXX_DEBUG (this->logger,
			   "B_" << j << "_"
			   << itr->first << "  : " << itr->second);

	    basisPolynomials_.back ().push_back (itr->second);
	  }
      }
  }

  template <int N>
  void BSpline<N>::setParameters (const vector_t& p)
  {
    assert (p.size () == this->parameters_.size ());
    this->parameters_ = p;
  }

  template <int N>
  void BSpline<N>::impl_compute (result_ref derivative, value_type t) const
  {
    t = detail::fixTime (t, *this);
    assert (this->timeRange ().first <= t && t <= this->timeRange ().second);
    this->derivative (derivative, t, 0);
  }

  template <int N>
  typename BSpline<N>::value_type
  BSpline<N>::Dt () const
  {
    assert (uniform_);
    return this->length () / ((value_type)nbp_ - static_cast<value_type> (order_));
  }

  template <int N>
  typename BSpline<N>::size_type
  BSpline<N>::interval (value_type t) const
  {
    t = detail::fixTime (t, *this);
    typedef boost::numeric::converter<size_type, value_type> Double2SizeType;

    size_type i = 0;
    size_type imin = order_;
    size_type imax = nbp_;

    // In the uniform case, we can access the interval directly
    if (uniform_)
      {
        double delta_t = (this->timeRange ().second
                          - this->timeRange ().first)
	  / (static_cast<double> (nbp_ - order_));

        i = imin + Double2SizeType::convert
          (std::floor ((t - this->timeRange ().first)/delta_t));
      }
    else
      {
	bool found = false;

	if (t == knots_[imax]) //We are exactly at the end of the spline
	{
		return imax;
	}
	while (!found)
	  {
	    i = Double2SizeType::convert
	      (std::floor (.5 * static_cast<double> (imin + imax) + .5));

	    if (t < knots_ [i])
	      {
		imax = i - 1;
	      }
	    else if (t >= knots_ [i + 1])
	      {
		if (t < knots_ [i + 2])
		  {
		    i = i + 1;
		    found = true;
		  }
		else
		  {
		    imin = i + 1;
		  }
	      }
	    else
	      {
		found = true;
	      }
	    assert (imin <= imax);
	  }
      }

    if (i > nbp_ - 1)
      i = nbp_ - 1;
    if (i < order_)
      i = order_;

    assert (knots_ [i] <= t);
    assert (t <= knots_ [i+1]);

    return i;
  }

  template <int N>
  typename BSpline<N>::vector_t
  BSpline<N>::basisFunctions (value_type t, size_type order) const
  {
    vector_t result (order_ + 1);

    t = detail::fixTime (t, *this);
    const size_type k = interval (t);
    const size_type n = this->outputSize ();

    for (size_type idx = 0; idx < order_ + 1; idx++)
      {
	const Polynomial<N>& B = basisPolynomials_[k - idx][idx];
	result (idx) =  B.derivative (t, order);
      }
  }

  template <int N> void
  BSpline<N>::impl_derivative (derivative_ref derivative, value_type t,
			       size_type order)
    const
  {

    t = detail::fixTime (t, *this);
    const size_type k = interval (t);
    const size_type n = this->outputSize ();

    derivative.setZero();

    ROBOPTIM_DEBUG_ONLY(value_type polynomial_sum = 0.);
    for (size_type idx = 0; idx < order_ + 1; ++idx)
      {
	const std::size_t k_ = static_cast<std::size_t> (k);
	const std::size_t idx_ = static_cast<std::size_t> (idx);

	const_vector_ref P_seg = this->parameters_.segment ((k - idx) * n, n);
	const Polynomial<N>& B = basisPolynomials_[k_ - idx_][idx_];

	ROBOPTIM_DEBUG_ONLY(polynomial_sum += B.derivative (t, order));

	derivative +=  B.derivative (t, order) * P_seg;
      }

    // this is true for any knot vector (FIXME: reference)
    if (order == 0)
      {
        // Note: if we deal with derived splines, then the
        // sum is 0.
        assert
          (std::abs (polynomial_sum - 1.) < 1e-8
           || std::abs (polynomial_sum) < 1e-8);
      }
  }

  template <int N> void
  BSpline<N>::impl_derivative (derivative_ref derivative,
			       StableTimePoint stp,
			       size_type order) const
  {
    this->impl_derivative (derivative,
			   stp.getTime (this->timeRange ()),
			   order);
  }

  template <int N>
  typename BSpline<N>::jacobian_t
  BSpline<N>::variationConfigWrtParam (value_type t) const
  {
    return variationDerivWrtParam (t, 0);
  }


  template <int N>
  typename BSpline<N>::jacobian_t
  BSpline<N>::variationDerivWrtParam (value_type t, size_type order)
    const
  {

    t = detail::fixTime (t, *this);
    const size_type k = interval (t);
    const size_type n = this->outputSize ();

    jacobian_t jac (n, nbp_ * n);
    jac.setZero ();
    matrix_t In (n, n);
    In.setIdentity ();

    for (std::size_t idx = 0; idx < static_cast<std::size_t> (order_ + 1); ++idx)
      {
	const std::size_t k_ = static_cast<std::size_t> (k);
	const Polynomial<N>& B = basisPolynomials_[k_ - idx][idx];

	jac.middleCols ((k - static_cast<size_type> (idx)) * n, n) =
	  B.derivative (t, order) * In;
      }

    return jac;
  }

  template <int N>
  typename BSpline<N>::jacobian_t
  BSpline<N>::variationConfigWrtParam (StableTimePoint stp)
    const
  {
    return this->variationConfigWrtParam (stp.getTime (this->timeRange ()));
  }


  template <int N>
  typename BSpline<N>::jacobian_t
  BSpline<N>::variationDerivWrtParam (StableTimePoint stp, size_type order)
    const
  {
    return this->variationDerivWrtParam
      (stp.getTime (this->timeRange ()), order);
  }

  template <int N>
  typename BSpline<N>::value_type
  BSpline<N>::singularPointAtRank (size_type rank) const
  {
    return (value_type)rank * this->length () / ((value_type)nbp_ - order_);
  }

  template <int N>
  typename BSpline<N>::vector_t
  BSpline<N>::derivBeforeSingularPoint (size_type rank, size_type order) const
  {
    return this->derivative (singularPointAtRank (rank), order);
  }

  template <int N>
  typename BSpline<N>::vector_t
  BSpline<N>::derivAfterSingularPoint (size_type rank, size_type order) const
  {
    return this->derivative (singularPointAtRank (rank), order);
  }

  template <int N> std::ostream&
  BSpline<N>::print (std::ostream& o) const
  {
    using roboptim::operator <<;

    o << "Order " << order_ << " B-spline:" << incindent
      << iendl << "Name: " << this->getName ()
      << iendl << "Number of parameters per spline function: " << nbp_
      << iendl << "Length: " << this->length ()
      << iendl << "Knot vector: " << knots_
      << iendl << "Parameters: " << this->parameters ()
      << decindent;
    return o;
  }

  template <int N>
  const typename BSpline<N>::vector_t& BSpline<N>::knotVector () const
  {
    return this->knots_;
  }

  template <int N>
  int BSpline<N>::dimension () const
  {
    return N;
  }

} // end of namespace trajectory.
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HXX
