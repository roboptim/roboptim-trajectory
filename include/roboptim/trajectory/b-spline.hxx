// Copyright (C) 2010 by Thomas Moulard, AIST, CNRS, INRIA.

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
  //FIXME: defined_lc_in has to be true (false untested).
  template <int N>
  BSpline<N>::BSpline (interval_t tr, size_type outputSize,
		       const vector_t& p,
		       std::string name)
    throw ()
    : Trajectory<N> (tr, outputSize, p, name),
      nbp_ (p.size () / outputSize), uniform_ (true)
  {
    //Parameter size should be a multiple of spline dimension
    assert (this->parameters_.size () % outputSize == 0);
    // number of control points should be at least order + 1 so that
    // there is at least one segment in spline
    assert (nbp_ >= order_ + 1);
    // Fill vector of regularly spaced knots
    size_type m = nbp_ + order_ + 1;

    //FIXME: check this
    double delta_t =
      (tr.second - tr.first)
      / static_cast<double> (m - order_ - 1 - order_);

    double ti = tr.first - order_*delta_t;
    knots_.resize(m);
    for (size_type i=0; i<m; i++)
      {
	knots_(i) = ti;
	ti += delta_t;
      }
    setParameters (p);
    computeBasisPolynomials ();
  }

  template <int N>
  BSpline<N>::BSpline (interval_t tr, size_type outputSize,
		       const vector_t& p, const vector_t& knots,
		       std::string name)
    throw ()
    : Trajectory<N> (tr, outputSize, p, name),
      nbp_ (p.size () / outputSize), knots_(knots), uniform_ (false)
  {
    //Parameter size should be a multiple of spline dimension
    assert (this->parameters_.size () % outputSize == 0);
    // number of control points should be at least 6.
    assert (nbp_ >= order_ + 1);
    // calculated number of control points must match the recieved
    // control point set
    assert (nbp_ + order_  + 1 == knots.size());
    // control points must be a monotonically increasing series
    for(size_type idx = 0; idx < knots.size () - 1; idx++)
      assert (knots[idx] <= knots[idx+1]);

    // not more then the order_ of knot points with the same value are allowed

    setParameters (p);
    computeBasisPolynomials ();
  }

  template <int N>
  BSpline<N>::BSpline (const BSpline<N>& spline) throw ()
    : Trajectory<N> (spline.timeRange (), spline.outputSize (),
		     spline.parameters_),
      nbp_ (spline.parameters_.size () / spline.outputSize ()),
      knots_ (spline.knots_), basisPolynomials_ (spline.basisPolynomials_),
      uniform_ (spline.uniform_)
  {
    // Parameter size should be a multiple of spline dimension
    assert (this->parameters_.size () % this->outputSize () == 0);
    // number of control points should be at least 6.
    assert (nbp_ >= order_ + 1);

    setParameters (spline.parameters ());
  }

  template <int N>
  Trajectory<N>* BSpline<N>::resize (interval_t timeRange) const throw ()
  {
    return new BSpline<N> (timeRange,
			   this->outputSize (),
			   this->parameters (),
			   this->knots_,
			   this->getName());
  }


  /**
   * \brief Generate base polynomial set
   * For the basic spline formula noted in the pdf from the docs section,
   * this implements (3), the recursion formula. It returns the results as
   * factors of b_j,0 where j is the index of the returned map.
   * \param j current knot index
   * \param n current basis function order
   * \return a std::map.
   * Call this with j as the index of the knot for which spline segment your
   * want to generate the basis function for, and n as the order of your spline.
   */
  template <int N>
  typename BSpline<N>::cox_map
  BSpline<N>::cox_de_boor(size_type j, size_type n) const
  {
    const bool debug=false;
    std::stringstream label;
    label << "[j=" << j << "/n=" << n << "]";

    if(n == 0) //end of recursion
      {
	cox_map map;
	Eigen::Matrix<double,N+1,1> temp_params;
	temp_params.setZero();
	temp_params[0]=1.;
	std::pair<cox_map_itr_t,bool> ptr =
	  map.insert (std::make_pair (j, polynomial_t (0., temp_params)));

	if(debug)
	  std::cout
	    << label.str() << "return(end)    : "
	    << ptr.first->second << std::endl;
	return map;
      }
    else
      {
	// t_{j}
	const double t0 = knots_[j+0];
	// t_{j+1}
	const double t1 = knots_[j+1];
	// t_{j+n}
	const double tn = knots_[j+n];
	// t_{j+n+1}
	const double tn1 = knots_[j+n+1];

	polynomial_t p_1_rat = 1./(tn-t0) * monomial_t(t0);
	// http://wolftype.com/ucsb/spatial/bspline.html, uniform b-splines
	if (std::isinf (p_1_rat.coefs_[1]))
	  {
	    //FIXME: this is probably not how it should work.
	    p_1_rat=monomial_t(t0);
	  }

	if(debug)
	  std::cout << label.str()
		    << "p_1_rat        : " << p_1_rat << std::endl;

	cox_map p_1_cox = cox_de_boor(j,n-1);
	cox_map p_1;

	for (cox_map_itr_t itr = p_1_cox.begin (); itr != p_1_cox.end(); itr++)
	  {
	    if(debug)
	      std::cout << label.str()
			<< "p_1_cox        : " << itr->first
			<< " : " << itr->second << std::endl;
	  polynomial_t p_prod = itr->second * p_1_rat;
	  if(debug)
	    std::cout << label.str() << "p_prod (1)     : "
		      << p_prod << std::endl;
	  p_1.insert (std::make_pair (itr->first,p_prod));

	}

	polynomial_t p_2_rat = 1. / (t1 - tn1) * monomial_t (tn1);

	// http://wolftype.com/ucsb/spatial/bspline.html, uniform b-splines
	if (std::isinf(p_2_rat.coefs_[1]))
	  {
	    //FIXME: this is probably not how it should work.
	    p_2_rat = Monomial<N>(tn1);
	  }

	if(debug)
	  std::cout << label.str() << "p_2_rat        : "
		    << p_2_rat << std::endl;

	cox_map p_2_cox = cox_de_boor(j+1,n-1);
	cox_map p_2;
	for(cox_map_itr_t itr=p_2_cox.begin();itr!=p_2_cox.end();itr++)
	  {
	    if(debug)
	      std::cout << label.str()
			<< "p_2_cox        : " << itr->first
			<< " : " << itr->second << std::endl;
	    polynomial_t p_prod = itr->second * p_2_rat;
	    if(debug)std::cout << label.str()
			       << "p_prod (2)     : " << p_prod << std::endl;
	    p_2.insert(std::make_pair(itr->first,p_prod));
	  }

	cox_map p;
	p.insert(p_1.begin(),p_1.end());
	for(cox_map_itr_t itr=p_2.begin();itr!=p_2.end();itr++)
	  {
	    cox_map_itr_t existing_element = p.find(itr->first);
	    if(existing_element!=p.end())
		existing_element->second = existing_element->second + itr->second;
	    else
	      p.insert(std::make_pair(itr->first,itr->second));
	  }

	if(debug)
	  std::cout << label.str() << "result of recursion branch" << std::endl;
	for(cox_map_itr_t itr=p.begin();itr!=p.end();itr++)
	  {
	    if(debug)
	      std::cout << label.str() << "B_" << j
			<< "_" << itr->first << "  : "
			<< itr->second << std::endl;
	  }
	return p;
      }
  }

  template <int N>
  void BSpline<N>::computeBasisPolynomials ()
  {
    const bool debug=false;
    basisPolynomials_.clear();
    for (size_type j=0; j<nbp_; j++)
      {
	//calculate basis polynomials for each interval of the knot vector
	basisPolynomials_.push_back (std::vector <polynomial_t> ());
	cox_map map = cox_de_boor (j, order_);

	for(cox_map_itr_t itr = map.begin(); itr!=map.end(); itr++)
	  {
	    if (debug)
	      std::cout << "B_" << j << "_"
			<< itr->first << "  : " << itr->second << std::endl;
	    basisPolynomials_.back ().push_back (itr->second);
	  }
      }
  }

  template <int N>
  void BSpline<N>::setParameters (const vector_t& p) throw ()
  {
    assert (p.size () == this->parameters_.size ());
    this->parameters_ = p;
  }

  template <int N>
  void BSpline<N>::impl_compute (result_t& derivative, double t) const throw ()
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
    return this->length () / ((value_type)nbp_ - static_cast<double> (order_));
  }

  template <int N>
  typename BSpline<N>::size_type
  BSpline<N>::interval (value_type t) const
  {
    t = detail::fixTime (t, *this);
    typedef boost::numeric::converter<size_type, double> Double2SizeType;
    // t_{order]
    double tmin = this->timeRange ().first;
    size_type imin = order_;
    // t_{m-?}
    double tmax = this->timeRange ().second;
    size_type imax = nbp_;

    unsigned int count = 0;
    bool found = false;
    size_type i = 1;
    size_type iPrev = 0;
    while (!found && iPrev != i)
      {
	i = Double2SizeType::convert
	  (std::floor (imin + (t - tmin) / (tmax - tmin) * (imax - imin)));
	if (t < knots_ [i])
	  {
	    tmax = knots_ [i-1];
	    imax = i-1;
	  }
	else if (t >= knots_ [i+1])
	  {
	    if (t < knots_ [i+2])
	      {
		i = i+1;
		found = true;
	      }
	    imin = i+1;
	    tmin = knots_ [i+1];
	  }
	else
	  {
	    found = true;
	  }
	count++;
	assert (count < 10000);
	iPrev = i;
      }
    if (i > nbp_-1)
      i = nbp_-1;
    if (i < order_)
      i = order_;
    return i;
  }

  template <int N>
  typename BSpline<N>::vector_t
  BSpline<N>::basisFunctions (value_type t, size_type order) const
  {
    vector_t result (order_+1);

    t = detail::fixTime (t, *this);
    const size_type k = interval (t);
    const size_type n = this->outputSize ();

    for(size_type idx=0;idx<order_+1;idx++)
      {
	const Polynomial<N>& B = basisPolynomials_[k-idx][idx];
	result(idx) =  B.derivative(t,order);
      }
  }

  template <int N> void
  BSpline<N>::impl_derivative (gradient_t& derivative, double t,
			       size_type order)
    const throw ()
  {

    t = detail::fixTime (t, *this);
    const size_type k = interval (t);
    const size_type n = this->outputSize ();

    derivative.setZero();

# ifndef NDEBUG
    double polynomial_sum = 0.;
# endif //! NDEBUG
    for(size_type idx=0;idx<order_+1;idx++)
      {
	const vector_t& P_seg = this->parameters_.segment((k - idx) * n, n);
	const Polynomial<N>& B = basisPolynomials_[k-idx][idx];

# ifndef NDEBUG
	polynomial_sum += B.derivative(t,order);
# endif //! NDEBUG

	derivative +=  B.derivative(t,order) * P_seg;
    }

    /* this is true for any knot vector (FIXME: reference) */
    if(order == 0)
      {
	assert
	  (std::abs( polynomial_sum - 1.)
	   < std::numeric_limits<double>::epsilon() * 1e6);
      }
  }

  template <int N> void
  BSpline<N>::impl_derivative (gradient_t& derivative,
			       StableTimePoint stp,
			       size_type order) const throw ()
  {
    this->impl_derivative (derivative,
			   stp.getTime (this->timeRange ()),
			   order);
  }

  template <int N>
  typename BSpline<N>::jacobian_t
  BSpline<N>::variationConfigWrtParam (double t) const throw ()
  {
    return variationDerivWrtParam (t, 0);
  }


  template <int N>
  typename BSpline<N>::jacobian_t
  BSpline<N>::variationDerivWrtParam (double t, size_type order)
    const throw ()
  {

    t = detail::fixTime (t, *this);
    const size_type k = interval (t);
    const size_type n = this->outputSize ();

    jacobian_t jac(n, nbp_ * n);
    jac.setZero();
    matrix_t In(n,n);
    In.setIdentity();

    for(size_type idx=0;idx<order_+1;idx++)
      {
	const Polynomial<N>& B = basisPolynomials_[k-idx][idx];

	jac.middleCols((k - order_) * n, n) =
	  B.derivative(t, order) * In;
      }

    return jac;
  }

  template <int N>
  typename BSpline<N>::jacobian_t
  BSpline<N>::variationConfigWrtParam (StableTimePoint stp)
    const throw ()
  {
    return this->variationConfigWrtParam (stp.getTime (this->timeRange ()));
  }


  template <int N>
  typename BSpline<N>::jacobian_t
  BSpline<N>::variationDerivWrtParam (StableTimePoint stp, size_type order)
    const throw ()
  {
    return this->variationDerivWrtParam
      (stp.getTime (this->timeRange ()), order);
  }

  template <int N>
  typename BSpline<N>::value_type
  BSpline<N>::singularPointAtRank (size_type rank) const
  {
    return (value_type)rank * this->length () / ((value_type)nbp_- order_);
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
  BSpline<N>::print (std::ostream& o) const throw ()
  {
    o << "order " << order_ << " BSpline" << incindent << std::endl
      << "Number of parameters per spline function: " << nbp_ << std::endl
      << "Length: " << this->length () << std::endl
      << "control points: ";
    for(int idx=0;idx<knots_.size();idx++)o << knots_[idx] << " ";
    o << "Parameters: " << this->parameters ()
      << decindent;
    return o;
  }

  template <int N>
  typename BSpline<N>::vector_t const & BSpline<N>::knots() const
  {
    return this->knots_;
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HXX
