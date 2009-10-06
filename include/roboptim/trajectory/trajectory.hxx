// Copyright (C) 2009 by Florent Lamiraux, Thomas Moulard, AIST, CNRS, INRIA.
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

#ifndef ROBOPTIM_TRAJECTORY_TRAJECTORY_HXX
# define ROBOPTIM_TRAJECTORY_TRAJECTORY_HXX
# include <boost/numeric/ublas/matrix_proxy.hpp>
# include <boost/numeric/ublas/vector_proxy.hpp>

namespace roboptim
{
  template <unsigned dorder>
  Trajectory<dorder>::Trajectory (interval_t tr,
				  size_type outputSize,
				  const vector_t& p,
				  std::string name)
    throw ()
    : parent_t (outputSize, name),
      timeRange_ (tr),
      parameters_ (p),
      singularPoints_ ()
  {
    //FIXME: can a trajectory be a single point?
    assert (tr.first <= tr.second);
  }


  template <unsigned dorder>
  Trajectory<dorder>::~Trajectory () throw ()
  {
  }


  template <unsigned dorder>
  const typename Trajectory<dorder>::vector_t&
  Trajectory<dorder>::parameters () const throw ()
  {
    return parameters_;
  }


  template <unsigned dorder>
  void
  Trajectory<dorder>::setParameters (const vector_t& p) throw ()
  {
    parameters_ = p;
  }


  template <unsigned dorder>
  typename Trajectory<dorder>::interval_t
  Trajectory<dorder>::timeRange () const throw ()
  {
    return timeRange_;
  }

  template <unsigned dorder>
  typename Trajectory<dorder>::value_type
  Trajectory<dorder>::length () const throw ()
  {
    return timeRange ().second - timeRange ().first;
  }


  template <unsigned dorder>
  typename Trajectory<dorder>::vector_t
  Trajectory<dorder>::state (double t, size_type order) const throw ()
  {
    using namespace boost::numeric::ublas;
    const value_type dimension = this->outputSize ();
    vector_t result ((order + 1) * dimension);

    for (size_type o = 0; o <= order; ++o)
      subrange (result, o * dimension, (o + 1) * dimension) =
	this->derivative (t, o);
    return result;
  }


  template <unsigned dorder>
  typename Trajectory<dorder>::jacobian_t
  Trajectory<dorder>::variationStateWrtParam (double t, size_type order)
    const throw ()
  {
    using namespace boost::numeric::ublas;
    const size_type dimension = this->outputSize ();
    const size_type parameterSize = parameters ().size ();
    jacobian_t result (dimension * (order + 1), parameterSize);

    for (size_type o = 0; o <= order; ++o)
      {
	range xrange (o * dimension, (o + 1) * dimension);
	range yrange (0, parameterSize);
	project (result, xrange, yrange) = this->variationDerivWrtParam (t, o);
      }
    return result;
  }


  template <unsigned dorder>
  typename Trajectory<dorder>::size_type
  Trajectory<dorder>::singularPoints () const throw ()
  {
    return singularPoints_;
  }


  template <unsigned dorder>
  void
  Trajectory<dorder>::impl_compute
  (typename Trajectory<dorder>::result_t& res , StableTimePoint stp) const throw ()
  {
    (*this) (res, stp.getTime (this->timeRange ()));
  }

  template <unsigned dorder>
  void
  Trajectory<dorder>::impl_derivative
  (typename Trajectory<dorder>::gradient_t& derivative,
   StableTimePoint stp,
   typename Trajectory<dorder>::size_type order) const throw ()
  {
    return this->impl_derivative (derivative,
				  stp.getTime (this->timeRange ()),
				  order);
  }

  template <unsigned dorder>
  typename Trajectory<dorder>::jacobian_t
  Trajectory<dorder>::variationConfigWrtParam (StableTimePoint stp)
    const throw ()
  {
    return this->variationConfigWrtParam (stp.getTime (this->timeRange ()));
  }


  template <unsigned dorder>
  typename Trajectory<dorder>::jacobian_t
  Trajectory<dorder>::variationDerivWrtParam (StableTimePoint stp,
					      size_type order)
    const throw ()
  {
    return this->variationDerivWrtParam
      (stp.getTime (this->timeRange ()), order);
  }

  template <unsigned dorder>
  bool
  Trajectory<dorder>::isValidTime (value_type t) const throw ()
  {
    value_type tmin = getLowerBound (this->timeRange ());
    value_type tmax = getUpperBound (this->timeRange ());
    return (tmin <= t && t <= tmax);
  }

  template <unsigned dorder>
  void
  Trajectory<dorder>::normalizeAngles (size_type index) throw ()
  {
    this->normalizeAngles (index, 0.);
  }

  template <unsigned dorder>
  void
  Trajectory<dorder>::normalizeAngles (size_type index, size_type offset)
    throw ()
  {
    value_type thetaPrev = 0.;
    for (unsigned i = 0; i < (parameters_.size () - offset) / this->outputSize (); ++i)
      {
	value_type& theta = this->parameters_[offset + i * this->outputSize () + index];
	if (theta - thetaPrev > M_PI)
	  theta -= M_PI * 2;
	else if (theta - thetaPrev < -M_PI)
	  theta += M_PI * 2;
	thetaPrev = theta;
      }
  }

  template <unsigned dorder>
  std::ostream&
  Trajectory<dorder>::print (std::ostream& o) const throw ()
  {
    o << "Generic (abstract) trajectory." << std::endl;
    return o;
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HXX
