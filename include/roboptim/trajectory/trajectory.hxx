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

namespace roboptim
{
namespace trajectory
{
  template <unsigned dorder>
  Trajectory<dorder>::Trajectory (interval_t tr,
				  size_type outputSize,
				  const vector_t& p,
				  std::string name)
    : parent_t (outputSize, name),
      timeRange_ (tr),
      parameters_ (p),
      singularPoints_ (),
      tolerance_ (1e-5)
  {
    //FIXME: can a trajectory be a single point?
    assert (tr.first <= tr.second);
  }


  template <unsigned dorder>
  Trajectory<dorder>::~Trajectory ()
  {
  }


  template <unsigned dorder>
  const typename Trajectory<dorder>::vector_t&
  Trajectory<dorder>::parameters () const
  {
    return parameters_;
  }


  template <unsigned dorder>
  void
  Trajectory<dorder>::setParameters (const vector_t& p)
  {
    parameters_ = p;
  }


  template <unsigned dorder>
  typename Trajectory<dorder>::interval_t
  Trajectory<dorder>::timeRange () const
  {
    return timeRange_;
  }

  template <unsigned dorder>
  typename Trajectory<dorder>::value_type
  Trajectory<dorder>::length () const
  {
    return timeRange ().second - timeRange ().first;
  }


  template <unsigned dorder>
  typename Trajectory<dorder>::vector_t
  Trajectory<dorder>::state (double t, size_type order) const
  {
    const size_type dimension = this->outputSize ();
    vector_t result ((order + 1) * dimension);

    for (size_type o = 0; o <= order; ++o)
      result.segment(o * dimension,dimension) =
	this->derivative (t, o);
    return result;
  }

  template <unsigned dorder>
  typename Trajectory<dorder>::vector_t
  Trajectory<dorder>::state
  (StableTimePoint stp, size_type order) const
  {
    const size_type dimension = this->outputSize ();
    vector_t result ((order + 1) * dimension);

    for (size_type o = 0; o <= order; ++o)
      result.segment(o * dimension, dimension) =
	this->derivative (stp, o);
    return result;
  }


  template <unsigned dorder>
  typename Trajectory<dorder>::jacobian_t
  Trajectory<dorder>::variationStateWrtParam (double t, size_type order)
    const
  {
    const size_type dimension = this->outputSize ();
    const size_type parameterSize = parameters ().size ();
    jacobian_t result (dimension * (order + 1), parameterSize);

    for (size_type o = 0; o <= order; ++o)
      result.block (o * dimension, 0, dimension, parameterSize)
	= this->variationDerivWrtParam (t, o);
    return result;
  }

  template <unsigned dorder>
  typename Trajectory<dorder>::jacobian_t
  Trajectory<dorder>::variationStateWrtParam
  (StableTimePoint stp, size_type order)
    const
  {
    const size_type dimension = this->outputSize ();
    const size_type parameterSize = parameters ().size ();
    jacobian_t result (dimension * (order + 1), parameterSize);

    for (size_type o = 0; o <= order; ++o)
      result.block (o * dimension, 0, dimension, parameterSize)
	= this->variationDerivWrtParam (stp, o);
    return result;
  }


  template <unsigned dorder>
  typename Trajectory<dorder>::size_type
  Trajectory<dorder>::singularPoints () const
  {
    return singularPoints_;
  }


  template <unsigned dorder>
  void
  Trajectory<dorder>::impl_compute (result_t& res , StableTimePoint stp)
    const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
      Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    (*this) (res, stp.getTime (this->timeRange ()));
  }

  template <unsigned dorder>
  bool
  Trajectory<dorder>::isValidTime (value_type t) const
  {
    value_type tmin = this->getLowerBound (this->timeRange ());
    value_type tmax = this->getUpperBound (this->timeRange ());
    return (tmin <= t && t <= tmax);
  }

  template <unsigned dorder>
  void
  Trajectory<dorder>::normalizeAngles (size_type index)
  {
    this->normalizeAngles (index, 0);
  }

  template <unsigned dorder>
  void
  Trajectory<dorder>::normalizeAngles (size_type index, size_type offset)
  {
    value_type thetaPrev = 0.;
    for (int i = 0; i < (parameters_.size () - offset) / this->outputSize (); ++i)
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
  void Trajectory<dorder>::tolerance (const double& tolerance)
  {
    tolerance_ = tolerance;
  }

  template <unsigned dorder>
  double Trajectory<dorder>::tolerance () const
  {
    return tolerance_;
  }

  template <unsigned dorder>
  std::ostream&
  Trajectory<dorder>::print (std::ostream& o) const
  {
    o << "Generic (abstract) trajectory." << std::endl;
    return o;
  }
} // end of namespace trajectory.
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HXX
