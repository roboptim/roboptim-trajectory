// Copyright (C) 2013 by Thomas Moulard, AIST, CNRS, INRIA.
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

#ifndef ROBOPTIM_TRAJECTORY_FILTER_VECTOR_INTERPOLATION_HXX
# define ROBOPTIM_TRAJECTORY_FILTER_VECTOR_INTERPOLATION_HXX
# include <boost/format.hpp>
# include <roboptim/trajectory/vector-interpolation.hh>

namespace roboptim
{
namespace trajectory
{
  inline
  VectorInterpolation::VectorInterpolation
  (const_vector_ref x, size_type outputSize, value_type dt)
    : roboptim::trajectory::Trajectory<3>
      (makeInterval (0., dt * static_cast<value_type> (x.size ())
                     / static_cast<value_type> (outputSize)),
       outputSize, x,
       "vectorInterpolation"),
      dt_ (dt)
  {
    setParameters (x);
  }

  inline
  VectorInterpolation::~VectorInterpolation ()
  {}

  inline
  VectorInterpolation::VectorInterpolation (const VectorInterpolation& vi)
    : roboptim::trajectory::Trajectory<3> (vi),
      dt_ (vi.dt_)
  {
    setParameters (vi.parameters ());
  }

  inline
  VectorInterpolation::size_type
  VectorInterpolation::numFrames () const
  {
    return
      static_cast<size_type> (this->parameters ().size ())
      / this->outputSize ();
  }

  inline boost::shared_ptr<VectorInterpolation>
  VectorInterpolation::trim (size_type start, size_type length) const
  {
    // check that start is not too high
    if (start >= this->numFrames ())
      throw std::runtime_error ("invalid starting frame");
    // if length is zero, compute automatically how many frames are
    // remaining
    if (length <= 0)
      length = this->numFrames () - start;
    // check again that the arguments are correct
    if (length < 1 || start + length > this->numFrames ())
      throw std::runtime_error ("invalid length");

    Eigen::VectorBlock<const vector_t> x = this->parameters ().segment
      (start * this->outputSize (), length * this->outputSize ());

    boost::shared_ptr<VectorInterpolation> result =
      boost::make_shared<VectorInterpolation> (x, this->outputSize (), this->dt_);
    return result;
  }

  inline
  void
  VectorInterpolation::impl_compute
  (result_ref result, double t)
    const
  {
    // Fix evaluation if we only have one frame
    if (parameters ().size () == this->outputSize ())
      {
	result = parameters ();
	return;
      }

    size_type before = static_cast<size_type> (std::floor (t / dt_));
    size_type after = static_cast<size_type> (std::ceil (t / dt_));

    // If we are outsize the domain of definition, return zero.
    if (before < 0)
      return;
    if (after >= (parameters ().size () / this->outputSize ()))
      return;

    result = parameters ().segment (before * this->outputSize (), this->outputSize ());
    // If we are exactly on a data point, just return.
    if (before == after)
      return;

    // Otherwise interpolate.
    // f(x) = (1. - alpha) * x_before + alpha * x_after
    double alpha = ((t / dt_) - std::floor (t / dt_)) * dt_;
    result *= 1. - alpha;
    result +=
      alpha * parameters ().segment (after * this->outputSize (), this->outputSize ());
  }

  inline
  void
  VectorInterpolation::setParameters (const vector_t& x)
  {
    if (x.size () % this->outputSize () != 0)
      {
	boost::format fmt
	  ("parameter vector size (%d) does not match outputSize (%d)");
	fmt % x.size () % this->outputSize ();
	throw std::runtime_error (fmt.str ());
      }

    Trajectory<3>::setParameters (x);
    singularPoints_ = this->parameters ().size () / this->outputSize ();
  }

  inline
  VectorInterpolation::jacobian_t
  VectorInterpolation::variationConfigWrtParam (double t) const
  {
    return variationDerivWrtParam (t, 0);
  }

  inline
  VectorInterpolation::jacobian_t
  VectorInterpolation::variationDerivWrtParam (double t, size_type order)
    const
  {
    VectorInterpolation::jacobian_t jacobian (this->outputSize (),
					      this->parameters ().size ());
    jacobian.setZero ();

    size_type before = static_cast<size_type> (std::floor (t / dt_));
    size_type after = static_cast<size_type> (std::ceil (t / dt_));

    // If we are outsize the domain of definition, return zero.
    if (before < 0)
      return jacobian;
    if (after >= (parameters ().size () / this->outputSize ()))
      return jacobian;

    if (order == 0)
      {
	double alpha = ((t / dt_) - std::floor (t / dt_)) * dt_;
	jacobian.block
	  (0, before * this->outputSize (),
	   this->outputSize (),
	   this->outputSize ()).diagonal ().setConstant (1 - alpha);
	jacobian.block
	  (0, after * this->outputSize (),
	   this->outputSize (),
	   this->outputSize ()).diagonal ().setConstant (alpha);
      }
    else if (order == 1)
      {
	jacobian.block
	  (0, before * this->outputSize (),
	   this->outputSize (),
	   this->outputSize ()).diagonal ().setConstant (- 1. / dt_);
	jacobian.block
	  (0, after * this->outputSize (),
	   this->outputSize (),
	   this->outputSize ()).diagonal ().setConstant (1. / dt_);
      }
    else
      jacobian.setZero ();

    return jacobian;
  }

  inline
  VectorInterpolation::value_type
  VectorInterpolation::singularPointAtRank (size_type rank) const
  {
    return static_cast<value_type> (rank) * dt_;
  }

  inline
  VectorInterpolation::vector_t
  VectorInterpolation::derivBeforeSingularPoint
  (size_type rank, size_type order) const
  {
    VectorInterpolation::vector_t gradient (this->gradientSize ());
    VectorInterpolation::value_type t =
      this->singularPointAtRank (rank) - (dt_ / 2.);
    this->derivative (gradient, t, order);
    return gradient;
  }

  inline
  VectorInterpolation::vector_t
  VectorInterpolation::derivAfterSingularPoint
  (size_type rank, size_type order) const
  {
    VectorInterpolation::vector_t gradient (this->gradientSize ());
    VectorInterpolation::value_type t =
      this->singularPointAtRank (rank) + (dt_ / 2.);
    this->derivative (gradient, t, order);
    return gradient;
  }

  inline
  VectorInterpolation::jacobian_t
  VectorInterpolation::variationConfigWrtParam (StableTimePoint tp) const
  {
    return variationDerivWrtParam (tp, 0);
  }

  inline
  VectorInterpolation::jacobian_t
  VectorInterpolation::variationDerivWrtParam
  (StableTimePoint stp, size_type order) const
  {
    return this->variationDerivWrtParam
      (stp.getTime (this->timeRange ()), order);
  }

  inline
  void
  VectorInterpolation::impl_derivative (gradient_ref gradient,
					double t,
					size_type order)
    const
  {
    size_type before = static_cast<size_type> (std::floor (t / dt_));
    size_type after = static_cast<size_type> (std::ceil (t / dt_));

    // If we are outsize the domain of definition, return zero.
    if (before < 0)
      return;
    if (after >= (parameters ().size () / this->outputSize ()))
      return;

    gradient.setZero ();

    if (order == 0)
      this->operator () (gradient, t);
    else if (order == 1)
      {
	// If we are exactly on a data point, just return.
	if (before == after)
	  return;

	gradient =
	  this->parameters ().segment (after * outputSize (), outputSize ())
	  - this->parameters ().segment (before * outputSize (), outputSize ());
	gradient /= dt_;
      }
    else
      gradient.setZero ();
  }

  inline
  void
  VectorInterpolation::impl_derivative (gradient_ref derivative,
					StableTimePoint stp,
					size_type order)
    const
  {
    this->impl_derivative
      (derivative,
       stp.getTime (this->timeRange ()),
       order);
  }

  inline Trajectory<3>*
  VectorInterpolation::resize (interval_t)
    const
  {
    throw std::runtime_error ("NOT IMPLEMENTED");
  }
} // end of namespace trajectory.
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_FILTER_VECTOR_INTERPOLATION_HXX
