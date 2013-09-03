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
# include <roboptim/trajectory/vector-interpolation.hh>

namespace roboptim
{
  inline
  VectorInterpolation::VectorInterpolation
  (const vector_t& x, size_type outputSize, value_type dt) throw ()
    : roboptim::Trajectory<3>
      (makeInterval (0., dt * x.size () / outputSize), outputSize, x,
       "vectorInterpolation"),
      dx_ (x.size ()),
      dt_ (dt)
  {
    setParameters (x);
  }

  inline
  VectorInterpolation::~VectorInterpolation () throw ()
  {}

  inline
  void
  VectorInterpolation::impl_compute
  (result_t& result, double t)
    const throw ()
  {
    size_type before = static_cast<size_type> (std::floor (t / dt_));
    size_type after = static_cast<size_type> (std::ceil (t / dt_));

    // If we are outsize the domain of definition, return zero.
    if (before < 0)
      return;
    if (after > parameters ().size ())
      return;

    result = parameters ().segment (before * this->outputSize (), this->outputSize ());
    // If we are exactly on a data point, just return.
    if (before == after)
      return;

    // Otherwise interpolate.
    // f(x) = (1. - alpha) * x_before + alpha * x_after
    double alpha = static_cast<value_type> (after) - (t / dt_);
    result *= 1. - alpha;
    result +=
      alpha * parameters ().segment (after * this->outputSize (), this->outputSize ());
  }

  inline
  void
  VectorInterpolation::setParameters (const vector_t& x) throw ()
  {
    if (x.size () % this->outputSize () != 0)
      throw std::runtime_error ("x and outputSize are not compatible");

    Trajectory<3>::setParameters (x);

    //FIXME: this should be configurable.
    for (size_type i = 0; i < x.size () - 1; ++i)
      dx_[i] = (x[i] - x[i + 1]) / (2. * dt_);
  }

  inline
  VectorInterpolation::jacobian_t
  VectorInterpolation::variationConfigWrtParam (double t) const throw ()
  {
    throw std::runtime_error ("NOT IMPLEMENTED");
  }

  inline
  VectorInterpolation::jacobian_t
  VectorInterpolation::variationDerivWrtParam (double t, size_type order)
    const throw ()
  {
    throw std::runtime_error ("NOT IMPLEMENTED");
  }

  inline
  VectorInterpolation::value_type
  VectorInterpolation::singularPointAtRank (size_type) const
  {
    return value_type (); // zero
  }

  inline
  VectorInterpolation::vector_t
  VectorInterpolation::derivBeforeSingularPoint (size_type rank, size_type order) const
  {}

  inline
  VectorInterpolation::vector_t
  VectorInterpolation::derivAfterSingularPoint (size_type rank, size_type order) const
  {}

  inline
  VectorInterpolation::jacobian_t
  VectorInterpolation::variationConfigWrtParam (StableTimePoint tp) const throw ()
  {
    throw std::runtime_error ("NOT IMPLEMENTED");
  }

  inline
  VectorInterpolation::jacobian_t
  VectorInterpolation::variationDerivWrtParam
  (StableTimePoint tp, size_type order) const throw ()
  {
    throw std::runtime_error ("NOT IMPLEMENTED");
  }

  inline
  void
  VectorInterpolation::impl_derivative (gradient_t& gradient,
					double t,
					size_type order)
    const throw ()
  {
    size_type before = static_cast<size_type> (std::floor (t / dt_));
    size_type after = static_cast<size_type> (std::ceil (t / dt_));

    // If we are outsize the domain of definition, return zero.
    if (before < 0)
      return;
    if (after > parameters ().size ())
      return;

    if (order == 0)
      this->operator () (gradient, t);
    else if (order == 1)
      {
	gradient =
	  dx_.segment (before * this->outputSize (), this->outputSize ());

	// If we are exactly on a data point, just return.
	if (before == after)
	  return;

	// Otherwise interpolate.
	// f(x) = (1. - alpha) * x_before + alpha * x_after
	double alpha = static_cast<value_type> (after) - (t  / dt_);
	gradient *= 1. - alpha;
	gradient +=
	  alpha * dx_.segment (after * this->outputSize (),
			       this->outputSize ());
      }
    else if (order == 2)
      throw std::runtime_error ("NOT IMPLEMENTED");
  }

  inline
  void
  VectorInterpolation::impl_derivative (gradient_t& gradient,
					StableTimePoint t,
					size_type functionId)
    const throw ()
  {
    throw std::runtime_error ("NOT IMPLEMENTED");
  }

  Trajectory<3>*
  VectorInterpolation::resize (interval_t timeRange)
    const throw ()
  {
    throw std::runtime_error ("NOT IMPLEMENTED");
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_FILTER_VECTOR_INTERPOLATION_HXX
