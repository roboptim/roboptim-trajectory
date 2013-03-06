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

#ifndef ROBOPTIM_TRAJECTORY_FREETIMETRAJECTORY_HXX
# define ROBOPTIM_TRAJECTORY_FREETIMETRAJECTORY_HXX
# include <cmath>

# include <boost/numeric/ublas/matrix.hpp>
# include <boost/numeric/ublas/matrix_proxy.hpp>

namespace roboptim
{
  namespace detail
  {
    inline Function::value_type
    unscaleTime (Function::value_type unscaled,
	       Function::value_type min,
	       Function::value_type scale)
    {
      return min + (unscaled - min) / scale;
    }

    template <unsigned dorder>
    Function::interval_t
    unscaleInterval (const Trajectory<dorder>& traj,
		     typename Function::value_type scale)
    {
      Function::value_type min =
	Function::getLowerBound (traj.timeRange ());
      Function::value_type max =
	unscaleTime (Function::getUpperBound (traj.timeRange ()), min, scale);
      return Function::makeInterval (min, max);
    }
  } // end of namespace detail

  template <typename T>
  FreeTimeTrajectory<T>::FreeTimeTrajectory
  (const fixedTimeTrajectory_t& traj, value_type s) throw ()
    : parent_t (detail::unscaleInterval (traj, s), traj.outputSize (),
		addScaleToParameters (traj.parameters (), s)),
      trajectory_ (traj.clone ())
  {
    assert (s != 0. && !std::isinf (s) && !std::isnan (s));
  }

  template <typename T>
  FreeTimeTrajectory<T>::FreeTimeTrajectory
  (const self_t& traj)
    throw ()
    : parent_t (traj.timeRange (), traj.outputSize (),
		traj.parameters ()),
      trajectory_ (traj.trajectory_->clone ())
  {
  }

  template <typename T>
  FreeTimeTrajectory<T>::~FreeTimeTrajectory () throw ()
  {
    delete trajectory_;
  }

  template <typename T>
  void
  FreeTimeTrajectory<T>::impl_compute (result_t& res , double t)
    const throw ()
  {
    (*trajectory_) (res, this->scaleTime (t));
  }

  template <typename T>
  void
  FreeTimeTrajectory<T>::impl_derivative (gradient_t& derivative,
					  double t,
					  size_type order) const throw ()
  {
    double scaled = this->scaleTime (t);
    trajectory_->derivative (derivative, scaled, order);
    derivative *= std::pow (this->timeScale (), 0. + order);
  }

  template <typename T>
  void
  FreeTimeTrajectory<T>::impl_derivative (gradient_t& derivative,
					       StableTimePoint stp,
					       size_type order) const throw ()
  {
    trajectory_->derivative (derivative, stp, order);
    derivative *= std::pow (this->timeScale (), 0. + order);
  }

  template <typename T>
  typename FreeTimeTrajectory<T>::jacobian_t
  FreeTimeTrajectory<T>::variationConfigWrtParam (double t) const throw ()
  {
    using namespace boost::numeric::ublas;
    value_type scaled = this->scaleTime (t);

    double tMin = this->getLowerBound (this->trajectory_->timeRange ());

    jacobian_t result (this->outputSize (),
		       this->parameters ().size ());
    result.setZero();

    // Compute variation w.r.t time scale (p_0)
    result.leftCols(1) = (t - tMin) * trajectory_->derivative (scaled, 1);

    // Fill 1..(n-1) lines with original jacobian.
    result.rightCols(result.cols()-1)
      = trajectory_->variationConfigWrtParam (scaled);

    return result;
  }

  template <typename T>
  typename FreeTimeTrajectory<T>::jacobian_t
  FreeTimeTrajectory<T>::variationDerivWrtParam (double t, size_type order)
    const throw ()
  {
    if (order == 0)
      return this->variationConfigWrtParam (t);

    using namespace boost::numeric::ublas;
    value_type scaled = this->scaleTime (t);
    value_type tmin = this->getLowerBound (this->trajectory_->timeRange ());

    jacobian_t result (this->outputSize (),
		       this->parameters ().size ());
    result.setZero();

    // Compute variation w.r.t time scale (p_0)
    result.leftCols(1) = trajectory_->derivative (scaled, order) * order;
    result.leftCols(1)
      += trajectory_->derivative (scaled, order + 1)
      * (this->timeScale () * (t - tmin));
    result.leftCols(1) *= std::pow (this->timeScale (), order - 1.);

    // Fill 1..(n-1) lines with original jacobian.
    result.rightCols(result.cols()-1)
      = trajectory_->variationDerivWrtParam(scaled, order)
      * std::pow (this->timeScale (), 0. + order);

    return result;
  }

  template <typename T>
  typename FreeTimeTrajectory<T>::jacobian_t
  FreeTimeTrajectory<T>::variationConfigWrtParam
  (StableTimePoint stp) const throw ()
  {
    using namespace boost::numeric::ublas;

    jacobian_t result (this->outputSize (),
		       this->parameters ().size ());
    result.setZero();

    // Compute variation w.r.t time scale (p_0)
    // ...is null do not do anything.

    // Fill 1..(n-1) lines with original jacobian.
    result.rightCols(result.cols()-1)
      = trajectory_->variationConfigWrtParam (stp);

    return result;
  }


  template <typename T>
  typename FreeTimeTrajectory<T>::jacobian_t
  FreeTimeTrajectory<T>::variationDerivWrtParam
  (StableTimePoint stp, size_type order) const throw ()
  {
    if (order == 0)
      return this->variationConfigWrtParam (stp);

    using namespace boost::numeric::ublas;

    jacobian_t result (this->outputSize (),
		       this->parameters ().size ());
    result.setZero();

    // Compute variation w.r.t time scale (p_0)
    // ...is null do not do anything.

    // Fill 1..(n-1) lines with original jacobian.
    result.rightCols(result.cols()-1)
      = trajectory_->variationDerivWrtParam (stp, order)
      * (order + 1) * (order + 1);

    return result;
  }

  template <typename T>
  typename FreeTimeTrajectory<T>::value_type
  FreeTimeTrajectory<T>::singularPointAtRank (size_type rank) const
  {
    double tMin = this->getLowerBound (this->timeRange ());
    return tMin + (trajectory_->singularPointAtRank (rank) - tMin)
      * this->timeScale ();
  }

  template <typename T>
  typename FreeTimeTrajectory<T>::vector_t
  FreeTimeTrajectory<T>::derivBeforeSingularPoint (size_type rank,
							size_type order)
    const
  {
    return trajectory_->derivBeforeSingularPoint (rank, order);
  }

  template <typename T>
  typename FreeTimeTrajectory<T>::vector_t
  FreeTimeTrajectory<T>::derivAfterSingularPoint (size_type rank, size_type order)
    const
  {
    return trajectory_->derivAfterSingularPoint (rank, order);
  }

  template <typename T>
  void
  FreeTimeTrajectory<T>::setParameters (const vector_t& p) throw ()
  {
    //FIXME: is this ok?
    vector_t p_ = p;
    if (p_[0] <= 0.)
      p_[0] = 1e-8;

    this->parameters_ = p_;
    this->timeRange_ =
      detail::unscaleInterval (*trajectory_, this->timeScale ());
    this->trajectory_->setParameters (removeScaleFromParameters (p_));
  }

  template <typename T>
  typename FreeTimeTrajectory<T>::value_type
  FreeTimeTrajectory<T>::timeScale () const throw ()
  {
    return this->parameters_[0];
  }

  template <typename T>
  double
  FreeTimeTrajectory<T>::scaleTime (double unscaled) const throw ()
  {
    value_type tMin = this->getLowerBound (this->timeRange ());
    value_type tmin = this->getLowerBound (this->trajectory_->timeRange ());
    value_type tmax = this->getUpperBound (this->trajectory_->timeRange ());

    unscaled = detail::fixTime (unscaled, *this);
    assert (this->isValidTime (unscaled));

    value_type res = tmin + (unscaled - tMin) * timeScale ();

    if (res > tmax)
      res = tmax;
    else if (res < tmin)
      res = tmin;
    return res;
  }

  template <typename T>
  double
  FreeTimeTrajectory<T>::unscaleTime (double scaled) const throw ()
  {
    value_type tMin = this->getLowerBound (this->timeRange ());
    value_type tMax = this->getUpperBound (this->timeRange ());
    value_type tmin = this->getLowerBound (this->trajectory_->timeRange ());

    scaled = detail::fixTime (scaled, *this);
    assert (trajectory_->isValidTime (scaled));

    value_type res = tMin + (scaled - tmin) / timeScale ();

    if (res > tMax)
      res = tMax;
    else if (res < tMin)
      res = tMin;
    return res;
  }

  template <typename T>
  std::ostream&
  FreeTimeTrajectory<T>::print (std::ostream& o) const throw ()
  {
    o << "Free time trajectory." << std::endl;
    return o;
  }

  template <typename T>
  void
  FreeTimeTrajectory<T>::normalizeAngles (size_type index) throw ()
  {
    this->normalizeAngles (index, 1);
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_FREETIMETRAJECTORY_HXX
