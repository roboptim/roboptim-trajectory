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
    scaleTime (Function::value_type unscaled,
	       Function::value_type min,
	       Function::value_type scale)
    {
      return min + (unscaled - min) / scale;
    }

    template <unsigned dorder>
    Function::interval_t
    scaleInterval (const Trajectory<dorder>& traj,
		   typename Function::value_type scale)
    {
      Function::value_type min =
	Function::getLowerBound (traj.timeRange ());
      Function::value_type max =
	scaleTime (Function::getUpperBound (traj.timeRange ()), min, scale);
      return Function::makeInterval (min, max);
    }
  } // end of namespace detail

  template <unsigned dorder>
  FreeTimeTrajectory<dorder>::FreeTimeTrajectory
  (const Trajectory<dorder>& traj, value_type s) throw ()
    : Trajectory<dorder> (detail::scaleInterval (traj, s), traj.outputSize (),
			  addScaleToParameters (traj.parameters (), s)),
      trajectory_ (traj.clone ())
  {
    assert (s != 0. && !std::isinf (s) && !std::isnan (s));
  }

  template <unsigned dorder>
  FreeTimeTrajectory<dorder>::FreeTimeTrajectory
  (const FreeTimeTrajectory<dorder>& traj)
    throw ()
    : Trajectory<dorder> (traj.timeRange (), traj.outputSize (),
			  traj.parameters ()),
      trajectory_ (traj.trajectory_->clone ())
  {
  }

  template <unsigned dorder>
  FreeTimeTrajectory<dorder>::~FreeTimeTrajectory () throw ()
  {
    delete trajectory_;
  }

  template <unsigned dorder>
  void
  FreeTimeTrajectory<dorder>::impl_compute
  (typename FreeTimeTrajectory<dorder>::result_t& res , double t) const throw ()
  {
    (*trajectory_) (res, this->scaleTime (t));
  }

  template <unsigned dorder>
  void
  FreeTimeTrajectory<dorder>::impl_derivative
  (typename FreeTimeTrajectory<dorder>::gradient_t& derivative,
   double t,
   typename FreeTimeTrajectory<dorder>::size_type order) const throw ()
  {
    assert (order >= 0);
    double scaled = this->scaleTime (t);
    trajectory_->derivative (derivative, scaled, order);
    derivative *= std::pow (this->timeScale (), 0. + order);
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::jacobian_t
  FreeTimeTrajectory<dorder>::variationConfigWrtParam (double t) const throw ()
  {
    using namespace boost::numeric::ublas;
    value_type scaled = this->scaleTime (t);

    double tMin = this->getLowerBound (this->trajectory_->timeRange ());

    jacobian_t result (this->outputSize (),
		       this->parameters ().size ());
    result.clear ();

    // Compute variation w.r.t time scale (p_0)
    column (result, 0) = (t - tMin) * trajectory_->derivative (scaled, 1);

    // Fill 1..(n-1) lines with original jacobian.
    project (result, range (0, result.size1 ()), range (1, result.size2 ()))
      = trajectory_->variationConfigWrtParam (scaled);

    return result;
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::jacobian_t
  FreeTimeTrajectory<dorder>::variationDerivWrtParam (double t, size_type order)
    const throw ()
  {
    if (order == 0)
      return this->variationConfigWrtParam (t);

    using namespace boost::numeric::ublas;
    value_type scaled = this->scaleTime (t);

    double tMin = this->getLowerBound (this->trajectory_->timeRange ());

    jacobian_t result (this->outputSize (),
		       this->parameters ().size ());
    result.clear ();

    // Compute variation w.r.t time scale (p_0)
    column (result, 0) = trajectory_->derivative (scaled, 1) * order;
    column (result, 0) += trajectory_->derivative (scaled, 1)
      * this->timeScale () * (t - tMin);
    column (result, 0) *= std::pow (this->timeScale (), order - 1.);

    // Fill 1..(n-1) lines with original jacobian.
    project (result, range (0, result.size1 ()), range (1, result.size2 ()))
      = trajectory_->variationDerivWrtParam (scaled, order)
      * std::pow (this->timeScale (), 0. + order);

    return result;
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::value_type
  FreeTimeTrajectory<dorder>::singularPointAtRank (size_type rank) const
  {
    double tMin = this->getLowerBound (this->timeRange ());
    return tMin + (trajectory_->singularPointAtRank (rank) - tMin)
      * this->timeScale ();
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::vector_t
  FreeTimeTrajectory<dorder>::derivBeforeSingularPoint (size_type rank,
							size_type order)
    const
  {
    return trajectory_->derivBeforeSingularPoint (rank, order);
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::vector_t
  FreeTimeTrajectory<dorder>::derivAfterSingularPoint (size_type rank, size_type order)
    const
  {
    return trajectory_->derivAfterSingularPoint (rank, order);
  }

  template <unsigned dorder>
  void
  FreeTimeTrajectory<dorder>::setParameters (const vector_t& p) throw ()
  {
    //FIXME: is this ok?
    vector_t p_ = p;
    if (p_[0] <= 0.)
      p_[0] = 1e-8;

    this->parameters_ = p_;
    this->timeRange_ = detail::scaleInterval (*trajectory_, this->timeScale ());
    this->trajectory_->setParameters (removeScaleFromParameters (p_));
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::value_type
  FreeTimeTrajectory<dorder>::timeScale () const throw ()
  {
    return this->parameters_[0];
  }

  template <unsigned dorder>
  double
  FreeTimeTrajectory<dorder>::scaleTime (double unscaled) const throw ()
  {
    value_type tMin = getLowerBound (this->timeRange ());
    value_type tmin = getLowerBound (this->trajectory_->timeRange ());
    value_type tmax = getUpperBound (this->trajectory_->timeRange ());

    unscaled = detail::fixTime (unscaled, *this);
    assert (this->isValidTime (unscaled));

    value_type res = tmin + (unscaled - tMin) * timeScale ();

    if (res > tmax)
      res = tmax;
    else if (res < tmin)
      res = tmin;
    return res;
  }

  template <unsigned dorder>
  double
  FreeTimeTrajectory<dorder>::unscaleTime (double scaled) const throw ()
  {
    value_type tMin = getLowerBound (this->timeRange ());
    value_type tMax = getUpperBound (this->timeRange ());
    value_type tmin = getLowerBound (this->trajectory_->timeRange ());

    scaled = detail::fixTime (scaled, *this);
    assert (trajectory_->isValidTime (scaled));

    value_type res = tMin + (scaled - tmin) / timeScale ();

    if (res > tMax)
      res = tMax;
    else if (res < tMin)
      res = tMin;
    return res;
  }

  template <unsigned dorder>
  std::ostream&
  FreeTimeTrajectory<dorder>::print (std::ostream& o) const throw ()
  {
    o << "Free time trajectory." << std::endl;
    return o;
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_FREETIMETRAJECTORY_HXX
