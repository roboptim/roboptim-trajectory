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

namespace roboptim
{
  namespace detail
  {
    Function::vector_t
    addScaleToParams (const Function::vector_t& p,
		      Function::value_type t = 1.);

    Function::vector_t
    removeScaleFromParams (const Function::vector_t& v);


    Function::vector_t
    addScaleToParams (const Function::vector_t& p,
		     Function::value_type t)
    {
      Function::vector_t res (p.size () + 1);
      res[0] = t;
      for (unsigned i = 0; i < p.size (); ++i)
	res[i + 1] = p[i];
      return res;
    }

    Function::vector_t
    removeScaleFromParams (const Function::vector_t& p)
    {
      Function::vector_t res (p.size () - 1);
      for (unsigned i = 1; i < p.size (); ++i)
	res[i - 1] = p[i];
      return res;
    }

    Function::value_type
    scaleTime (Function::value_type unscaled,
	       Function::value_type min,
	       Function::value_type scale)
    {
      return min + scale * (unscaled - min);
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
  (const Trajectory<dorder>& traj, double s)
    throw ()
    : Trajectory<dorder> (detail::scaleInterval (traj, s), traj.outputSize (),
			  detail::addScaleToParams (traj.parameters (), s)),
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
    (*trajectory_) (scaleTime (t));
  }

  template <unsigned dorder>
  void
  FreeTimeTrajectory<dorder>::impl_derivative
  (typename FreeTimeTrajectory<dorder>::gradient_t& derivative,
   double t,
   typename FreeTimeTrajectory<dorder>::size_type order) const throw ()
  {
    double scaled = scaleTime (t);
    trajectory_->derivative (derivative, scaled, order);

    double coef = 1.;
    for (size_t i = 0; i < order; ++i)
      coef *= this->timeScale (); //FIXME: timeScale?
    derivative *= coef;
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::jacobian_t
  FreeTimeTrajectory<dorder>::variationConfigWrtParam (double t) const throw ()
  {
    double t_min = this->getLowerBound (this->timeRange ());
    double scaled = scaleTime (t);
    jacobian_t jac = trajectory_->variationConfigWrtParam (scaled);

    // Last column corresponds to derivative wrt lambda_{p+1}
    vector_t dGamma0_dt = trajectory_->derivative (scaled, 1);
    dGamma0_dt *= (t - t_min);

    // Fill last column of Jacobian with dGamma0_dt
    size_type timeScalingIndex = getTimeScalingIndex ();
    for (size_type i = 0; i < jac.size1(); ++i)
      jac (i, timeScalingIndex) = dGamma0_dt (i);
    return jac;
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::jacobian_t
  FreeTimeTrajectory<dorder>::variationDerivWrtParam (double t, size_type order)
    const throw ()
  {
    if (!order)
      return variationConfigWrtParam (t);

    double t_min = this->getLowerBound (this->timeRange ());

    double scaled = scaleTime (t);
    jacobian_t jac =
      trajectory_->variationDerivWrtParam (scaled, order);

    double lambda_p_plus_1_pow_n_minus_1 = 1.;
    for (size_type i = 0; i< order - 1; ++i)
      lambda_p_plus_1_pow_n_minus_1 *= this->timeScale ();

    jac *= lambda_p_plus_1_pow_n_minus_1 * this->timeScale ();

    // Fill last column of Jacobian
    vector_t lastColumn1 = trajectory_->derivative (scaled, order);
    lastColumn1 *= order;

    vector_t lastColumn2 = trajectory_->derivative (scaled, order);
    lastColumn2 *= this->timeScale () * (t-t_min);

    vector_t lastColumn =
      (lastColumn1 + lastColumn2) * lambda_p_plus_1_pow_n_minus_1;

    size_type timeScalingIndex = getTimeScalingIndex ();
    for (size_type i = 0; i < jac.size1(); ++i)
      jac (i, timeScalingIndex) = lastColumn (i);
    return jac;
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::value_type
  FreeTimeTrajectory<dorder>::singularPointAtRank (size_type rank) const
  {
    double t_min = this->getLowerBound (this->timeRange ());
    return t_min + (trajectory_->singularPointAtRank (rank) - t_min)
      / this->timeScale ();
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
    this->parameters_ = p;
    this->trajectory_->setParameters (detail::removeScaleFromParams (p));
  }

  template <unsigned dorder>
  double
  FreeTimeTrajectory<dorder>::timeScale () const throw ()
  {
    return this->parameters_[0];
  }

  template <unsigned dorder>
  double
  FreeTimeTrajectory<dorder>::scaleTime (double unscaled) const throw ()
  {
    return detail::scaleTime (unscaled,
			      getLowerBound (this->timeRange ()),
			      this->timeScale ());
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
