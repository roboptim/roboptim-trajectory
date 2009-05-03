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
  template <unsigned dorder>
  FreeTimeTrajectory<dorder>::FreeTimeTrajectory
  (const Trajectory<dorder>& traj, double s)
    throw ()
    : Trajectory<dorder>
      (std::make_pair (traj.timeRange ().first,
		       scaleTime (traj.timeRange ().second)),
       traj.m, traj.parameters ()),
      trajectory_ (traj),
      timeScale_ (s)
  {
    assert (timeScale_ != 0.
	    && std::isinf (timeScale_)
	    && std::isnan (timeScale_));
  }


  template <unsigned dorder>
  FreeTimeTrajectory<dorder>::~FreeTimeTrajectory () throw ()
  {
  }


  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::vector_t
  FreeTimeTrajectory<dorder>::operator () (double t) const throw ()
  {
    return trajectory_ (scaleTime (t));
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::vector_t
  FreeTimeTrajectory<dorder>::derivative (double t, size_type order) const throw ()
  {
    double scaled = scaleTime (t);
    vector_t d = trajectory_.deriv (scaled, order);

    double coef = 1.;
    for (size_t i = 0; i < order; ++i)
      coef *= ;
    outDerivative *= coef;
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::jacobian_t
  FreeTimeTrajectory<dorder>::variationConfigWrtParam (double t) const throw ()
  {
    double t_min = timeRange ().first;
    double scaled = scaleTime (t);
    jacobian_t jac = trajectory_.variationConfigWrtParam (scale);

    // Last column corresponds to derivative wrt lambda_{p+1}
    vector_t dGamma0_dt = trajectory_.derivative (scaled, 1);
    dGamma0_dt *= (t - t_min);

    // Fill last column of Jacobian with dGamma0_dt
    size_type timeScalingIndex = getTimeScalingIndex ();
    for (size_type i = 0; i < jac.size1(); ++i)
      jac (i, timeScalingIndex) = dGamma0_dt (i);
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::jacobian_t
  FreeTimeTrajectory<dorder>::variationDerivWrtParam (double t, size_type order)
    const throw ()
  {
    if (!order)
      return variationConfigWrtParam (t);

    double t_min = timeRange ().first;

    double scaled = scaleTime (t);
    jacobian_t jac =
      trajectory_.variationDerivWrtParam (scaled, order);

    double lambda_p_plus_1_pow_n_minus_1 = 1.;
    for (size_type i = 0; i< order - 1; ++i)
      lambda_p_plus_1_pow_n_minus_1 *= timeScale_;
  }
  outJacobian *= lambda_p_plus_1_pow_n_minus_1 * timeScale_;

  // Fill last column of Jacobian
  vector_t lastColumn1 = trajectory_.derivative (scaled, order, lastColumn1);
  lastColumn1 *= order;

  vector_t lastColumn2 = trajectory_.derivative (scaled, order, lastColumn2);
  lastColumn2 *= timeScale_ * (t-t_min);

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
    double t_min = timeRange ().first;
    return t_min + (trajectory_.singularPointAtRank (rank) - t_min)
      / timeScale_;
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::vector_t
  FreeTimeTrajectory<dorder>::derivBeforeSingularPoint (size_type rank,
							size_type order)
    const
  {
    return trajectory_.derivBeforeSingularPoint (rank, order);
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::vector_t
  FreeTimeTrajectory<dorder>::derivAfterSingularPoint (size_type rank, size_type order)
    const
  {
    return trajectory_.derivAfterSingularPoint (rank, order);
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::double
  FreeTimeTrajectory<dorder>::timeScale () const throw ()
  {
    return timeScale_;
  }

  template <unsigned dorder>
  typename FreeTimeTrajectory<dorder>::double
  FreeTimeTrajectory<dorder>::scaleTime (double unscaled) const throw ()
  {
    double min = timeRange ().first;
    double scaled = min + timeScale_ * (unscaled - min);
    return scaled;
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
