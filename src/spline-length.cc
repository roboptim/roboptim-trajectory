// Copyright (C) 2009 by Thomas Moulard, AIST, CNRS, INRIA.
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

#include "debug.hh"
#include <roboptim/trajectory/sys.hh>

#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/optional.hpp>

#include <roboptim/trajectory/spline-length.hh>

namespace roboptim
{
  namespace
  {
    struct SumLength
    {
      SumLength (const CubicBSpline& traj, double& res)
	: traj_ (traj),
	  res_ (res)
      {}

      void operator () (const double& t)
      {
	res_ += traj_.derivative (t, 1).dot(traj_.derivative (t, 1));
      }

    private:
      const CubicBSpline& traj_;
      double& res_;
    };

    struct SumLengthGrad
    {
      SumLengthGrad (const CubicBSpline& traj,
		     SplineLength::gradient_t& grad)
	: traj_ (traj),
	  grad_ (grad)
      {}

      void operator () (const double& t)
      {
	grad_ += traj_.derivative (t, 1).adjoint ()
	  * traj_.variationDerivWrtParam (t, 1);
      }

    private:
      const CubicBSpline& traj_;
      SplineLength::gradient_t& grad_;
    };

  }

  SplineLength::SplineLength (const CubicBSpline& spline,
			      size_type nDiscretizationPoints,
			      boost::optional<interval_t> interval)
    : TrajectoryCost<CubicBSpline> (spline, "spline length"),
      interval_ (interval ? *interval : spline.timeRange ()),
      nDiscretizationPoints_ (nDiscretizationPoints)
  {
  }

  SplineLength::~SplineLength () throw()
  {
  }

  void
  SplineLength::impl_compute (result_t& res, const argument_t& p)
    const throw()
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
      Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    trajectory_t traj = trajectory_;
    traj.setParameters (p);

    SumLength sumlength (traj, res[0]);
    foreach (interval_, nDiscretizationPoints_, sumlength);

    const value_type delta =
      getUpperBound (interval_) - getLowerBound (interval_);
    res[0] *= delta / (value_type)nDiscretizationPoints_;
    res[0] /= 2.;
  }

  void
  SplineLength::impl_gradient (gradient_t& grad, const argument_t& p,
			       size_type ONLY_DEBUG (i))
    const throw()
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
      Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    assert (i == 0);
    grad.setZero ();

    trajectory_t traj = trajectory_;
    traj.setParameters (p);

    SumLengthGrad sumlengthgrad (traj, grad);
    foreach (interval_, nDiscretizationPoints_, sumlengthgrad);
    const value_type delta =
      getUpperBound (interval_) - getLowerBound (interval_);
    grad *= delta / (value_type)nDiscretizationPoints_;
  }
} // end of namespace roboptim.
