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

#include <boost/numeric/ublas/vector_expression.hpp>
#include <roboptim/trajectory/spline-length.hh>

namespace roboptim
{
  namespace
  {
    struct SumLength
    {
      SumLength (const Spline& traj, double& res)
	: traj_ (traj),
	  res_ (res)
      {}

      void operator () (const double& t)
      {
	using namespace boost::numeric::ublas;
	res_ += prec_inner_prod (traj_.derivative (t, 1),
				 traj_.derivative (t, 1));
      }

    private:
      const Spline& traj_;
      double& res_;
    };

    struct SumLengthGrad
    {
      SumLengthGrad (const Spline& traj,
		     SplineLength::gradient_t& grad)
	: traj_ (traj),
	  grad_ (grad)
      {}

      void operator () (const double& t)
      {
	using namespace boost::numeric::ublas;
	noalias (grad_) += prec_prod (traj_.derivative (t, 1),
				      traj_.variationDerivWrtParam (t, 1));
      }

    private:
      const Spline& traj_;
      SplineLength::gradient_t& grad_;
    };

  }

  SplineLength::SplineLength (const Spline& spline,
			      discreteInterval_t interval)
    throw ()
    : TrajectoryCost<Spline> (spline, "spline length"),
      interval_ (interval)
  {
  }

  SplineLength::~SplineLength () throw ()
  {
  }

  void
  SplineLength::impl_compute (result_t& res, const argument_t& p)
    const throw ()
  {
    trajectory_t traj = trajectory_;
    traj.setParameters (p);

    SumLength sumlength (traj, res[0]);
    foreach (interval_, sumlength);
    res[0] *= getUpperBound (interval_) - getLowerBound (interval_);
    res[0] /= 2.;
  }

  void
  SplineLength::impl_gradient (gradient_t& grad, const argument_t& p,
			       size_type i)
    const throw ()
  {
    assert (i == 0);
    grad.clear ();

    trajectory_t traj = trajectory_;
    traj.setParameters (p);

    SumLengthGrad sumlengthgrad (traj, grad);
    foreach (interval_, sumlengthgrad);
    grad *= getUpperBound (interval_) - getLowerBound (interval_);
  }
} // end of namespace roboptim.
