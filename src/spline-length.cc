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

#include <roboptim/trajectory/spline-length.hh>

namespace roboptim
{
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

    using namespace boost;
    using namespace boost::numeric::ublas;

    for (value_type t = get<0> (interval_); t <= get<1> (interval_);
	 t += get<2> (interval_))
      {
	double tmp = norm_2 (traj.derivative (t, 2));
	res[0] += tmp * tmp;
      }
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

    using namespace boost;
    using namespace boost::numeric::ublas;

    for (value_type t = get<0> (interval_); t <= get<1> (interval_);
	 t += get<2> (interval_))
      noalias (grad) += prod (traj.derivative (t, 2),
			      traj.variationDerivWrtParam (t, 2));
  }
} // end of namespace roboptim.
