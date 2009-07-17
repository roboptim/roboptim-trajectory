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

#include <boost/numeric/ublas/vector.hpp>
#include <boost/scoped_ptr.hpp>

#include <roboptim/core/finite-difference-gradient.hh>

#include <roboptim/trajectory/trajectory.hh>

#include "roboptim/trajectory/limit-speed.hh"

namespace roboptim
{
  LimitSpeed::LimitSpeed (StableTimePoint timePoint,
			  const GenericTrajectory& spline) throw ()
    : DerivableFunction (1 + spline.parameters ().size (), 1, "limit speed"),
      timePoint_ (timePoint),
      trajectory_ (spline)
  {}

  LimitSpeed::~LimitSpeed () throw ()
  {}

  void
  LimitSpeed::impl_compute (result_t& res, const argument_t& p) const throw ()
  {
    using namespace boost::numeric::ublas;
    res.clear ();

    boost::scoped_ptr<GenericTrajectory> updatedSpline (trajectory_.clone ());
    updatedSpline->setParameters (p);

    //FIXME: should be done in updatedSpline!
    value_type t =
      timePoint_.getTime (getUpperBound (updatedSpline->timeRange ()));
    res[0] = norm_2 (updatedSpline->derivative (t));
  }

  void
  LimitSpeed::impl_gradient (gradient_t& grad, const argument_t& p, size_type i)
    const throw ()
  {
    assert (i == 0);
    grad.clear ();

    //FIXME: compute gradient analytically.
    FiniteDifferenceGradient fdfunction (*this);
    fdfunction.gradient (grad, p, 0);
  }
} // end of namespace roboptim.
