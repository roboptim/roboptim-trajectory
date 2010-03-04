// Copyright (C) 2010 by Thomas Moulard, AIST, CNRS, INRIA.
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

#ifndef ROBOPTIM_TRAJECTORY_CUBIC_B_SPLINE_HXX
# define ROBOPTIM_TRAJECTORY_CUBIC_B_SPLINE_HXX
# include <roboptim/core/numeric-linear-function.hh>

namespace roboptim
{
  template <typename P>
  void
  CubicBSpline::freezeCurveStart (P& problem, size_type offset) const throw ()
  {
    using boost::shared_ptr;
    const size_type paramSize = parameters ().size ();
    const Function::interval_t interval = Function::makeInterval (0, 0);
    assert (paramSize >= offset + 3 * outputSize ());

    for (size_type i = 0; i < outputSize (); ++i)
      {
	Function::matrix_t A (1, offset + paramSize);
	Function::vector_t b (1);
	A.clear ();
	b.clear ();

	A (0, offset + i + 0. * outputSize ()) = 1. / 6.;
	A (0, offset + i + 1. * outputSize ()) = 2. / 3.;
	A (0, offset + i + 2. * outputSize ()) = 1. / 6.;
	b (0) = -parameters ()[offset + i];
	NumericLinearFunction* boundaryCond = new NumericLinearFunction (A, b);
	shared_ptr<LinearFunction> boundaryCondShPtr (boundaryCond);
	problem.addConstraint (boundaryCondShPtr, interval);
      }
  }

  template <typename P>
  void
  CubicBSpline::freezeCurveEnd (P& problem, size_type offset) const throw ()
  {
    using boost::shared_ptr;
    const size_type paramSize = parameters ().size ();
    const Function::interval_t interval = Function::makeInterval (0, 0);
    assert (paramSize >= offset + 3 * outputSize ());

    for (size_type i = 0; i < outputSize (); ++i)
      {
	Function::matrix_t A (1, offset + paramSize);
	Function::vector_t b (1);
	A.clear ();
	b.clear ();

	A (0, offset + paramSize - 1 - i - 0. * outputSize ()) = 1. / 6.;
	A (0, offset + paramSize - 1 - i - 1. * outputSize ()) = 2. / 3.;
	A (0, offset + paramSize - 1 - i - 2. * outputSize ()) = 1. / 6.;
	b (0) = -parameters ()[paramSize - 1 - i];
	NumericLinearFunction* boundaryCond = new NumericLinearFunction (A, b);
	shared_ptr<LinearFunction> boundaryCondShPtr (boundaryCond);
	problem.addConstraint (boundaryCondShPtr, interval);
      }
  }
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HXX
