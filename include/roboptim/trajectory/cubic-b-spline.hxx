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

namespace roboptim
{
  template <typename P>
  void
  CubicBSpline::frezzeCurveStart (P& problem)
  {
    using boost::numeric::ublas::zero_matrix;
    using boost::numeric::ublas::zero_vector;

    Function::matrix_t A[3];
    Function::vector_t b[3];

    A[0] = zero_matrix<double>(1, paramSize);
    A[1] = zero_matrix<double>(1, paramSize);
    A[2] = zero_matrix<double>(1, paramSize);

    b[0] = zero_vector<double>(1);
    b[1] = zero_vector<double>(1);
    b[2] = zero_vector<double>(1);
    A[0](0,1) = A[1](0,2) = A[2](0,3)=1./6.;
    A[0](0,4) = A[1](0,5) = A[2](0,6)=2./3.;
    A[0](0,7) = A[1](0,8) = A[2](0,9)=1./6.;
    b[0](0) = -freeTimeTraj_.parameters ()[1+0];
    b[1](0) = -freeTimeTraj_.parameters ()[1+1];
    b[2](0) = -freeTimeTraj_.parameters ()[1+2];
    NumericLinearFunction* boundaryCond;
    boundaryCond = new NumericLinearFunction(A[0], b[0]);
    boost::shared_ptr<LinearFunction>	boundaryCond0ShPtr(boundaryCond);
    boundaryCond = new NumericLinearFunction(A[1], b[1]);
    boost::shared_ptr<LinearFunction>	boundaryCond1ShPtr(boundaryCond);
    boundaryCond = new NumericLinearFunction(A[2], b[2]);
    boost::shared_ptr<LinearFunction>	boundaryCond2ShPtr(boundaryCond);

    problem.addConstraint(boundaryCond0ShPtr, interval);
    problem.addConstraint(boundaryCond1ShPtr, interval);
    problem.addConstraint(boundaryCond2ShPtr, interval);
  }

  template <typename P>
  void
  CubicBSpline::frezzeCurveEnd (P& problem)
  {
    using boost::numeric::ublas::zero_matrix;
    using boost::numeric::ublas::zero_vector;

    Function::matrix_t A[3];
    Function::vector_t b[3];

    A[0] = zero_matrix<double>(1, paramSize);
    A[1] = zero_matrix<double>(1, paramSize);
    A[2] = zero_matrix<double>(1, paramSize);
    b[0] = zero_vector<double>(1);
    b[1] = zero_vector<double>(1);
    b[2] = zero_vector<double>(1);

    A[0](0,paramSize-9) = A[1](0,paramSize-8) = A[2](0,paramSize-7)=1./6.;
    A[0](0,paramSize-6) = A[1](0,paramSize-5) = A[2](0,paramSize-4)=2./3.;
    A[0](0,paramSize-3) = A[1](0,paramSize-2) = A[2](0,paramSize-1)=1./6.;
    b[0](0) = -freeTimeTraj_.parameters ()[paramSize-3];
    b[1](0) = -freeTimeTraj_.parameters ()[paramSize-2];
    b[2](0) = -freeTimeTraj_.parameters ()[paramSize-1];
    boundaryCond = new NumericLinearFunction(A[0], b[0]);
    boost::shared_ptr<LinearFunction>	boundaryCond3ShPtr(boundaryCond);
    boundaryCond = new NumericLinearFunction(A[1], b[1]);
    boost::shared_ptr<LinearFunction>	boundaryCond4ShPtr(boundaryCond);
    boundaryCond = new NumericLinearFunction(A[2], b[2]);
    boost::shared_ptr<LinearFunction>	boundaryCond5ShPtr(boundaryCond);

    Function::interval_t interval = Function::makeInterval(0., 0.);
    problem.addConstraint(boundaryCond3ShPtr, interval);
    problem.addConstraint(boundaryCond4ShPtr, interval);
    problem.addConstraint(boundaryCond5ShPtr, interval);
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HXX
