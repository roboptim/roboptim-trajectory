// Copyright (C) 2013 by Alexander Werner, DLR.
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

#ifndef ROBOPTIM_TRAJECTORY_CONSTRAINT_B_SPLINE_HXX
# define ROBOPTIM_TRAJECTORY_CONSTRAINT_B_SPLINE_HXX
# include <boost/shared_ptr.hpp>
# include <boost/numeric/conversion/converter.hpp>
# include <roboptim/core/numeric-linear-function.hh>

# include <Eigen/SVD>

namespace roboptim
{
  template<int N>
  ConstraintBSpline<N>::ConstraintBSpline (interval_t timeRange,
					   size_type dimension,
					   const vector_t& parameters,
					   const std::string name) throw()
    : BSpline<N> (timeRange, dimension, parameters, name),
      constraints_ (0, parameters.rows()), tunables_ (parameters)
  {}

  template<int N>
  ConstraintBSpline<N>::ConstraintBSpline (interval_t timeRange,
					   size_type dimension,
					   const vector_t& parameters,
					   const vector_t& knots,
					   const std::string name)throw()
    : BSpline<N> (timeRange, dimension, parameters, knots, name),
      constraints_ (0, parameters.rows()), tunables_ (parameters)
  {}

  template<int N>
  ConstraintBSpline<N>::~ConstraintBSpline() throw() {}

  template<int N>
  void ConstraintBSpline<N>::addfixedConstraint (double t, size_type dimension,
						 value_type value,
						 size_type order) throw()
  {
    // constraint_values = constraints_ * parameters_ (1)
    const int existing_constraints = constraints_.rows();
    //insert constraint into new row of (1)
    constraint_values_.conservativeResize (existing_constraints + 1);
    constraint_values_ (existing_constraints) = value;

    constraints_.conservativeResize (existing_constraints + 1, Eigen::NoChange);
    constraints_.bottomRows (1).setZero();

    t = detail::fixTime (t, *this);
    const size_type k = this->interval (t);
    const size_type n = this->outputSize ();

    for (size_type idx = 0; idx < this->order_ + 1; idx++)
      {
        constraints_ (existing_constraints, (k - idx) * n + dimension)
	  = this->basisPolynomials_[k - idx][idx].derivative (t, order);
      }
    updateProjector();
  }


  template<int N>
  void ConstraintBSpline<N>::addCoupledConstraint (
						   value_type t_1, size_type dimension_1,
						   value_type t_2, size_type dimension_2,
						   size_type derivative, value_type factor) throw()
  {
    // constraint_values = constraints_ * parameters_ (1)
    const int existing_constraints = constraints_.rows();
    //insert constraint into new row of (1)
    constraint_values_.conservativeResize (existing_constraints + 1);
    constraint_values_ (existing_constraints) = 0.;

    constraints_.conservativeResize (existing_constraints + 1, Eigen::NoChange);
    constraints_.bottomRows (1).setZero();

    t_1 = detail::fixTime (t_1, *this);
    t_2 = detail::fixTime (t_2, *this);
    const size_type k_1 = this->interval (t_1);
    const size_type k_2 = this->interval (t_2);
    const size_type n = this->outputSize ();

    for (size_type idx = 0; idx < this->order_ + 1; idx++)
      {
        constraints_ (existing_constraints, (k_1 - idx) * n + dimension_1)
	  = this->basisPolynomials_[k_1 - idx][idx].derivative (t_1, derivative);
        constraints_ (existing_constraints, (k_2 - idx) * n + dimension_2)
	  += factor *
	  this->basisPolynomials_[k_2 - idx][idx].derivative (t_2, derivative);
      }
    updateProjector();
  }


  template<int N>
  typename ConstraintBSpline<N>::vector_t const&
  ConstraintBSpline<N>::parameters() const throw()
  {
    return tunables_;
  }

  template<int N>
  void ConstraintBSpline<N>::setParameters (const vector_t& parameters_in)
    throw ()
  {
    assert (tunables_.rows() == parameters_in.rows());
    tunables_ = parameters_in;

    this->parameters_ = projector_offset_;
    this->parameters_.noalias() += projector_ * tunables_;
  }

  template<int N>
  void ConstraintBSpline<N>::updateProjector()
  {
    Eigen::JacobiSVD<matrix_t> svd (constraints_.rows(), constraints_.cols());
    svd.compute (constraints_, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const int null_space_dim = svd.matrixV().rows() - svd.singularValues().rows();
    assert (null_space_dim > 0);	//otherwise no parameters to optimize are left
    projector_ = svd.matrixV().rightCols (null_space_dim);
    projector_offset_ = svd.solve (constraint_values_);
    //preserve parameters as good as possible ( not violating the constraints)
    //projector should be invertible at any time and it should also be a base => transpose
    //this also updated the size of tunables_ to the correct value
    tunables_ = projector_.transpose() * (this->parameters_ - projector_offset_);
    this->parameters_ = projector_offset_;
    this->parameters_.noalias() += projector_ * tunables_;
  }

  template<int N>
  Trajectory<N>*  ConstraintBSpline<N>::resize
  (interval_t timeRange) const throw ()
  {
    //create a non-constraint spline
    ConstraintBSpline<N>* c = new ConstraintBSpline<N> (
							timeRange,
							this->outputSize(),
							this->parameters_,
							this->knots(),
							this->getName()
							);
    //copy constraints from this instance
    c->constraints_ = this->constraints_;
    c->constraint_values_ = this->constraint_values_;
    c->tunables_ = this->tunables_;
    return c;
  }

  template<int N>
  typename ConstraintBSpline<N>::jacobian_t
  ConstraintBSpline<N>::variationDerivWrtParam
  (double t, size_type order) const throw ()
  {
    const size_type n = this->outputSize ();
    jacobian_t jac_basispolynomials (n, this->nbp_ * n);
    jac_basispolynomials.setZero();
    const size_type k = this->interval (t);
    for (size_type idx = 0; idx < this->order_ + 1; idx++)
      {
        const Polynomial<N>& B = this->basisPolynomials_[k - idx][idx];
        jac_basispolynomials.middleCols ((k - idx)*n, n).diagonal().
	  setConstant (
		       B.derivative (t, order)
		       );
      }
    const int variables = tunables_.rows();
    jacobian_t jac (n, variables);
    jac = jac_basispolynomials * this->projector_;
    return jac;
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_CONSTRAINT_B_SPLINE_HXX
