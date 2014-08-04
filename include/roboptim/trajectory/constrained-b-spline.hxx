// Copyright (C) 2013 by Alexander Werner, DLR.
// Copyright (C) 2014 by Benjamin Chrétien, CNRS-LIRMM.
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

#ifndef ROBOPTIM_TRAJECTORY_CONSTRAINED_B_SPLINE_HXX
# define ROBOPTIM_TRAJECTORY_CONSTRAINED_B_SPLINE_HXX

# include <boost/shared_ptr.hpp>

# include <Eigen/SVD>

namespace roboptim
{
namespace trajectory
{
  template <int N>
  ConstrainedBSpline<N>::ConstrainedBSpline (interval_t timeRange,
					     size_type dimension,
					     const vector_t& parameters,
					     const std::string name)
    : BSpline<N> (timeRange, dimension, parameters, name),
      constraints_ (0, parameters.rows ()),
      constraint_values_ (),
      tunables_ (parameters),
      projector_ (),
      projector_offset_ ()
  {}

  template <int N>
  ConstrainedBSpline<N>::ConstrainedBSpline (interval_t timeRange,
					     size_type dimension,
					     const vector_t& parameters,
					     const vector_t& knots,
					     const std::string name)
    : BSpline<N> (timeRange, dimension, parameters, knots, name),
      constraints_ (0, parameters.rows ()),
      constraint_values_ (),
      tunables_ (parameters),
      projector_ (),
      projector_offset_ ()
  {}

  template <int N>
  ConstrainedBSpline<N>::~ConstrainedBSpline () {}

  template <int N>
  void ConstrainedBSpline<N>::addFixedConstraint (double t,
                                                  size_type dimension,
                                                  value_type value,
                                                  size_type order)
  {
    // Constraint_values = constraints_ * parameters_ (1)
    const typename vector_t::Index existing_constraints = constraints_.rows ();
    // Insert constraint into new row of (1)
    constraint_values_.conservativeResize (existing_constraints + 1);
    constraint_values_ (existing_constraints) = value;

    constraints_.conservativeResize (existing_constraints + 1,
                                     Eigen::NoChange);
    constraints_.bottomRows (1).setZero ();

    t = detail::fixTime (t, *this);
    const size_type k = this->interval (t);
    const std::size_t k_ = static_cast<std::size_t> (k);
    const size_type n = this->outputSize ();

    for (size_type idx = 0; idx < this->order_ + 1; idx++)
      {
	const std::size_t idx_ = static_cast<std::size_t> (idx);
        constraints_ (existing_constraints, (k - idx) * n + dimension)
	  = this->basisPolynomials ()[k_ - idx_][idx_].derivative (t, order);
      }
    updateProjector ();
  }


  template <int N>
  void ConstrainedBSpline<N>::addCoupledConstraint
  (value_type t_1, size_type dimension_1,
   value_type t_2, size_type dimension_2,
   size_type derivative, value_type factor)
  {
    // Constraint_values = constraints_ * parameters_ (1)
    const typename vector_t::Index existing_constraints = constraints_.rows ();

    // Insert constraint into new row of (1)
    constraint_values_.conservativeResize (existing_constraints + 1);
    constraint_values_ (existing_constraints) = 0.;

    constraints_.conservativeResize (existing_constraints + 1,
                                     Eigen::NoChange);
    constraints_.bottomRows (1).setZero ();

    t_1 = detail::fixTime (t_1, *this);
    t_2 = detail::fixTime (t_2, *this);
    const size_type k_1 = this->interval (t_1);
    const std::size_t k_1_ = static_cast<std::size_t> (k_1);
    const size_type k_2 = this->interval (t_2);
    const std::size_t k_2_ = static_cast<std::size_t> (k_2);
    const size_type n = this->outputSize ();

    for (size_type idx = 0; idx < this->order_ + 1; idx++)
      {
	const std::size_t idx_ = static_cast<std::size_t> (idx);
        constraints_ (existing_constraints, (k_1 - idx) * n + dimension_1)
	  = this->basisPolynomials ()
	  [k_1_ - idx_][idx_].derivative (t_1, derivative);
        constraints_ (existing_constraints, (k_2 - idx) * n + dimension_2)
	  -= factor *
	  this->basisPolynomials ()[k_2_ - idx_][idx_].derivative (t_2, derivative);
      }
    updateProjector ();
  }


  template <int N>
  const typename ConstrainedBSpline<N>::vector_t&
  ConstrainedBSpline<N>::parameters () const
  {
    return tunables_;
  }

  template <int N>
  void ConstrainedBSpline<N>::setParameters (const vector_t& parameters_in)
  {
    assert (tunables_.rows () == parameters_in.rows ());
    tunables_ = parameters_in;

    this->parameters_ = projector_offset_;
    this->parameters_.noalias () += projector_ * tunables_;
  }

  template <int N>
  void ConstrainedBSpline<N>::updateProjector ()
  {
    Eigen::JacobiSVD<matrix_t> svd (constraints_.rows (),
                                    constraints_.cols ());
    svd.compute (constraints_, Eigen::ComputeFullU | Eigen::ComputeFullV);

    // Get the dimension of the null space
    const typename vector_t::Index null_space_dim = svd.matrixV ().rows ()
      - svd.singularValues ().rows ();

    // Otherwise no parameters to optimize are left
    assert (null_space_dim > 0);

    projector_ = svd.matrixV ().rightCols (null_space_dim);
    projector_offset_ = svd.solve (constraint_values_);

    // Preserve parameters as good as possible (not violating the constraints)
    // Projector should be invertible at any time and it should also be a base
    // => transpose
    // This also updated the size of tunables_ to the correct value
    tunables_ = projector_.transpose ()
      * (this->parameters_ - projector_offset_);
    this->parameters_ = projector_offset_;
    this->parameters_.noalias () += projector_ * tunables_;
  }

  template <int N>
  Trajectory<N>*  ConstrainedBSpline<N>::resize (interval_t timeRange)
    const
  {
    // Create an unconstrained B-spline
    ConstrainedBSpline<N>* c
      = new ConstrainedBSpline<N> (timeRange,
                                   this->outputSize (),
                                   this->parameters_,
                                   this->knotVector (),
                                   this->getName ());

    // Copy constraints from this instance
    c->constraints_ = this->constraints_;
    c->constraint_values_ = this->constraint_values_;
    c->tunables_ = this->tunables_;

    return c;
  }

  template <int N>
  typename ConstrainedBSpline<N>::jacobian_t
  ConstrainedBSpline<N>::variationDerivWrtParam
  (double t, size_type order) const
  {
    const size_type n = this->outputSize ();
    jacobian_t jac_basispolynomials (n, this->getNumberControlPoints () * n);
    jac_basispolynomials.setZero ();
    const size_type k = this->interval (t);
    const std::size_t k_ = static_cast<std::size_t> (k);
    for (size_type idx = 0; idx < this->order_ + 1; idx++)
      {
	const std::size_t idx_ = static_cast<std::size_t> (idx);
        const Polynomial<N>& B = this->basisPolynomials ()[k_ - idx_][idx_];
        jac_basispolynomials.middleCols ((k - idx) * n, n).diagonal ().
	  setConstant (B.derivative (t, order));
      }
    const typename vector_t::Index variables = tunables_.rows ();
    jacobian_t jac (n, variables);
    jac = jac_basispolynomials * this->projector_;
    return jac;
  }

} // end of namespace trajectory.
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_CONSTRAINED_B_SPLINE_HXX
