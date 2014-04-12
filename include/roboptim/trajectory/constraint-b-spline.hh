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

#ifndef ROBOPTIM_TRAJECTORY_CONSTRAINT_B_SPLINE_HH
# define ROBOPTIM_TRAJECTORY_CONSTRAINT_B_SPLINE_HH
# include <roboptim/trajectory/sys.hh>
# include <roboptim/trajectory/deprecated.hh>

# include <roboptim/trajectory/b-spline.hh>

# include <Eigen/Sparse>

namespace roboptim
{
  /// \addtogroup roboptim_function
  /// @{

  /// B-Spline trajectory.

  /// Implement a B-Spline as a trajectory as described in
  /// doc/quintic-b-spline.tex
  template<int N>
  class ConstraintBSpline : public BSpline<N>
  {
  public:
    typedef BSpline<N> parent_t;
    typedef typename parent_t::interval_t interval_t;
    typedef typename parent_t::size_type size_type;
    typedef typename parent_t::value_type value_type;
    typedef typename parent_t::vector_t vector_t;
    typedef typename parent_t::matrix_t matrix_t;
    typedef typename parent_t::jacobian_t jacobian_t;

    /// \brief see B-Spline constructors for documentation
    ConstraintBSpline (interval_t timeRange, size_type dimension,
		       const vector_t& parameters,
		       const std::string name = "Constrain B-Spline") throw ();

    /// \brief see B-Spline constructors for documentation
    ConstraintBSpline (interval_t tr, size_type dimension,
		       const vector_t& parameters,
		       const vector_t& knots,
		       std::string name = "Constraint B-Spline") throw ();

    virtual ~ConstraintBSpline () throw ();

    /**
     * Creates a constraint on the basic spline.
     * This reduces the number of parameter by one.
     * The constraint equations is:
     * value = \frac{\partial^derivative}{\partial t^derivative} f_{dimension}(t)
     * \param t time in the spline to constrain.
     * \param value Desired spline value at t.
     * \param dimension A zero based info which dimension of the spline to constrain.
     * \param derivative Which derivative of the split to constrain.
     */
    void addfixedConstraint (double t, size_type dimension,
			     value_type value,
			     size_type derivative=0) throw ();

    /**
     * Creates a constraint against another part of the spline.
     * 0 = \frac{\partial^derivative}{\partial t^derivative} f_{dimension_1}(t_1) - factor * \frac{\partial^derivative}{\partial t^derivative} f_{dimension_2}(t_2)
     */
    void addCoupledConstraint
    (value_type t_1, size_type dimension_1,
     value_type t_2, size_type dimension_2,
     size_type derivative=0, value_type factor=1.) throw ();

    /**
     * overloaded parameters method from Trajectory<N>.
     * Returns only the tunable parameters.
     */
    const vector_t& parameters () const throw ();

    /**
     * overloaded setParameters method from Trajectory<N>.
     * Given the free parameters it calculates the spline
     * parameters.
     */
    void setParameters (const vector_t&) throw ();

    virtual Trajectory<N>* resize (interval_t timeRange) const throw ();
    jacobian_t variationDerivWrtParam (double t, size_type order) const throw ();
  protected:
    /**
     * Updates the projector matrix. Called after adding
     * a constraint.
     */
    void updateProjector();

    vector_t constraint_values_;
    matrix_t constraints_;
    vector_t tunables_;
    matrix_t projector_;
    vector_t projector_offset_;
  };

  /// @}

} // end of namespace roboptim.

# include <roboptim/trajectory/constraint-b-spline.hxx>
#endif //! ROBOPTIM_TRAJECTORY_CONSTRAINT_B_SPLINE_HH
