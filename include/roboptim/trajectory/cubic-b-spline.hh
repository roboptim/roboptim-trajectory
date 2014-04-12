// Copyright (C) 2009, 2010 by Thomas Moulard, AIST, CNRS, INRIA.
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

#ifndef ROBOPTIM_TRAJECTORY_CUBIC_B_SPLINE_HH
# define ROBOPTIM_TRAJECTORY_CUBIC_B_SPLINE_HH
# include <roboptim/trajectory/sys.hh>
# include <roboptim/trajectory/deprecated.hh>

# include <roboptim/trajectory/trajectory.hh>
# include <roboptim/trajectory/polynomial-3.hh>
# include <roboptim/trajectory/fwd.hh>


namespace roboptim
{
  /// \addtogroup roboptim_function
  /// @{

  /// Cubic B-Spline trajectory.

  /// Implement a B-Spline as a trajectory as described in doc/cubic-b-spline.tex
  class CubicBSpline : public Trajectory<3>
  {
  public:
    /// \brief Instantiate a cubic B-Spline from its definition.
    ///
    /// \param timeRange spline time range: $\f$[t_3,t_n]\f$
    /// \param dimension spline dimension: \f$n\f$
    /// \param parameters vector of parameters defining control points
    /// \param name function title
    ///
    /// Number of control points is inferred from dimension of dimenion of
    /// parameter vector.
    CubicBSpline (interval_t timeRange, size_type dimension,
		  const vector_t& parameters,
		  const std::string name = "cubic B-Spline") throw ();

    /// \brief Copy constructor.
    /// \param spline spline that will be copied
    CubicBSpline (const CubicBSpline& spline) throw ();

    virtual ~CubicBSpline () throw ();

    /// \brief Modify spline parameters.
    virtual void setParameters (const vector_t&) throw ();

    virtual jacobian_t variationConfigWrtParam (double t) const throw ();
    virtual jacobian_t variationDerivWrtParam (double t, size_type order)
      const throw ();
    virtual value_type singularPointAtRank (size_type rank) const;
    virtual vector_t derivBeforeSingularPoint (size_type rank, size_type order)
      const;
    virtual vector_t derivAfterSingularPoint (size_type rank, size_type order)
      const;

    ROBOPTIM_IMPLEMENT_CLONE(CubicBSpline)

    virtual Trajectory<derivabilityOrder>* resize (interval_t timeRange)
    const throw ()
    {
      return new CubicBSpline (timeRange, this->outputSize (), this->parameters ());
    }

    /// \brief Display the function on the specified output stream.
    ///
    /// \param o output stream used for display
    /// \return output stream
    virtual std::ostream& print (std::ostream& o) const throw ();

    jacobian_t variationConfigWrtParam (StableTimePoint tp) const throw ();
    jacobian_t variationDerivWrtParam (StableTimePoint tp, size_type order)
      const throw ();

    /// \brief Add a constraint to a problem in order to freeze the B-spline
    /// at its start.
    /// \param problem problem to which the constraint will be added.
    /// \param offset offset of the B-spline parameters in the problem's
    /// parameter list.
    template <typename P>
    void freezeCurveStart (P& problem, size_type offset = 0) const throw ();

    /// \brief Add a constraint to a problem in order to freeze the B-spline
    /// at its end.
    /// \param problem problem to which the constraint will be added.
    /// \param offset offset of the B-spline parameters in the problem's
    /// parameter list.
    template <typename P>
    void freezeCurveEnd (P& problem, size_type offset = 0) const throw ();

    /// \brief Regular spacing between B-spline knots. This is only valid for
    /// uniform B-splines.
    /// \return regular spacing between B-spline knots.
    value_type Dt () const ROBOPTIM_TRAJECTORY_DEPRECATED;

    /// \brief Translate the basis polynomials to a given time t1.
    /// \param t1 new center time, i.e. P = sum(a_i*(t-t1)^i, i={0,3})
    ///
    /// This method can be useful when one needs to have all the polynomials
    /// expressed in the same basis (e.g. t1 = 0).
    void translateBasisPolynomials (double t1) throw ();

    /// \brief Return the polynomial expression of the cubic B-spline on each
    /// time interval.
    void toPolynomials (std::vector<Polynomial3>& res) const throw ();

    /// \brief Constant getter for the basis polynomials of the cubic B-spline.
    /// \return constant reference to the basis polynomials.
    ///
    /// Note: computeBasisPolynomials() needs to be called beforehand (which is
    /// done in the CubicBSpline constructor).
    const std::vector <std::vector <Polynomial3> >&
    basisPolynomials() const
    {
      return basisPolynomials_;
    }

    /// \brief Get the number of control points of the spline.
    /// \return Number of control points of the spline.
    size_type getNumberControlPoints() const
    {
      return nbp_;
    }

    /// \brief Return the knot vector of the spline.
    /// \return knot vector of the spline (const).
    const std::vector <value_type>& knot_vector () const
    {
      return knots_;
    }

  protected:
    using Trajectory<3>::impl_compute;
    void impl_compute (result_t&, double) const throw ();
    void impl_derivative (gradient_t& g, double x, size_type order)
      const throw ();
    void impl_derivative (gradient_t& g, StableTimePoint, size_type order)
      const throw ();

    /// \brief Find the index of the interval in which t is.
    /// \param t instant considered.
    /// \return index of the interval in which t is.
    size_type interval (value_type t) const;

    /// \brief Compute the basis polynomials for the cubic B-spline.
    void computeBasisPolynomials ();

    /// \brief Compute the basis functions for a given instant t.
    /// \param t instant considered.
    /// \param order order of the basis functions.
    /// \return basis functions evaluated at t.
    vector_t basisFunctions (value_type t, size_type order) const
      ROBOPTIM_TRAJECTORY_DEPRECATED;

  private:
    /// \brief Number of control points.
    size_type nbp_;

    /// \brief Vector of knots.
    std::vector <value_type> knots_;

    /// \brief Basis polynomials.
    /// basisPolynomials_[i][j] = B_{i,i+j}
    std::vector <std::vector <Polynomial3> > basisPolynomials_;

    /// \brief Whether the B-spline is uniform.
    /// Note: for backward compatibility only.
    bool uniform_;
  };

  /// @}

} // end of namespace roboptim.

# include <roboptim/trajectory/cubic-b-spline.hxx>
#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HH
