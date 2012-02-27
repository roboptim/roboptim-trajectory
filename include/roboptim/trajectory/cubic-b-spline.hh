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

    template <typename P>
    void freezeCurveStart (P& problem, size_type offset = 0) const throw ();
    template <typename P>
    void freezeCurveEnd (P& problem, size_type offset = 0) const throw ();

    value_type Dt () const ROBOPTIM_TRAJECTORY_DEPRECATED;

  protected:
    using Trajectory<3>::impl_compute;
    void impl_compute (result_t&, double) const throw ();
    void impl_derivative (gradient_t& g, double x, size_type order)
      const throw ();
    void impl_derivative (gradient_t& g, StableTimePoint, size_type order)
      const throw ();

    size_type interval (value_type t) const;
    vector_t basisFunctions (value_type t, size_type order) const
      ROBOPTIM_TRAJECTORY_DEPRECATED;
    void computeBasisPolynomials ();

  private:
    /// \brief Number of control points.
    size_type nbp_;
    /// Vector of knots
    std::vector <value_type> knots_;
    /// basisPolynomials_[i][j] = B_{i,i+j}
    std::vector <std::vector <Polynomial3> > basisPolynomials_;
    /// For backward compatibility only
    bool uniform_;
  };

  /// @}

} // end of namespace roboptim.

# include <roboptim/trajectory/cubic-b-spline.hxx>
#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HH
