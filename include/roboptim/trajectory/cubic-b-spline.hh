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

#ifndef ROBOPTIM_TRAJECTORY_CUBIC_B_SPLINE_HH
# define ROBOPTIM_TRAJECTORY_CUBIC_B_SPLINE_HH
# include <roboptim/trajectory/sys.hh>

# include <roboptim/trajectory/trajectory.hh>
# include <roboptim/trajectory/fwd.hh>


namespace roboptim
{
  /// \addtogroup roboptim_function
  /// @{

  /** \brief Cubic B-Spline trajectory.
  
      Implement a B-Spline as a trajectory as described below: given 
      \li a number \f$m\geq 4\f$ of control points, 
      \li regularly spaced time points: \f$t_0 < t_1 < \cdots < t_{m}\f$, \f$\forall i\in\{0,...,m-1\}\f$, \f$t_{i+1}-t_i = \Delta t\f$
      \li \f$m\f$ control points \f$ P_0,\cdots, P_{m-1}\f$ in \f$\textbf{R}^n\f$ ,
      the cubic B-spline of control points \f$P_0,\cdots, P_{m-1}\f$ is defined over \f$[t_3,t_{m}]\f$ by
      \f[ B(t) = \sum_{i=0}^{m-1} P_i b_{i,3}(t) \f]
      where basis functions \f$b_{i,3}\f$ are defined by:
      \f{eqnarray*}{
      b_{i,3}(t)=& \frac{(t-t_i)^3}{6\Delta t^3} & \mbox{ if } t_{i} \leq t < t_{i+1} \\
      & \frac{(t-t_i)^2(t_{i+2}-t)+(t-t_i)(t_{i+3}-t)(t-t_{i+1})+(t_{i+4}-t)(t-t_{i+1})^2}{6\Delta t^3}& \mbox{ if } t_{i+1} \leq t < t_{i+2} \\  
      & \frac{(t-t_i)(t_{i+3}-t)^2+(t_{i+4}-t)(t-t_{i+1})(t_{i+3}-t)+(t_{i+4}-t)^2(t-t_{i+2})}{6\Delta t^3}& \mbox{ if } t_{i+2} \leq t < t_{i+3} \\  
      & \frac{(t_{i+4}-t)^3}{6\Delta t^3}& \mbox{ if } t_{i+3} \leq t < t_{i+4} 
      \f}
  */
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

  protected:
    void impl_compute (result_t&, double) const throw ();
    void impl_derivative (gradient_t& g, double x, size_type order)
      const throw ();
    void impl_derivative (gradient_t& g, StableTimePoint, size_type order)
      const throw ();

  private:
    /// \brief Number of control points.
    unsigned int nbp_;
  };

  /// @}

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HH
