// Copyright (C) 2009 by Florent Lamiraux, Thomas Moulard, AIST, CNRS, INRIA.
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

#ifndef ROBOPTIM_TRAJECTORY_TRAJECTORY_HH
# define ROBOPTIM_TRAJECTORY_TRAJECTORY_HH
# include <roboptim/core/sys.hh>

# include <stdexcept>
# include <utility>

# include <roboptim/core/n-times-derivable-function.hh>
# include <roboptim/trajectory/fwd.hh>
# include <roboptim/trajectory/stable-time-point.hh>

# define ROBOPTIM_IMPLEMENT_CLONE(C)		\
  virtual C* clone () const			\
  {						\
    return new C (*this);			\
  }

namespace roboptim
{
namespace trajectory
{
  namespace detail
  {
    template <typename T>
    double fixTime (double t, const T& trajectory);

    template <typename T>
    double fixTime (double t, const T& trajectory)
    {
      double epsilon = trajectory.tolerance ();
      const double tmin = Function::getLowerBound (trajectory.timeRange ());
      const double tmax = Function::getUpperBound (trajectory.timeRange ());
      if (t < tmin && (t + epsilon) >= tmin)
        return tmin;
      if (t > tmax && (t - epsilon) <= tmax)
        return tmax;
      return t;
    }
  } // end of namespace detail.

  /// \addtogroup roboptim_meta_function
  /// @{

  /// \brief Abstract trajectory
  ///
  /// A trajectory is a piecewise smooth mapping \f$\Gamma\f$ from
  /// - the Cartesian product of a definition interval and a vector space
  ///   of parameters \f$\textbf{R}^m\f$
  /// - to a vector space \f$\textbf{R}^n\f$:
  /**   \f[
        \begin{array}{llll}
        \Gamma: & [t_{min}, t_{max}] \times \textbf{R}^m & \rightarrow & \textbf{R}^n \\
        & (t, \textbf{p}) & \rightarrow & \Gamma_{\textbf{p}}(t)
        \end{array}
        \f]*/
  ///
  /// \tparam DerivabilityOrder derivability order
  template <unsigned DerivabilityOrder>
  class Trajectory : public NTimesDerivableFunction<DerivabilityOrder>
  {
  public:
    /// \brief Parent type.
    typedef NTimesDerivableFunction<DerivabilityOrder> parent_t;

    using parent_t::operator ();
    using parent_t::derivative;
    using parent_t::impl_compute;
    using parent_t::impl_derivative;

    /// \brief Import value type.
    typedef typename parent_t::value_type value_type;
    /// \brief Import size type.
    typedef typename parent_t::size_type size_type;
    /// \brief Import vector type.
    typedef typename parent_t::vector_t         vector_t;
    typedef typename parent_t::vector_ref       vector_ref;
    typedef typename parent_t::const_vector_ref const_vector_ref;
    /// \brief Import result type.
    typedef typename parent_t::result_t         result_t;
    typedef typename parent_t::result_ref       result_ref;
    typedef typename parent_t::const_result_ref const_result_ref;
    /// \brief Import gradient type.
    typedef typename parent_t::gradient_t         gradient_t;
    typedef typename parent_t::gradient_ref       gradient_ref;
    typedef typename parent_t::const_gradient_ref const_gradient_ref;
    /// \brief Import jacobian type.
    typedef typename parent_t::jacobian_t jacobian_t;
    /// \brief Import interval type.
    typedef typename parent_t::interval_t interval_t;

    virtual ~Trajectory ();

    /// \name Accessing parameters, and state.
    /// \{

    const_vector_ref parameters () const;

    /// \brief Set parameters.
    /// \param vector_t parameters.
    /// \throw std::runtime_error
    virtual void setParameters (const_vector_ref);

    interval_t timeRange () const;
    value_type length () const;

    /// \brief Get state along trajectory
    ///
    /// \param t value \f$t\f$ in the definition interval.
    /// \param order the higher order \f$r\f$ of the required derivative
    /// \return the state defined as the vector containing the
    /// config and first derivatives:
    /** \f[
        \textbf{X}(t) = \left(\Gamma_{\textbf{p}}(t),
        \frac{d\Gamma_{\textbf{p}}}{dt}(t),\cdots,
        \frac{d^{r}\Gamma_{\textbf{p}}}{dt^{r}}(t)\right)
        \f]*/
    /// The configuration and derivatives are concatenated in one vector.
    virtual vector_t state (double t, size_type order) const;
    virtual vector_t state (StableTimePoint t, size_type order) const;

    /// \}

    /// \name Accessing parameters and gradients
    /// \{

    /// \brief Get the variation of a configuration with respect to parameter
    /// vector.
    /// \param t value \f$t\f$ in the definition interval.
    /// \return Jacobian:
    /// \f[\frac{\partial\Gamma_{\textbf{p}}(t)}{\partial\textbf{p}}\f]
    virtual jacobian_t variationConfigWrtParam (double t) const = 0;

    /// \brief Get the variation of a derivative with respect to parameter
    /// vector.
    ///
    /// \param t value \f$t\f$ in the definition interval.
    /// \param order order \f$r\f$ of the derivative.
    /// \return jacobian
    /** \f[
        \frac{\partial}{\partial\textbf{p}}
        \left(\frac{d^r\Gamma_{\textbf{p}}}{dt^r}(t)\right)
        \f]*/
    virtual jacobian_t variationDerivWrtParam (double t, size_type order)
      const = 0;

    /// \brief Get the variation of the state with respect to parameter vector
    ///
    /// \param t value \f$t\f$ in the definition interval.
    /// \param order order \f$r\f$ of the derivative.
    /// \return jacobian
    /** \f[
        \left(\begin{array}{c}\frac{\partial}{\partial \lambda}
        \Gamma_{\textbf{p}}(t) \\
        \vdots \\
        \frac{\partial}{\partial \lambda}\frac{d\Gamma_{\textbf{p}}^r}{dt^r}(t)
        \end{array}\right)
        \f]**/

    jacobian_t variationStateWrtParam (double t, size_type order) const;
    jacobian_t variationStateWrtParam
    (StableTimePoint stp, size_type order) const;

    /// \}


    /// \name Singular points
    /// \{

    /// \brief Get number of singular points
    size_type singularPoints () const;

    /// \brief Get singular point at given rank.
    virtual value_type singularPointAtRank (size_type rank) const = 0;

    /// \brief Get left limit value of derivative at given singular point
    ///
    /// \param rank rank of the singular points.
    /// \param order order of derivation.
    /// \return Limit of the derivative at singular point
    /// for increasing parameter values.
    virtual vector_t derivBeforeSingularPoint (size_type rank, size_type order)
      const = 0;

    /// \brief Get right limit value of derivative at given singular point
    /// \param rank rank of the singular points.
    /// \param order order of derivation.
    /// \retval derivative Limit of the derivative at singular point for
    /// decreasing parameter values.
    virtual vector_t derivAfterSingularPoint (size_type rank, size_type order)
      const = 0;
    /// \}


    result_t operator () (StableTimePoint argument) const
    {
      result_t result (this->outputSize ());
      result.setZero();
      (*this) (result, argument);
      return result;
    }

    void operator () (result_ref result, StableTimePoint argument) const
    {
      assert (this->isValidResult (result));
      this->impl_compute (result, argument);
      assert (this->isValidResult (result));
    }

    gradient_t derivative (StableTimePoint argument, size_type order = 1) const
    {
      gradient_t derivative (this->derivativeSize ());
      derivative.setZero ();
      this->derivative (derivative, argument, order);
      return derivative;
    }

    void derivative (gradient_ref derivative,
		     StableTimePoint argument,
		     size_type order = 1) const
    {
      assert (order <= Trajectory<DerivabilityOrder>::derivabilityOrder
	      && this->isValidDerivative (derivative));
      this->impl_derivative (derivative, argument, order);
      assert (this->isValidDerivative (derivative));
    }

    virtual jacobian_t
    variationConfigWrtParam (StableTimePoint tp)
      const = 0;
    virtual jacobian_t
    variationDerivWrtParam (StableTimePoint tp, size_type order)
      const = 0;

    bool isValidTime (value_type t) const;

    /// \brief Normalize angles in parameters array.
    ///
    /// Make sure angles are continuous.
    /// \param index Angles index in parameter array.
    virtual void normalizeAngles (size_type index);

    virtual Trajectory<DerivabilityOrder>* clone () const = 0;

    /// \brief Clone and resize a trajectory.
    virtual Trajectory<DerivabilityOrder>* resize (interval_t timeRange)
      const = 0;

    /// \name Tolerance for inclusion of parameter in interval of definition.
    /// \{
    /// Set tolerance for inclusion of parameter in interval of definition.
    void tolerance (const double& tolerance);
    /// Get tolerance for inclusion of parameter in interval of definition.
    double tolerance () const;
    /// \}

    virtual std::ostream& print (std::ostream&) const;
  protected:
    void impl_compute (result_ref, StableTimePoint) const;
    virtual void
    impl_derivative (gradient_ref g, StableTimePoint, size_type order)
      const = 0;

    Trajectory (interval_t, size_type, const_vector_ref,
		std::string name = std::string ());

    /// \brief Internal version of normalizeAngles allowing an optional offset.
    ///
    /// Used to factorize code between trajectories and free time trajectories.
    /// \param index Angles index in parameter array.
    /// \param offset Index of the first control point in the parameter vector.
    virtual void normalizeAngles (size_type index, size_type offset);

    interval_t timeRange_;
    vector_t parameters_;
    size_type singularPoints_;
  private:
    double tolerance_;
  };

  /// @}

} // end of namespace trajectory.
} // end of namespace roboptim.

# include <roboptim/trajectory/trajectory.hxx>
#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HH
