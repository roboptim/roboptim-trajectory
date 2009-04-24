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

/**
 * \brief Class Trajectory declaration.
 */

#ifndef ROBOPTIM_TRAJECTORY_TRAJECTORY_HH
# define ROBOPTIM_TRAJECTORY_TRAJECTORY_HH
# include <utility>

# include <roboptim-trajectory/fwd.hh>
# include <roboptim-core/n-times-derivable-function.hh>

namespace roboptim
{
  /**
     \brief Abstract trajectory
     A trajectory is a piecewise smooth mapping \f$\Gamma\f$ from
     <ul>
     <li>the Cartesian product of a definition interval and a vector space
     of parameters \f$\textbf{R}^m\f$</li>
     <li>to a vector space \f$\textbf{R}^n\f$:
     \f[
     \begin{array}{llll}
     \Gamma : & [t_{min}, t_{max}] \times \textbf{R}^m & \rightarrow &
     \textbf{R}^n \\
     & (t, \textbf{p}) & \rightarrow & \Gamma_{\textbf{p}}(t)
     \end{array}
     \f]
     </li>
     </ul>
  */
  template <unsigned DerivabilityOrder>
  class Trajectory : public NTimesDerivableFunction<DerivabilityOrder>
  {
  public:
    typedef NTimesDerivableFunction<DerivabilityOrder> parent_t;

    typedef typename parent_t::value_type value_type;
    typedef typename parent_t::size_type size_type;
    typedef typename parent_t::vector_t vector_t;
    typedef typename parent_t::jacobian_t jacobian_t;

    typedef std::pair<value_type, value_type> bound_t;


    Trajectory (bound_t, size_type, const vector_t&) throw ();
    virtual ~Trajectory () throw ();

    /// \name Accessing parameters, and state.
    /// \{

    vector_t& parameters () throw ();
    const vector_t& parameters () const throw ();

    bound_t timeRange () const throw ();
    value_type length () const throw ();

    /// \brief Get state along trajectory
    ///
    /// \param t value \f$t\f$ in the definition interval.
    /// \param order the higher order \f$r\f$ of the required derivative
    /// \return the state defined as the vector containing the
    /// config and first derivatives:
    /// \f[
    /// \textbf{X}(t) = \left(\Gamma_{\textbf{p}}(t),
    /// \frac{d\Gamma_{\textbf{p}}}{dt}(t),\cdots,
    /// \frac{d^{r}\Gamma_{\textbf{p}}}{dt^{r}}(t)\right)
    /// \f]
    /// The configuration and derivatives are concatenated in one vector.
    virtual vector_t state (double t, size_type order) const throw ();

    /// \}

    /// \name Accessing parameters and gradients
    /// \{

    /// \brief Get the variation of a configuration with respect to parameter
    /// vector.
    /// \param t value \f$t\f$ in the definition interval.
    /// \return Jacobian:
    /// \f[\frac{\partial\Gamma_{\textbf{p}}(t)}{\partial\textbf{p}}\f]
    virtual jacobian_t variationConfigWrtParam (double t) const throw () = 0;

    /// \brief Get the variation of a derivative with respect to parameter
    /// vector.
    ///
    /// \param t value \f$t\f$ in the definition interval.
    /// \param order order \f$r\f$ of the derivative.
    /// \return Jacobian:
    /// \f[\frac{\partial}{\partial\textbf{p}}
    /// \left(\frac{d^r\Gamma_{\textbf{p}}}{dt^r}(t)\right)\f]
    virtual jacobian_t variationDerivWrtParam (double t, size_type order)
      const throw () = 0;

    /**
       \brief Get the variation of the state with respect to parameter vector

       \param t value \f$t\f$ in the definition interval.
       \param order order \f$r\f$ of the derivative.
       \return Jacobian:
       \f[\left(\begin{array}{c}\frac{\partial}{\partial \lambda}
       \Gamma_{\textbf{p}}(t) \						\
       \vdots \								\
       \frac{\partial}{\partial \lambda}\frac{d\Gamma_{\textbf{p}}^r}{dt^r}(t)
       \end{array}\right)\f]
    */
    jacobian_t variationStateWrtParam (double t, size_type order) const throw ();

    /// \}


    /// \name Singular points
    /// \{

    /// \brief Get number of singular points
    size_type singularPoints () const throw ();

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

    virtual std::ostream& print (std::ostream&) const throw ();
  protected:
    bound_t timeRange_;
    vector_t parameters_;
    size_type singularPoints_;
  };
} // end of namespace roboptim.

# include <roboptim-trajectory/trajectory.hxx>
#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HH
