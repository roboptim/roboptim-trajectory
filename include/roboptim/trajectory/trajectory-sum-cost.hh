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

#ifndef ROBOPTIM_TRAJECTORY_TRAJECTORY_SUM_COST_HH
# define ROBOPTIM_TRAJECTORY_TRAJECTORY_SUM_COST_HH
# include <boost/shared_ptr.hpp>

# include <roboptim/trajectory/fwd.hh>
# include <roboptim/core/derivable-function.hh>

namespace roboptim
{
  /// \addtogroup roboptim_meta_function
  /// @{

  /// \brief Trajectory cost function defined by sum of state evaluations at parameters.
  ///
  /// The state along a trajectory is defined as the vector containing the
  /// configuration and derivatives up to order \f$r\f$ of the
  /// configuration.
  /**\f[
\textbf{Cost}(\Gamma) = \sum_{i=1}{n} cost
\left({\Gamma(t_i)}, {\dot{\Gamma}(t_i)},\cdots,\frac{d^{r}\Gamma}{dt^{r}}(t_i)
\right)
     \f]*/
  /// where
  /// - \f$\textbf{Cost}\f$ is the trajectory cost,
  /// - \f$cost\f$ is the state cost,
  /// - \f$t_i\ \ i=1\cdots n\f$ are the parameters along the trajectory where the cost is evaluated,
  /// - \f$r\f$ is called the order of the state.
  ///
  /// \tparam T trajectory type
  template <typename T>
  class TrajectorySumCost : public DerivableFunction
  {
  public:
    /// \brief Parent type.
    typedef DerivableFunction parent_t;
    /// \brief Trajectory type.
    typedef T trajectory_t;

    /// \brief Import vector type.
    typedef typename parent_t::vector_t vector_t;
    /// \brief Import gradient type.
    typedef typename parent_t::gradient_t gradient_t;
    /// \brief Import discrete interval type.
    typedef typename parent_t::discreteInterval_t discreteInterval_t;

    /// \brief Constructor.
    ///
    /// \param gamma Trajectory \f$\Gamma\f$ along which the state is evaluated.
    /// \param cost state cost: \f$cost\f$.
    /// \param tpt parameter \f$t\f$ where the state is evaluated.
    /// \param order order \f$r\f$ of derivation.
    TrajectorySumCost (const trajectory_t& gamma,
		       boost::shared_ptr<DerivableFunction> cost,
		       const discreteInterval_t& interval,
		       size_type order = 1) throw ();

    virtual ~TrajectorySumCost () throw ();

    size_type order () const throw ();

  protected:
    void impl_compute (result_t&, const argument_t&) const throw ();
    void impl_gradient (gradient_t&, const argument_t&, size_type) const throw ();

  private:
    const trajectory_t& trajectory_;
    boost::shared_ptr<DerivableFunction> function_;
    discreteInterval_t interval_;
    size_type order_;
  };

  /// @}

} // end of namespace roboptim.

# include <roboptim/trajectory/trajectory-sum-cost.hxx>
#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_SUM_COST_HH
