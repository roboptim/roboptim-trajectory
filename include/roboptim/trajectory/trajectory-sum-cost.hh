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
# include <roboptim/trajectory/sys.hh>

# include <boost/shared_ptr.hpp>
# include <boost/tuple/tuple.hpp>

# include <roboptim/trajectory/fwd.hh>
# include <roboptim/core/derivable-function.hh>
# include <roboptim/trajectory/stable-time-point.hh>

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


    /// \name Discrete interval of stable time points.
    /// \{

    /// \brief Types representing a discrete interval of stable time points.
    /// A discrete interval is a triplet of values:
    /// - lower bound,
    /// - upper bound,
    /// - step.
    typedef boost::tuple<StableTimePoint,
			 StableTimePoint,
			 StableTimePoint> discreteStableTimePointInterval_t;


    /// \brief Construct a discrete interval.
    ///
    /// \param min miminum value of the interval
    /// \param max maxinum value of the interval
    /// \param step discretization step
    static discreteStableTimePointInterval_t
    makeDiscreteInterval (StableTimePoint min,
			  StableTimePoint max,
			  StableTimePoint step)
    {
      return discreteStableTimePointInterval_t (min, max, step);
    }

    /// \brief Construct a discrete interval.
    ///
    /// \param interval continuous interval
    /// \param step discretization step
    static discreteStableTimePointInterval_t
    makeDiscreteInterval (interval_t interval,
			  StableTimePoint step)
    {
      return discreteStableTimePointInterval_t (getLowerBound (interval),
						getUpperBound (interval),
						step);
    }

    /// \brief Get the lower bound of a discrete interval
    ///
    /// \param interval accessed discrete interval
    /// \return lower bound of the discrete interval
    static StableTimePoint
    getLowerBound (const discreteStableTimePointInterval_t& interval) throw ()
    {
      return boost::get<0> (interval);
    }

    /// \brief Get the upper bound of a discrete interval
    ///
    /// \param interval accessed discrete interval
    /// \return upper bound of the discrete interval
    static StableTimePoint
    getUpperBound (const discreteStableTimePointInterval_t& interval) throw ()
    {
      return boost::get<1> (interval);
    }

    /// \brief Get the upper step of a discrete interval
    ///
    /// \param interval accessed discrete interval
    /// \return upper step of the discrete interval
    static StableTimePoint
    getStep (const discreteStableTimePointInterval_t& interval) throw ()
    {
      return boost::get<2> (interval);
    }

    /// \}



    /// \brief Constructor.
    ///
    /// \param gamma Trajectory \f$\Gamma\f$ along which the state is evaluated.
    /// \param cost state cost: \f$cost\f$.
    /// \param interval Stable time points list where cost will be evaluated.
    /// \param order order \f$r\f$ of derivation.
    TrajectorySumCost (const trajectory_t& gamma,
		       boost::shared_ptr<DerivableFunction> cost,
		       const discreteStableTimePointInterval_t& interval,
		       size_type order = 1) throw ();

    virtual ~TrajectorySumCost () throw ();

    size_type order () const throw ();

  protected:
    void impl_compute (result_t&, const argument_t&) const throw ();
    void impl_gradient (gradient_t&, const argument_t&, size_type) const throw ();

  private:
    const trajectory_t& trajectory_;
    boost::shared_ptr<DerivableFunction> function_;
    discreteStableTimePointInterval_t interval_;
    size_type order_;
  };

  /// @}

} // end of namespace roboptim.

# include <roboptim/trajectory/trajectory-sum-cost.hxx>
#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_SUM_COST_HH
