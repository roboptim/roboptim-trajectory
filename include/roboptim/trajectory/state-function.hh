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

#ifndef ROBOPTIM_TRAJECTORY_STATE_FUNCTION_HH
# define ROBOPTIM_TRAJECTORY_STATE_FUNCTION_HH
# include <roboptim/trajectory/sys.hh>

# include <stdexcept>

# include <boost/shared_ptr.hpp>
# include <boost/static_assert.hpp>
# include <boost/type_traits/is_base_of.hpp>

# include <roboptim/trajectory/fwd.hh>
# include <roboptim/core/derivable-function.hh>
# include <roboptim/trajectory/stable-time-point.hh>

namespace roboptim
{
namespace trajectory
{
  /// \addtogroup roboptim_meta_function
  /// @{

  /// \brief Trajectory cost function defined by state evaluation at parameter.
  ///
  /// The state along a trajectory is defined as the vector containing the
  /// configuration and derivatives up to order \f$r\f$ of the
  /// configuration.
  /**\f[
\textbf{Cost}(\Gamma) = cost
\left({\Gamma(t)}, {\dot{\Gamma}(t)},\cdots,\frac{d^{r}\Gamma}{dt^{r}}(t)
\right)
     \f]*/
  /// where
  /// - \f$\textbf{Cost}\f$ is the trajectory cost,
  /// - \f$cost\f$ is the state cost,
  /// - \f$t\f$ is the parameter along the trajectory where the cost is
  /// evaluated (fixed at construction),
  /// - \f$r\f$ is called the order of the state.
  ///
  /// \tparam T trajectory type

  template <typename T>
  class StateFunction : public DerivableFunction
  {
    /// Check that T is a trajectory type.
    BOOST_STATIC_ASSERT((boost::is_base_of
			 <Trajectory<T::derivabilityOrder>, T>::value));
  public:
    /// \brief Trajectory type.
    typedef T trajectory_t;

    /// \brief Constructor.
    ///
    /// \param gamma Trajectory \f$\Gamma\f$ along which the state is evaluated.
    /// \param cost state cost: \f$cost\f$.
    /// \param tpt parameter \f$t\f$ where the state is evaluated.
    /// \param order order \f$r\f$ of derivation.
    /// \throw std::runtime_error
    StateFunction (const trajectory_t& gamma,
		   boost::shared_ptr<DerivableFunction> cost,
		   const StableTimePoint tpt,
		   size_type order = 1);

    virtual ~StateFunction ();

    size_type order () const;

    template <typename F, typename CLIST>
    static void addToProblem (const trajectory_t& trajectory,
			      boost::shared_ptr<DerivableFunction> function,
			      unsigned order,
			      Problem<F, CLIST>& problem,
			      typename Function::interval_t bounds,
			      unsigned nConstraints)
    {
      using namespace boost;

      for (unsigned i = 0; i < nConstraints; ++i)
	{
	  const value_type t = (i + 1.) / (nConstraints + 1.);
	  assert (t > 0. && t < 1.);
	  shared_ptr<DerivableFunction> constraint
	    (new StateFunction (trajectory, function, t * tMax, order));
	  problem.addConstraint (constraint, bounds);
	}
    }

    template <typename F, typename CLIST>
    static void addToProblem
    (const trajectory_t& trajectory,
     boost::shared_ptr<DerivableFunction> function,
     unsigned order,
     Problem<F, CLIST>& problem,
     const typename Function::intervals_t& bounds,
     const std::vector<typename Function::value_type>& scales,
     unsigned nConstraints)
    {
      using namespace boost;

      for (unsigned i = 0; i < nConstraints; ++i)
	{
	  const value_type t = (i + 1.) / (nConstraints + 1.);
	  assert (t > 0. && t < 1.);
	  shared_ptr<DerivableFunction> constraint
	    (new StateFunction (trajectory, function, t * tMax, order));
	  problem.addConstraint (constraint, bounds, scales);
	}
    }


  protected:
    void impl_compute (result_t&, const argument_t&) const;
    void impl_gradient (gradient_t&, const argument_t&, size_type)
      const;

  private:
    const trajectory_t& trajectory_;
    boost::shared_ptr<DerivableFunction> function_;
    StableTimePoint tpt_;
    size_type order_;
  };

  /// @}

} // end of namespace trajectory.
} // end of namespace roboptim.

# include <roboptim/trajectory/state-function.hxx>
#endif //! ROBOPTIM_TRAJECTORY_STATE_FUNCTION_HH
