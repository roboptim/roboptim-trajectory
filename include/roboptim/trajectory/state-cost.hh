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

#ifndef ROBOPTIM_TRAJECTORY_STATE_COST_HH
# define ROBOPTIM_TRAJECTORY_STATE_COST_HH
# include <boost/shared_ptr.hpp>

# include <roboptim/trajectory/fwd.hh>
# include <roboptim/core/derivable-function.hh>
# include <roboptim/trajectory/stable-time-point.hh>

namespace roboptim
{
  /// \addtogroup roboptim_meta_function
  /// @{

  /// \brief Cost function taking a state as its input.
  ///
  /// The state of a system is defined as the vector containing the
  /// configuration and derivatives up to order \f$r\f$ of the
  /// configuration:
  /**\f[
\textbf{X} =
\left(\textbf{q}, \frac{d\textbf{q}}{dt},\cdots,\frac{d^{r}\textbf{q}}{dt^{r}}(t)
\right)
     \f]*/
  /// \f$r\f$ is called the order of the state.
  /// The cost of a state is a real valued function defined over the state space:
  /// \f[C(\textbf{X})\in \textbf{R}\f]
  template <typename T>
  class StateCost : public DerivableFunction
  {
  public:
    /// \brief Trajectory type.
    typedef T trajectory_t;

    /// \brief Concrete class should call this constructor.
    StateCost (const trajectory_t&,
	       const DerivableFunction&,
	       const StableTimePoint tpt,
	       size_type order = 1) throw ();

    virtual ~StateCost () throw ();

    size_type order () const throw ();

    template <typename F, typename CLIST>
    static void addToProblem (const T& trajectory,
			      const DerivableFunction& function,
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
	    (new StateCost (trajectory, function, t * tMax, order));
	  problem.addConstraint (constraint, bounds);
	}
    }

  protected:
    void impl_compute (result_t&, const argument_t&) const throw ();
    void impl_gradient (gradient_t&, const argument_t&, size_type) const throw ();

  private:
    const trajectory_t& trajectory_;
    const DerivableFunction& function_;
    StableTimePoint tpt_;
    size_type order_;
  };

  /// @}

} // end of namespace roboptim.

# include <roboptim/trajectory/state-cost.hxx>
#endif //! ROBOPTIM_TRAJECTORY_STATE_COST_HH
