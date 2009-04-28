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
# include <roboptim-trajectory/fwd.hh>
# include <roboptim-core/derivable-function.hh>

namespace roboptim
{

  /**
     \brief Cost of a state

     The state of a system is defined as the vector containing the configuration and derivatives up to order \f$r\f$ of the configuration:
     \f[
     \textbf{X} = \left(\textbf{q}, \frac{d\textbf{q}}{dt},\cdots,\frac{d^{r}\textbf{q}}{dt^{r}}(t)\right)
     \f]
     \f$r\f$ is called the order of the state.

     The cost of a state is a real valued function defined over the state space:
     \f[
     C(\textbf{X})\in \textbf{R}
     \f]
  */
  template <typename T>
  class StateCost : public DerivableFunction
  {
  public:
    typedef T trajectory_t;
    StateCost (size_type m);

    virtual ~StateCost();
  };
} // end of namespace roboptim.

# include <roboptim-trajectory/state-cost.hxx>
#endif //! ROBOPTIM_TRAJECTORY_STATE_COST_HH
