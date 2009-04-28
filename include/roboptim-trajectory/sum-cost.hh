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

#ifndef ROBOPTIM_TRAJECTORY_SUM_COST_HH
# define ROBOPTIM_TRAJECTORY_SUM_COST_HH
# include <roboptim-trajectory/fwd.hh>
# include <roboptim-core/derivable-function.hh>

namespace roboptim
{
  template <typename T>
  class SumCost : TrajectoryCost<T>
  {
  public:
    typedef T trajectory_t;
    typedef StateCost<T> stateCost_t;

    SumCost (const trajectory_t&, const stateCost_t&, const vector_t&) throw ();
    //FIXME: implement automatic discretization?
    //SumCost (DiscreteRange)
    //Range = min + max + step

    virtual ~SumCost () throw ();

    virtual vector_t operator () (const vector_t&) const throw ();
    virtual gradient_t gradient (const vector_t&, int) const throw ();

  protected:
    stateCost_t stateCost_;
    const vector_t& points_;
  };

} // end of namespace roboptim.

# include <roboptim-trajectory/sum-cost.hxx>
#endif //! ROBOPTIM_TRAJECTORY_SUM_COST_HH
