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

#ifndef ROBOPTIM_TRAJECTORY_ANTHROPOMORPHIC_COST_FUNCTION_HH
# define ROBOPTIM_TRAJECTORY_ANTHROPOMORPHIC_COST_FUNCTION_HH
# include <boost/optional.hpp>

# include <roboptim/core/derivable-function.hh>
# include <roboptim/trajectory/fwd.hh>
# include <roboptim/trajectory/stable-time-point.hh>

namespace roboptim
{
  /// Cost function from ``An optimal control model unifying holonomic
  /// and nonholonomic walking'' Katja Mombaur, Jean-Paul Laumond,
  /// Eiichi Yoshida
  /// (2008 8th IEEE-RAS Interational Conference on Humanoid Robots).
  template <typename T>
  class AnthropomorphicCostFunction : public DerivableFunction
  {
  public:
    AnthropomorphicCostFunction (const T& trajectory,
				 const typename boost::optional<vector_t> alpha
				 = typename boost::optional<vector_t> (),
				 const typename boost::optional<value_type> ksi1
				 = typename boost::optional<value_type> (),
				 const typename boost::optional<value_type> ksi2
				 = typename boost::optional<value_type> ())
      throw ();
    ~AnthropomorphicCostFunction () throw ();

    static vector_t defaultAlpha ()
    {
      vector_t res (4);
      res[0] = 1.;
      res[1] = 10.;
      res[2] = 10.;
      res[3] = 5.;
      return res;
    }
    
    static value_type defaultKsi1 ()
    {
      return M_PI / 18.;
    }

    static value_type defaultKsi2 ()
    {
      return .5;
    }

  protected:
    void impl_compute (result_t& res, const argument_t& p) const throw ();
    void impl_gradient (gradient_t& grad, const argument_t& p, size_type i)
      const throw ();
  private:
    value_type alpha3 (value_type deltaTheta, value_type dsquare) const throw ();

    const T& trajectory_;
    const vector_t alpha_;
    value_type ksi1_;
    value_type ksi2_;
    value_type alpha3_;
  };
} // end of namespace roboptim.

# include <roboptim/trajectory/anthropomorphic-cost-function.hxx>
#endif //! ROBOPTIM_TRAJECTORY_ANTHROPOMORPHIC_COST_FUNCTION_HH
