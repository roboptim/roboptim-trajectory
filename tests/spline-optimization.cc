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

#include <boost/numeric/ublas/io.hpp>

#include <roboptim/core/solver-factory.hh>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/spline.hh>
#include <roboptim/trajectory/trajectory-sum-cost.hh>
#include <roboptim/trajectory/state-cost.hh>

#include "common.hh"

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;

typedef boost::variant<const DerivableFunction*,
		       const LinearFunction*> constraint_t;
typedef Solver<DerivableFunction, constraint_t> solver_t;


struct MyStateCost : public StateCost<Spline>
{
  MyStateCost (size_type n)
    : StateCost<Spline> (n)
  {
  }

  virtual vector_t operator () (const vector_t&) const throw ()
  {
    vector_t res (m);
    res.clear ();
    return res;
  }

  virtual gradient_t gradient (const vector_t&, int) const throw ()
  {
    gradient_t grad (n);
    grad.clear ();
    return grad;
  }
};

int run_test ()
{
  Spline::vector_t params (8);

  // Initial position.
  params[0] = 0.,  params[1] = 0.;
  // Control point 1.
  params[2] = 25.,  params[3] = 50.;
  // Control point 2.
  params[4] = 50.,  params[5] = 25.;
  // Final position.
  params[6] = 100., params[7] = 100.;

  Spline spline (std::make_pair (0., 5.), 2, params);
  discreteInterval_t interval (0., 5., 0.01);

  std::cout
    << "# Values:" << std::endl
    << "# " << spline (0.) << std::endl
    << "# " << spline (2.5) << std::endl
    << "# " << spline (5.) << std::endl

    << "# 1st derivative:" << std::endl
    << "# " << spline.derivative (0., 1) << std::endl
    << "# " << spline.derivative (2.5, 1) << std::endl
    << "# " << spline.derivative (5., 1) << std::endl

    << "# 2nd derivative:" << std::endl
    << "# " << spline.derivative (0., 2) << std::endl
    << "# " << spline.derivative (2.5, 2) << std::endl
    << "# " << spline.derivative (5., 2) << std::endl

    << "# variationConfigWrtParam:" << std::endl
    << "# " << spline.variationConfigWrtParam (0.) << std::endl
    << "# " << spline.variationConfigWrtParam (2.5) << std::endl
    << "# " << spline.variationConfigWrtParam (5.) << std::endl;

  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();
  gnuplot << plot_xy (spline, interval);

  // Optimize.
  discreteInterval_t costInterval (0., 5., 0.5);

  MyStateCost statecost (2);
  TrajectorySumCost<Spline> cost (spline, statecost, costInterval);

  solver_t::problem_t problem (cost);

  SolverFactory<solver_t> factory ("cfsqp", problem);
  solver_t& solver = factory ();

  solver_t::result_t res = solver.minimum ();
  Result& result = boost::get<Result> (res);

  spline.setParameters (result.x);

  std::cout << (gnuplot << plot_xy (spline, interval));
  return 0;
}

GENERATE_TEST ()
