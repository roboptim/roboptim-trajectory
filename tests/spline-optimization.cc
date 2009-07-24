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

#include <boost/assign/list_of.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <roboptim/core/finite-difference-gradient.hh>
#include <roboptim/core/solver-factory.hh>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

#include <roboptim/trajectory/freeze.hh>
#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/spline-length.hh>
#include <roboptim/trajectory/spline.hh>
#include <roboptim/trajectory/trajectory-cost.hh>
#include <roboptim/trajectory/visualization/trajectory.hh>

#include "common.hh"

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;

typedef boost::mpl::vector<DerivableFunction, LinearFunction> constraint_t;
typedef Solver<DerivableFunction, constraint_t> solver_t;

int run_test ()
{
  using namespace boost::assign;
  Spline::vector_t params (8);

  // Initial position.
  params[0] = 0.,  params[1] = 0.;
  // Control point 1.
  params[2] = 25.,  params[3] = 100.;
  // Control point 2.
  params[4] = 75.,  params[5] = 0.;
  // Final position.
  params[6] = 100., params[7] = 100.;

  Spline::interval_t timeRange = Spline::makeInterval (0., 4.);

  Spline spline (timeRange, 2, params, "before");
  discreteInterval_t interval (0., 4., 0.01);

  std::cout
    << "# Values:" << std::endl
    << "# " << spline (0.) << std::endl
    << "# " << spline (2.5) << std::endl
    << "# " << spline (4.) << std::endl

    << "# 1st derivative:" << std::endl
    << "# " << spline.derivative (0., 1) << std::endl
    << "# " << spline.derivative (2.5, 1) << std::endl
    << "# " << spline.derivative (4., 1) << std::endl

    << "# 2nd derivative:" << std::endl
    << "# " << spline.derivative (0., 2) << std::endl
    << "# " << spline.derivative (2.5, 2) << std::endl
    << "# " << spline.derivative (4., 2) << std::endl

    << "# variationConfigWrtParam:" << std::endl
    << "# " << spline.variationConfigWrtParam (0.) << std::endl
    << "# " << spline.variationConfigWrtParam (2.5) << std::endl
    << "# " << spline.variationConfigWrtParam (4.) << std::endl;

  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();
  gnuplot
    << set ("multiplot layout 1,2")
    << set ("grid")
    << plot_xy (spline);

  // Optimize.
  discreteInterval_t costInterval (0., 4., 0.5);
  SplineLength cost (spline, costInterval);

  // Check cost gradient.
  try
  {
    Function::vector_t x (params.size ());
    x.clear ();
    checkGradientAndThrow (cost, 0, x, 2e-4);

    x = params;
    checkGradientAndThrow (cost, 0, x, 2e-4);
  }
  catch (BadGradient& bg)
    {
      std::cout << bg << std::endl;
      return 1;
    }

  solver_t::problem_t problem (cost);
  problem.startingPoint () = params;

  typedef Freeze<DerivableFunction, constraint_t, LinearFunction> freeze_t;
  freeze_t freeze (problem,
		   list_of <freeze_t::frozenArgument_t>
		   (0, params[0])
		   (1, params[1])
		   (params.size () - 2, params[params.size () - 2])
		   (params.size () - 1, params[params.size () - 1]));
  freeze ();

  SolverFactory<solver_t> factory ("cfsqp", problem);
  solver_t& solver = factory ();

  std::cerr << "Cost function (before): " << cost (params) << std::endl;
  std::cerr << "Parameters (before): " << params << std::endl;

  solver_t::result_t res = solver.minimum ();

  std::cerr << solver << std::endl;

  switch (res.which ())
    {
    case GenericSolver::SOLVER_VALUE:
      {
	Result& result = boost::get<Result> (res);
	Spline optimizedSpline (timeRange, 2, result.x, "after");
	params = result.x;
	gnuplot << plot_xy (optimizedSpline);
	break;
      }

    case GenericSolver::SOLVER_NO_SOLUTION:
      {
	std::cerr << "No solution" << std::endl;
	return 1;
      }
    case GenericSolver::SOLVER_VALUE_WARNINGS:
      {
	ResultWithWarnings& result = boost::get<ResultWithWarnings> (res);
	Spline optimizedSpline (timeRange, 2, result.x, "after");
	params = result.x;
	std::cerr << result << std::endl;
	gnuplot << plot_xy (optimizedSpline);
	break;
      }

    case GenericSolver::SOLVER_ERROR:
      {
	SolverError& result = boost::get<SolverError> (res);
	std::cerr << result << std::endl;
      return 1;
      }
    }

  std::cerr << "Parameters (after): " << params << std::endl;

  // Check cost gradient (final).
  try
    {
      checkGradientAndThrow (cost, 0, params, 2e-3);
    }
  catch (BadGradient& bg)
    {
      std::cout << bg << std::endl;
      return 1;
    }

  std::cout << (gnuplot << unset ("multiplot"));
  return 0;
}

GENERATE_TEST ()
