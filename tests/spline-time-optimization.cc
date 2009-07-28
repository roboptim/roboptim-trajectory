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

#include <fstream>

#include <boost/assign/list_of.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <roboptim/core/finite-difference-gradient.hh>
#include <roboptim/core/solver-factory.hh>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

#include <roboptim/trajectory/free-time-trajectory.hh>
#include <roboptim/trajectory/freeze.hh>
#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/limit-speed.hh>
#include <roboptim/trajectory/spline-length.hh>
#include <roboptim/trajectory/spline.hh>
#include <roboptim/trajectory/trajectory-cost.hh>

#include <roboptim/trajectory/visualization/limit-speed.hh>


#include "common.hh"

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;

typedef boost::mpl::vector<DerivableFunction, LinearFunction> constraint_t;
typedef Solver<DerivableFunction, constraint_t> solver_t;
typedef FreeTimeTrajectory<Spline::derivabilityOrder> freeTime_t;


// Problem parameters.
const unsigned nControlPoints = 31;
const double vMax = 200.;

int run_test ()
{
  using namespace boost;
  using namespace boost::assign;

  const double finalPos = 200.;
  Spline::vector_t params (nControlPoints);

  for (unsigned i = 0; i < nControlPoints; ++i)
    params[i] = finalPos / (nControlPoints - 1) * i;

  // Make trajectories.
  Spline::interval_t timeRange = Spline::makeInterval (0., 4.);
  Spline spline (timeRange, 1, params, "before");
  freeTime_t freeTimeTraj (spline, 1.);

  // Define cost.
  Function::matrix_t a (1, freeTimeTraj.parameters ().size ());
  a.clear ();
  a (0, 0) = 1.;
  Function::vector_t b (1);
  b.clear ();
  roboptim::NumericLinearFunction cost (a, b);

  // Create problem.
  solver_t::problem_t problem (cost);
  problem.startingPoint () = freeTimeTraj.parameters ();

  {
    // Be sure that scale is positive.
    Function::matrix_t a (1, problem.function ().inputSize ());
    Function::vector_t b (1);
    a.clear (), b.clear ();
    a (0, 0) = 1.;
    shared_ptr<NumericLinearFunction> speedPositivity
      (new NumericLinearFunction (a, b));
    problem.addConstraint
      (static_pointer_cast<LinearFunction> (speedPositivity),
       Function::makeLowerInterval (1e-3));
  }

  const freeTime_t::vector_t freeTimeParams = freeTimeTraj.parameters ();
  const unsigned freeTimeParamsSize = freeTimeParams.size ();

  typedef Freeze<DerivableFunction, constraint_t, DerivableFunction> freeze_t;
  freeze_t freeze (problem,
		   list_of <freeze_t::frozenArgument_t>
		   (1, freeTimeParams[1])
		   (freeTimeParamsSize - 1,
		    freeTimeParams[freeTimeParamsSize - 1]));
  freeze ();

  Function::interval_t vRange (0., 2 * vMax * vMax);
  LimitSpeed<FreeTimeTrajectory<Spline::derivabilityOrder> >::addToProblem
    (freeTimeTraj, problem, vRange, nControlPoints * 10);

  std::ofstream limitSpeedStream ("limit-speed.gp");
  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();

  gnuplot
    << set ("multiplot layout 1,2 title "
	    "'variation of speed before and after optimization'")
    << set ("grid");
  gnuplot << plot_limitSpeed (freeTimeTraj, vMax);


  SolverFactory<solver_t> factory ("cfsqp", problem);
  solver_t& solver = factory ();

  std::cout << solver << std::endl;

  solver_t::result_t res = solver.minimum ();
  std::cout << res << std::endl;

  FreeTimeTrajectory<Spline::derivabilityOrder> optimizedTrajectory =
    freeTimeTraj;

  switch (solver.minimumType ())
    {
    case GenericSolver::SOLVER_VALUE:
      {
	const Result& result = solver.getMinimum<Result> ();
	optimizedTrajectory.setParameters (result.x);
	break;
      }

    case GenericSolver::SOLVER_VALUE_WARNINGS:
      {
	const ResultWithWarnings& result =
	  solver.getMinimum<ResultWithWarnings> ();
	optimizedTrajectory.setParameters (result.x);
	break;
      }

    case GenericSolver::SOLVER_NO_SOLUTION:
    case GenericSolver::SOLVER_ERROR:
      return 1;
    }

  gnuplot << plot_limitSpeed (optimizedTrajectory, vMax);
  limitSpeedStream << (gnuplot << unset ("multiplot"));
  return 0;
}

GENERATE_TEST ()
