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

#undef NDEBUG

#include <roboptim/trajectory/sys.hh>

#include <fstream>

#include <boost/assign/list_of.hpp>
#include <boost/mpl/vector.hpp>

#include <roboptim/core/io.hh>
#include <roboptim/core/finite-difference-gradient.hh>
#include <roboptim/core/solver-factory.hh>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

#include <roboptim/trajectory/free-time-trajectory.hh>
#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/limit-speed.hh>
#include <roboptim/trajectory/spline-length.hh>
#include <roboptim/trajectory/cubic-b-spline.hh>
#include <roboptim/trajectory/trajectory-cost.hh>

#include <roboptim/trajectory/visualization/limit-speed.hh>


#include <roboptim/core/plugin/cfsqp.hh>


#include "shared-tests/common.hh"

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;

typedef CFSQPSolver::problem_t::constraints_t constraint_t;
typedef CFSQPSolver solver_t;
typedef FreeTimeTrajectory<CubicBSpline> freeTime_t;


// Problem parameters.
const unsigned nControlPoints = 15;
const unsigned nConstraintsPerCtrlPts = 10;
const double vMax = 85.;

int run_test ()
{
  using namespace boost;
  using namespace boost::assign;

  const double finalPos = 200.;
  CubicBSpline::vector_t params (nControlPoints);

  params[0] = 0;
  params[1] = 0;
  for (unsigned i = 0; i < nControlPoints-4; ++i)
    params[i+2] = finalPos / (nControlPoints - 5) * i;
  params[nControlPoints-2] = finalPos;
  params[nControlPoints-1] = finalPos;

  // Make trajectories.
  CubicBSpline::interval_t timeRange = CubicBSpline::makeInterval (0., 4.);
  CubicBSpline spline (timeRange, 1, params, "before");
  freeTime_t freeTimeTraj (spline, 1.);

  // Define cost.
  Function::matrix_t a (1, freeTimeTraj.parameters ().size ());
  a.setZero ();
  a (0, 0) = -1.;
  Function::vector_t b (1);
  b.setZero ();
  roboptim::NumericLinearFunction cost (a, b);

  // Create problem.
  solver_t::problem_t problem (cost);
  problem.startingPoint () = freeTimeTraj.parameters ();

  // Scale has to remain positive.
  problem.argumentBounds ()[0] = Function::makeLowerInterval (0.);

  const freeTime_t::vector_t freeTimeParams = freeTimeTraj.parameters ();

  spline.freezeCurveStart (problem, 1);
  spline.freezeCurveEnd (problem, 1);

  Function::interval_t vRange = Function::makeUpperInterval (.5 * vMax * vMax);
  LimitSpeed<FreeTimeTrajectory<CubicBSpline> >::addToProblem
    (freeTimeTraj, problem, vRange, nControlPoints * nConstraintsPerCtrlPts);

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
  std::cerr << res << std::endl;

  FreeTimeTrajectory<CubicBSpline> optimizedTrajectory =
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
