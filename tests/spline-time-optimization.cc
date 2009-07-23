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

int run_test ()
{
  using namespace boost;
  using namespace boost::assign;
  Spline::vector_t params (4);

  // Initial position.
  params[0] = 0.;
  // Control point 1.
  params[1] = 25.;
  // Control point 2.
  params[2] = 75.;
  // Final position.
  params[3] = 100.;

  // Make trajectories.
  Spline::interval_t timeRange = Spline::makeInterval (0., 4.);
  Spline spline (timeRange, 1, params, "before");
  FreeTimeTrajectory<Spline::derivabilityOrder> freeTimeTraj (spline, 1.);

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
    a(0, 0) = 1.;
    shared_ptr<NumericLinearFunction> speedPositivity
      (new NumericLinearFunction (a, b));
    problem.addConstraint
      (static_pointer_cast<DerivableFunction> (speedPositivity),
       Function::makeLowerInterval (0.));
  }


  typedef Freeze<DerivableFunction, constraint_t, LinearFunction> freeze_t;
  freeze_t freeze (problem,
		   list_of <freeze_t::frozenArgument_t>
		   (1, params[0])
		   (params.size (), params[params.size () - 1]));
  freeze ();


  Function::interval_t vRange (0., 100.);
  LimitSpeed<FreeTimeTrajectory<Spline::derivabilityOrder> >::addToProblem (freeTimeTraj, problem, vRange, 10);

  SolverFactory<solver_t> factory ("cfsqp", problem);
  solver_t& solver = factory ();


  std::cout << "Cost function (before): " << cost (freeTimeTraj.parameters ())
	    << std::endl
	    << "Parameters (before): " << freeTimeTraj.parameters ()
	    << std::endl;

  std::cout << solver << std::endl;

  solver_t::result_t res = solver.minimum ();

  switch (res.which ())
    {
    case GenericSolver::SOLVER_VALUE:
      {
	Result& result = boost::get<Result> (res);
	FreeTimeTrajectory<Spline::derivabilityOrder> optimizedTrajectory =
	  freeTimeTraj;
	optimizedTrajectory.setParameters (result.x);
	std::cout << "Parameters (after): " << optimizedTrajectory.parameters ()
		  << std::endl;
	std::ofstream limitSpeedStream ("limit-speed.gp");
	limitSpeedStream
	  << (Gnuplot::make_interactive_gnuplot ()
	      << plot_limitSpeed (optimizedTrajectory));
	break;
      }

    case GenericSolver::SOLVER_NO_SOLUTION:
      {
	std::cout << "No solution" << std::endl;
	return 1;
      }
    case GenericSolver::SOLVER_VALUE_WARNINGS:
      {
	ResultWithWarnings& result = boost::get<ResultWithWarnings> (res);
	FreeTimeTrajectory<Spline::derivabilityOrder> optimizedTrajectory =
	  freeTimeTraj;
	optimizedTrajectory.setParameters (result.x);
	std::cout << result << std::endl;
	std::cout << "Parameters (after): " << optimizedTrajectory.parameters ()
		  << std::endl;
	std::ofstream limitSpeedStream ("limit-speed.gp");
	limitSpeedStream
	  << (Gnuplot::make_interactive_gnuplot ()
	      << plot_limitSpeed (optimizedTrajectory));
	break;
      }

    case GenericSolver::SOLVER_ERROR:
      {
	SolverError& result = boost::get<SolverError> (res);
	std::cout << result << std::endl;
      return 1;
      }
    }

  return 0;
}

GENERATE_TEST ()
