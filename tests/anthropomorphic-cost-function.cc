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

#include <roboptim/trajectory/sys.hh>

#include <fstream>

#include <boost/assign/list_of.hpp>
#include <boost/mpl/vector.hpp>

#include <roboptim/core/io.hh>
#include <roboptim/core/finite-difference-gradient.hh>
#include <roboptim/core/solver-factory.hh>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>

#include <roboptim/trajectory/anthropomorphic-cost-function.hh>
#include <roboptim/trajectory/free-time-trajectory.hh>
#include <roboptim/trajectory/freeze.hh>
#include <roboptim/trajectory/frontal-speed.hh>
#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/limit-omega.hh>
#include <roboptim/trajectory/orthogonal-speed.hh>
#include <roboptim/trajectory/cubic-b-spline.hh>
#include <roboptim/trajectory/stable-point-state-function.hh>
#include <roboptim/trajectory/trajectory-cost.hh>

#include <roboptim/trajectory/visualization/trajectory.hh>
#include <roboptim/trajectory/visualization/speed.hh>

#include <roboptim/core/plugin/cfsqp.hh>

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;

typedef CFSQPSolver::problem_t::constraints_t constraint_t;
typedef CFSQPSolver solver_t;
typedef FreeTimeTrajectory<Spline> freeTime_t;

// Problem parameters.
const unsigned configurationSpaceSize = 3;
const unsigned nControlPoints = 11;

const unsigned nConstraintsPerCtrlPts = 10;
const double vMax = 0.4;

int optimize (double initialX,
	      double initialY,
	      double initialTheta,
	      double finalX,
	      double finalY,
	      double finalTheta,
	      bool setStartingPoint);

int optimize (double initialX,
	      double initialY,
	      double initialTheta,
	      double finalX,
	      double finalY,
	      double finalTheta,
	      bool setStartingPoint)
{
  using namespace boost;
  using namespace boost::assign;

  CubicBSpline::vector_t params (nControlPoints * configurationSpaceSize);

  const double dx = finalX - initialX;
  const double dy = finalY - initialY;
  const double dtheta = finalTheta - initialTheta;

  for (unsigned i = 0; i < nControlPoints; ++i)
    {
      params[i * configurationSpaceSize] = initialX;
      params[i * configurationSpaceSize + 1] = initialY;
      params[i * configurationSpaceSize + 2] = initialTheta;

      params[i * configurationSpaceSize] +=
	dx / (nControlPoints - 1) * i;
      params[i * configurationSpaceSize + 1] +=
	dy / (nControlPoints - 1) * i;
      params[i * configurationSpaceSize + 2] +=
	dtheta / (nControlPoints - 1) * i;
    }

  // Make trajectories.
  CubicBSpline::interval_t timeRange = CubicBSpline::makeInterval (0., 16.);
  CubicBSpline spline (timeRange, configurationSpaceSize, params, "before");
  freeTime_t freeTimeTraj (spline, 1.);

  // Define cost.
  AnthropomorphicCostFunction<freeTime_t> cost (freeTimeTraj);

  // Create problem.
  solver_t::problem_t problem (cost);
  problem.startingPoint () = freeTimeTraj.parameters ();

  // Scale has to remain positive.
  problem.argumentBounds ()[0] = Function::makeLowerInterval (.5);

  const freeTime_t::vector_t freeTimeParams = freeTimeTraj.parameters ();

  std::vector<Function::size_type> indices;
  indices.push_back (1);
  indices.push_back (2);
  indices.push_back (3);
  indices.push_back (freeTimeParams.size () - 3);
  indices.push_back (freeTimeParams.size () - 2);
  indices.push_back (freeTimeParams.size () - 1);
  makeFreeze (problem) (indices, freeTimeParams);

  // Add constraints on speeds.
  // Frontal
  boost::shared_ptr<DerivableFunction> frontalSpeed (new FrontalSpeed ());
  Function::interval_t vRangeFrontal = Function::makeInterval (0., vMax);
  StablePointStateFunction<freeTime_t>::addToProblem
    (freeTimeTraj, frontalSpeed, 1, problem, vRangeFrontal,
     nControlPoints * nConstraintsPerCtrlPts);

  // // Orthogonal
  // boost::shared_ptr<DerivableFunction> orthogonalSpeed (new OrthogonalSpeed ());
  // Function::interval_t vRangeOrthogonal = Function::makeInterval (-vMax, vMax);
  // StablePointStateFunction<freeTime_t>::addToProblem
  //   (freeTimeTraj, orthogonalSpeed, 1, problem, vRangeOrthogonal,
  //    nControlPoints * nConstraintsPerCtrlPts);

  // Omega (theta dot)
  Function::interval_t vRangeOmega = Function::makeInterval (-.5, .5);
  LimitOmega<freeTime_t>::addToProblem
    (freeTimeTraj, problem, vRangeOmega,
     nControlPoints * nConstraintsPerCtrlPts);

  // Create the solver plug-in.
  SolverFactory<solver_t> factory ("cfsqp", problem);
  solver_t& solver = factory ();

  std::cerr << solver << std::endl;

  std::ofstream trajectoryStream ("trajectory.gp");
  Gnuplot gnuplotTraj = Gnuplot::make_interactive_gnuplot ();
  gnuplotTraj
    << set ("multiplot layout 2, 3 title 'trajectory'")
    << set ("grid")
    << plot_xy (spline)
    << plot_xytheta (spline)
    << plot_speeds (spline);

  solver_t::result_t res = solver.minimum ();
  std::cerr << res << std::endl;

  FreeTimeTrajectory<CubicBSpline> optimizedTrajectory = freeTimeTraj;

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

  CubicBSpline updatedSpline(optimizedTrajectory().timeRange (),
			     optimizedTrajectory().outputSize (),
			     optimizedTrajectory().paramters (),
			     "after");

  std::cerr
    << "Final time range: " << updatedSpline.timeRange () << std::endl
    << "Final length: " << updatedSpline.length () << std::endl;

  gnuplotTraj << plot_xy (updatedSpline)
	      << plot_xytheta (updatedSpline)
	      << plot_speeds (updatedSpline)
	      << unset ("multiplot");
  trajectoryStream << gnuplotTraj;

  return 0;
}
