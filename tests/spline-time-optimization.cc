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

#include "shared-tests/fixture.hh"

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

using namespace roboptim;
using namespace roboptim::trajectory;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;
using namespace roboptim::trajectory::visualization::gnuplot;

// Solver type different from the test suite's solver_t
typedef Solver<DifferentiableFunction,
	       boost::mpl::vector<LinearFunction, DifferentiableFunction> >
test_solver_t;

typedef CubicBSpline spline_t;
typedef test_solver_t::problem_t::constraints_t constraint_t;
typedef FreeTimeTrajectory<spline_t> freeTime_t;
typedef spline_t::size_type size_type;

typedef Eigen::MatrixXd matrix_t;

// Limited data => use full precision serialization
typedef matrix_t store_t;

// Problem parameters.
const unsigned nControlPoints = 6;
const unsigned nConstraintsPerCtrlPts = 5;
const double vMax = 85.;

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (trajectory_spline_time_optimization)
{
  using namespace boost;
  using namespace boost::assign;

  // Fix clash with std::set
  using roboptim::visualization::gnuplot::set;

  boost::shared_ptr<boost::test_tools::output_test_stream>
    output = retrievePattern ("spline-time-optimization");

  size_type row_iter = 0;

  matrix_t reference = readMatrix<store_t> ("spline-time-optimization.dat");
  if (reference.size () == 0) {
    throw std::runtime_error
      ("failed to open file \"spline-time-optimization.dat\"");
  }

  double tol = 5e-5;

  const double finalPos = 200.;
  spline_t::vector_t params (nControlPoints);

  params[0] = 0;
  params[1] = 0;
  for (unsigned i = 0; i < nControlPoints-4; ++i)
    params[i+2] = finalPos / (nControlPoints - 5) * i;
  params[nControlPoints-2] = finalPos;
  params[nControlPoints-1] = finalPos;

  // Make trajectories.
  spline_t::interval_t timeRange = spline_t::makeInterval (0., 4.);
  spline_t spline (timeRange, 1, params, "before");
  freeTime_t freeTimeTraj (spline, 1.);

  // Define cost.
  Function::matrix_t a (1, freeTimeTraj.parameters ().size ());
  a.setZero ();
  a (0, 0) = -1.;
  Function::vector_t b (1);
  b.setZero ();
  roboptim::NumericLinearFunction cost (a, b);

  // Create problem.
  test_solver_t::problem_t problem (cost);
  problem.startingPoint () = freeTimeTraj.parameters ();

  // Scale has to remain positive.
  problem.argumentBounds ()[0] = Function::makeLowerInterval (0.);

  const freeTime_t::vector_t freeTimeParams = freeTimeTraj.parameters ();

  spline.freezeCurveStart (problem, 1);
  spline.freezeCurveEnd (problem, 1);

  Function::interval_t vRange = Function::makeUpperInterval (.5 * vMax * vMax);
  LimitSpeed<FreeTimeTrajectory<spline_t> >::addToProblem
    (freeTimeTraj, problem, vRange, nControlPoints * nConstraintsPerCtrlPts);

  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();

  gnuplot
    << set ("multiplot layout 1,2 title "
	    "'variation of speed before and after optimization'")
    << set ("grid");
  gnuplot << plot_limitSpeed (freeTimeTraj, vMax);

  SolverFactory<test_solver_t> factory (TESTSUITE_SOLVER, problem);
  test_solver_t& solver = factory ();

  // Ipopt-specific parameters
  // WARNING: these parameters may not be relevant! These are only set to
  // prevent hour-long unit testing...
  solver.parameters()["ipopt.linear_solver"].value = IPOPT_LINEAR_SOLVER;
  solver.parameters()["ipopt.tol"].value = 1e-3;
  solver.parameters()["ipopt.acceptable_tol"].value = 5e-2;
  solver.parameters()["ipopt.mu_strategy"].value = "adaptive";

  (*output) << solver;

  test_solver_t::result_t res = solver.minimum ();
  std::cerr << res << std::endl;

  FreeTimeTrajectory<spline_t> optimizedTrajectory =
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
      BOOST_CHECK(false);
    }

  gnuplot << plot_limitSpeed (optimizedTrajectory, vMax);
  // Do not embed Gnuplot data for comparison because of numerical
  // differences
  std::cout << (gnuplot << unset ("multiplot"));

  std::cout << output->str () << std::endl;

  BOOST_CHECK (output->match_pattern ());

  // Loop over the interval of definition
  double eval_step = 0.01;
  size_type n_eval = static_cast<size_type> (std::floor ((timeRange.second - timeRange.first)/eval_step)+1);
  matrix_t eval_data (n_eval, 2);
  eval_data.setZero ();

  for (double t = timeRange.first; t < timeRange.second;
       t += eval_step)
    {
      double res = optimizedTrajectory (t)[0];
      BOOST_CHECK_SMALL (t - reference (row_iter, 0), tol);
      BOOST_CHECK_SMALL (res - reference (row_iter, 1), tol);
      eval_data.row (row_iter) << t, res;
      ++row_iter;
    }

  // Generate the archive
  //writeMatrix<matrix_t> ("spline-time-optimization.dat", eval_data);

  // Compare data
  BOOST_CHECK (allclose (eval_data, reference, tol));
}

BOOST_AUTO_TEST_SUITE_END ()
