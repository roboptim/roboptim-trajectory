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

#include <boost/assign/list_of.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <roboptim/core/finite-difference-gradient.hh>
#include <roboptim/core/solver-factory.hh>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/spline-length.hh>
#include <roboptim/trajectory/cubic-b-spline.hh>
#include <roboptim/trajectory/trajectory-cost.hh>
#include <roboptim/trajectory/visualization/trajectory.hh>

#include "shared-tests/common.hh"

using namespace roboptim;
using namespace roboptim::trajectory;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;

// Solver type different from the test suite's solver_t
typedef Solver<DifferentiableFunction,
	       boost::mpl::vector<LinearFunction, DifferentiableFunction> >
test_solver_t;

typedef test_solver_t::problem_t::constraints_t constraint_t;

void showSpline (const CubicBSpline& spline);

void showSpline (const CubicBSpline& spline)
{
  std::cout
    << "# Values:" << std::endl
    << "# " << normalize (spline (0.)) << std::endl
    << "# " << normalize (spline (2.5)) << std::endl
    << "# " << normalize (spline (4.)) << std::endl

    << "# 1st derivative:" << std::endl
    << "# " << normalize (spline.derivative (0., 1)) << std::endl
    << "# " << normalize (spline.derivative (2.5, 1)) << std::endl
    << "# " << normalize (spline.derivative (4., 1)) << std::endl

    << "# 2nd derivative:" << std::endl
    << "# " << normalize (spline.derivative (0., 2)) << std::endl
    << "# " << normalize (spline.derivative (2.5, 2)) << std::endl
    << "# " << normalize (spline.derivative (4., 2)) << std::endl

    << "# variationConfigWrtParam:" << std::endl
    << "# " << normalize (spline.variationConfigWrtParam (0.)) << std::endl
    << "# " << normalize (spline.variationConfigWrtParam (2.5)) << std::endl
    << "# " << normalize (spline.variationConfigWrtParam (4.)) << std::endl;
}

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (trajectory_spline_optimization)
{
  using namespace boost::assign;
  CubicBSpline::vector_t params (16);

  // Initial position.
  params[0] = 0.,  params[1] = 0.;
  params[2] = 0.,  params[3] = 0.;
  params[4] = 0.,  params[5] = 0.;
  // Control point 3.
  params[6] = 25.,  params[7] = 100.;
  // Control point 4.
  params[8] = 75.,  params[9] = 0.;
  // Final position.
  params[10] = 100., params[11] = 100.;
  params[12] = 100., params[13] = 100.;
  params[14] = 100., params[15] = 100.;

  CubicBSpline::interval_t timeRange = CubicBSpline::makeInterval (0., 4.);

  CubicBSpline spline (timeRange, 2, params, "before");
  discreteInterval_t interval (0., 4., 0.01);

  showSpline (spline);

  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();
  gnuplot
    << set ("multiplot layout 1,2")
    << set ("grid")
    << plot_xy (spline, interval);

  // Optimize.
  SplineLength cost (spline);

  // Check cost gradient.
  try
    {
      Function::vector_t x (params.size ());
      x.setZero ();
      checkGradientAndThrow (cost, 0, x, 2e-3);

      x = params;
      checkGradientAndThrow (cost, 0, x, 2e-3);
    }
  catch (BadGradient<EigenMatrixDense>& bg)
    {
      std::cerr << bg << std::endl;
      BOOST_CHECK(false);
    }

  test_solver_t::problem_t problem (cost);
  problem.startingPoint () = params;

  spline.freezeCurveStart (problem);
  spline.freezeCurveEnd (problem);

  SolverFactory<test_solver_t> factory (TESTSUITE_SOLVER, problem);
  test_solver_t& solver = factory ();

  // Set optional log file for debugging
  SET_LOG_FILE(solver);

  std::cerr << "Cost function (before): " << cost (params) << std::endl;
  std::cerr << "Parameters (before): " << params << std::endl;

  std::cerr << solver << std::endl;

  test_solver_t::result_t res = solver.minimum ();

  switch (res.which ())
    {
    case GenericSolver::SOLVER_VALUE:
      {
	Result& result = boost::get<Result> (res);
	CubicBSpline optimizedSpline (timeRange, 2, result.x, "after");
	showSpline (optimizedSpline);
	params = result.x;
	gnuplot << plot_xy (optimizedSpline, interval);
	break;
      }

    case GenericSolver::SOLVER_NO_SOLUTION:
      {
	std::cerr << "No solution" << std::endl;
	BOOST_CHECK(false);
      }
    case GenericSolver::SOLVER_VALUE_WARNINGS:
      {
	ResultWithWarnings& result = boost::get<ResultWithWarnings> (res);
	CubicBSpline optimizedSpline (timeRange, 2, result.x, "after");
	showSpline (optimizedSpline);
	params = result.x;
	std::cerr << result << std::endl;
	gnuplot << plot_xy (optimizedSpline, interval);
	break;
      }

    case GenericSolver::SOLVER_ERROR:
      {
	SolverError& result = boost::get<SolverError> (res);
	std::cerr << result << std::endl;
	BOOST_CHECK(false);
      }
    }

  std::cerr << "Parameters (after): " << params << std::endl;

  // Check cost gradient (final).
  try
    {
      checkGradientAndThrow (cost, 0, params, 2e-3);
    }
  catch (BadGradient<EigenMatrixDense>& bg)
    {
      std::cerr << bg << std::endl;
      BOOST_CHECK(false);
    }

  std::cout << (gnuplot << unset ("multiplot"));
}

BOOST_AUTO_TEST_SUITE_END ()
