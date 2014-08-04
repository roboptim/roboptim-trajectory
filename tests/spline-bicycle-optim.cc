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

#include <roboptim/trajectory/freeze.hh>
#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/spline-length.hh>
#include <roboptim/trajectory/cubic-b-spline.hh>
#include <roboptim/trajectory/trajectory-sum-cost.hh>
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

typedef TrajectorySumCost<CubicBSpline>::discreteStableTimePointInterval_t
discreteStableTimePointInterval_t;

/*
  Parameter of the cost function
*/
//#define BICYCLE_COST_FUNCTION
static double m=.5;

class PositiveCostVar : public DerivableFunction
{
public:
  PositiveCostVar() : DerivableFunction(4, 1, "positive variation of height")
  {
  }

  void impl_compute(result_t& res, const argument_t& x) const
  {
#ifdef BICYCLE_COST_FUNCTION
    res[0] = std::max(0., m*x[3]);
#else
    double y0 = x[1];
    double x1 = x[2];
    double y1 = x[3];
    double v = 1+m*y0;
    res[0] = .5*v*(x1*x1+y1*y1);
#endif
  }

  void impl_gradient(gradient_t& grad, const argument_t& x,
		     size_type functionId = 0) const
  {
    assert (functionId == 0);
    grad.setZero ();
#ifdef BICYCLE_COST_FUNCTION
    if (x[3] > 0) {
      grad[3] =  m;
    } else if (x[3] == 0) {
      grad[3] = .5*m;
    } else {
      grad[3] = 0;
    }
#else
    double y0 = x[1];
    double x1 = x[2];
    double y1 = x[3];
    double v = 1+m*y0;

    grad[1] = .5*m*(x1*x1+y1*y1);
    grad[2] = v*x1;
    grad[3] = v*y1;
#endif
  }
};


BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (trajectory_spline_bibycle_optim)
{
  using namespace boost::assign;
  CubicBSpline::vector_t params (22);

  // Initial position.
  params[0] = 0.,  params[1] = 0.;
  // Control point .
  params[2] = .1, params[3] = -2.0;
  // Control point .
  params[4] = .2, params[5] = -2.0;
  // Control point .
  params[6] = .3, params[7] = -2.0;
  // Control point .
  params[8] = .4, params[9] = -2.0;
  // Control point .
  params[10] = .5, params[11] = -2.0;
  // Control point .
  params[12] = .6, params[13] = -2.0;
  // Control point .
  params[14] = .7, params[15] = -2.0;
  // Control point .
  params[16] = .8, params[17] = -2.0;
  // Control point .
  params[18] = .9, params[19] = -2.0;
  // Final position.
  params[20] = 1.,  params[21] = .1;

  CubicBSpline::interval_t timeRange = CubicBSpline::makeInterval (0., 1.);
  std::vector<roboptim::Function::interval_t> defDomain;

  for (unsigned int i=0; i<11; i++) {
    defDomain.push_back(Function::makeInterval(0.0, 1.0));
    defDomain.push_back(Function::makeInterval(-2.0, 10.0));
  }

  CubicBSpline spline (timeRange, 2, params, "before");
  discreteStableTimePointInterval_t interval
    (0. * tMax, 1. * tMax, 0.01 * tMax);

  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();
  gnuplot
    << set ("multiplot layout 1,2")
    << set ("grid")
    << plot_xy (spline);

  // Optimize.
  boost::shared_ptr<PositiveCostVar>
    positiveCostVarShPtr(new PositiveCostVar()) ;

  TrajectorySumCost < CubicBSpline > sumCost (spline, positiveCostVarShPtr,
					      interval, 1);
  // Check cost gradient.
  try
    {
      Function::vector_t x (params.size ());
      x.setZero ();
      checkGradientAndThrow (sumCost, 0, x, 2e-3);

      x = params;
      checkGradientAndThrow (sumCost, 0, x, 2e-3);
    }
  catch (BadGradient<EigenMatrixDense>& bg)
    {
      std::cout << bg << std::endl;
      BOOST_CHECK(false);
    }

  test_solver_t::problem_t problem (sumCost);
  problem.startingPoint () = params;
  problem.argumentBounds() = defDomain;

  std::vector<Function::size_type> indices;
  indices.push_back (0);
  indices.push_back (1);
  indices.push_back (params.size () - 2);
  indices.push_back (params.size () - 1);
  makeFreeze (problem) (indices, params);

  SolverFactory<test_solver_t> factory (TESTSUITE_SOLVER, problem);
  test_solver_t& solver = factory ();

  // Set optional log file for debugging
  SET_LOG_FILE(solver);

  std::cerr << "Cost function (before): " << sumCost (params) << std::endl;
  std::cerr << "Parameters (before): " << params << std::endl;

  std::cerr << solver << std::endl;

  test_solver_t::result_t res = solver.minimum ();

  switch (res.which ())
    {
    case GenericSolver::SOLVER_VALUE:
      {
	Result& result = boost::get<Result> (res);
	CubicBSpline optimizedSpline (timeRange, 2, result.x, "after");
	params = result.x;
	gnuplot << plot_xy (optimizedSpline);
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
	params = result.x;
	std::cerr << result << std::endl;
	gnuplot << plot_xy (optimizedSpline);
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
      checkGradientAndThrow (sumCost, 0, params, 2e-3);
    }
  catch (BadGradient<EigenMatrixDense>& bg)
    {
      std::cout << bg << std::endl;

      BOOST_CHECK(false);
    }

  std::cout << (gnuplot << unset ("multiplot"));
}

BOOST_AUTO_TEST_SUITE_END ()
