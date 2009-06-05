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

#include <roboptim/core/finite-difference-gradient.hh>

#include <roboptim/core/solver-factory.hh>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/spline.hh>
#include <roboptim/trajectory/trajectory-cost.hh>
#include <roboptim/trajectory/spline-length.hh>

#include "common.hh"

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;

typedef boost::variant<const DerivableFunction*,
		       const LinearFunction*> constraint_t;
typedef Solver<DerivableFunction, constraint_t> solver_t;


struct FixStartEnd : public DerivableFunction
{
  FixStartEnd (const vector_t& parameters) throw ()
    : DerivableFunction (parameters.size (), 1),
      parameters_ (parameters)
  {
  }

  ~FixStartEnd () throw () {}

  void
  impl_compute (result_t& res, const argument_t& x) const throw ()
  {
    res (0) = 0.;

    res (0) += (x[0] - parameters_[0]) * (x[0] - parameters_[0]);
    res (0) += (x[1] - parameters_[1]) * (x[1] - parameters_[1]);

    int n = parameters_.size () - 2;
    res (0) += (x[n] - parameters_[n]) * (x[n] - parameters_[n]);

    ++n;
    res (0) += (x[n] - parameters_[n]) * (x[n] - parameters_[n]);
    res (0) -= 1;
  }

  void
  impl_gradient (gradient_t& grad, const argument_t& x, int i) const throw ()
  {
    grad.clear ();

    grad (0) += 2 * (x[0] - parameters_[0]);
    grad (1) += 2 * (x[1] - parameters_[1]);

    int n = parameters_.size () - 2;
    grad (n) += 2 * (x[n] - parameters_[n]);

    ++n;
    grad (n) += 2 * (x[n] - parameters_[n]);
  }

  const vector_t parameters_;
};

int run_test ()
{
  Spline::vector_t params (8);

  // Initial position.
  params[0] = 0.,  params[1] = 0.;
  // Control point 1.
  params[2] = 25.,  params[3] = 100.;
  // Control point 2.
  params[4] = 75.,  params[5] = 0.;
  // Final position.
  params[6] = 100., params[7] = 100.;

  Spline spline (std::make_pair (0., 4.), 2, params);
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
    << plot_xy (spline, interval);

  // Optimize.
  discreteInterval_t costInterval (0., 4., 0.5);
  SplineLength cost (spline, costInterval);

  // Check cost gradient.
  {
    Function::vector_t x (params.size ());
    x.clear ();
    if (! checkGradient (cost, 0, x))
      {
	FiniteDifferenceGradient fdfunction (cost);
	DerivableFunction::gradient_t grad = cost.gradient (x, 0);
	DerivableFunction::gradient_t fdgrad = fdfunction.gradient (x, 0);

	std::cerr << "Cost: " << cost (x) << std::endl
		  << "Grad: " << grad << std::endl
		  << "FD Grad: " << fdgrad << std::endl
		  << "Difference: " << grad - fdgrad << std::endl;
	assert (0 && "Bad gradient");
      }

    x = params;
    if (! checkGradient (cost, 0, x))
      {
	FiniteDifferenceGradient fdfunction (cost);
	DerivableFunction::gradient_t grad = cost.gradient (x, 0);
	DerivableFunction::gradient_t fdgrad = fdfunction.gradient (x, 0);

	std::cerr << "Cost: " << cost (x) << std::endl
		  << "Grad: " << grad << std::endl
		  << "FD Grad: " << fdgrad << std::endl
		  << "Difference: " << grad - fdgrad << std::endl;
	assert (0 && "Bad gradient");
      }
  }

  solver_t::problem_t problem (cost);
  problem.startingPoint () = params;

  FixStartEnd fse (params);

  // Check constraint gradient.
  {
    Function::vector_t x (params.size ());
    x.clear ();
    if (! checkGradient (fse, 0, x))
      {
	FiniteDifferenceGradient fdfunction (fse);
	DerivableFunction::gradient_t grad = fse.gradient (x, 0);
	DerivableFunction::gradient_t fdgrad = fdfunction.gradient (x, 0);

	std::cerr << fse (x) << std::endl
		  << grad << std::endl
		  << fdgrad << std::endl;
	assert (0 && "Bad gradient");
      }

    x = params;
    assert (checkGradient (fse, 0, x));
  }


  problem.addConstraint (&fse, Function::makeBound (0., 0.));

  SolverFactory<solver_t> factory ("cfsqp", problem);
  solver_t& solver = factory ();

  std::cerr << "Cost function (before): " << cost (params) << std::endl;
  std::cerr << "Parameters (before): " << params << std::endl;

  solver_t::result_t res = solver.minimum ();

  std::cerr << "Parameters (after): " << params << std::endl;

  std::cerr << solver << std::endl;

  switch (res.which ())
    {
    case GenericSolver::SOLVER_VALUE:
      {
	Result& result = boost::get<Result> (res);
	spline.setParameters (result.x);
	params = result.x;
	gnuplot << plot_xy (spline, interval);
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
	spline.setParameters (result.x);
	params = result.x;
	std::cerr << result << std::endl;
	gnuplot << plot_xy (spline, interval);
	break;
      }

    case GenericSolver::SOLVER_ERROR:
      {
	SolverError& result = boost::get<SolverError> (res);
	std::cerr << result << std::endl;
      return 1;
      }
    }


  // Check cost gradient (final).
  {
    if (! checkGradient (cost, 0, params))
      {
	FiniteDifferenceGradient fdfunction (cost);
	DerivableFunction::gradient_t grad = cost.gradient (params, 0);
	DerivableFunction::gradient_t fdgrad = fdfunction.gradient (params, 0);

	std::cerr << "Cost: " << cost (params) << std::endl
		  << "Grad: " << grad << std::endl
		  << "FD Grad: " << fdgrad << std::endl
		  << "Difference: " << grad - fdgrad << std::endl;
	assert (0 && "Bad gradient");
      }
  }


  std::cout << (gnuplot << unset ("multiplot"));
  return 0;
}

GENERATE_TEST ()
