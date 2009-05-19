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
#include <roboptim/trajectory/trajectory-cost.hh>

#include "common.hh"

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;

typedef boost::variant<const DerivableFunction*,
		       const LinearFunction*> constraint_t;
typedef Solver<DerivableFunction, constraint_t> solver_t;


struct LengthCost : public TrajectoryCost<Spline>
{
  LengthCost (const Spline& spline, discreteInterval_t interval)
    : TrajectoryCost<Spline> (spline),
      interval_ (interval)
  {
  }

  virtual vector_t operator () (const vector_t& x) const throw ()
  {
    vector_t res (m);

    trajectory_t traj = trajectory_;
    traj.setParameters (x);

    using namespace boost;
    using namespace boost::numeric::ublas;

    for (value_type i = get<0> (interval_); i <= get<1> (interval_);
	 i += get<2> (interval_))
      {
	double tmp = norm_1 (traj.derivative (i, 1));
	res[0] += tmp * tmp;
      }
    return res;
  }

  virtual gradient_t gradient (const vector_t& x, int) const throw ()
  {
    gradient_t grad (n);

    trajectory_t traj = trajectory_;
    traj.setParameters (x);

    using namespace boost;
    using namespace boost::numeric::ublas;

    for (value_type i = get<0> (interval_); i <= get<1> (interval_);
	 i += get<2> (interval_))
      {
	double tmp = norm_1 (traj.variationDerivWrtParam (i, 2));
	grad[0] += tmp * tmp;
      }
    return grad;
  }

    discreteInterval_t interval_;
};

struct FixStartEnd : public DerivableFunction
{
  FixStartEnd (const vector_t& parameters)
    : DerivableFunction (1, 1),
      parameters_ (parameters)
  {
  }

  virtual vector_t operator () (const vector_t& x) const throw ()
  {
    vector_t res (m);
    res (0) = 0.;

    res (0) += (x[0] - parameters_[0]) * (x[0] - parameters_[0]);
    res (0) += (x[1] - parameters_[1]) * (x[1] - parameters_[1]);

    int n = parameters_.size () - 2;
    res (0) += (x[n] - parameters_[n]) * (x[n] - parameters_[n]);

    ++n;
    res (0) += (x[n] - parameters_[n]) * (x[n] - parameters_[n]);

    return res;
  }

  virtual gradient_t gradient (const vector_t& x, int) const throw ()
  {
    gradient_t grad (n);
    grad (0) = 0.;

    grad (0) += (parameters_[0] - x[0]);
    grad (0) += (parameters_[1] - x[1]);

    int n = parameters_.size () - 2;
    grad (0) += (parameters_[n] - x[n]);

    ++n;
    grad (0) += (parameters_[n] - x[n]);

    grad (0) *= 2.;
    return grad;
  }

  const vector_t& parameters_;
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
  gnuplot
    << set ("multiplot")
    << plot_xy (spline, interval);

  // Optimize.
  discreteInterval_t costInterval (0., 5., 0.5);
  LengthCost cost (spline, costInterval);

  solver_t::problem_t problem (cost);
  problem.startingPoint () = params;

  problem.argBounds ()[0] = Function::makeBound (params[0], params[0]);
  problem.argBounds ()[1] = Function::makeBound (params[1], params[1]);

  const int n = params.size ();
  problem.argBounds ()[n - 2] = Function::makeBound (params[n - 2], params[n-2]);
  problem.argBounds ()[n - 1] = Function::makeBound (params[n - 1], params[n-1]);

//   FixStartEnd fse (params);
//   problem.addConstraint (&fse, Function::makeBound (0., 0.));

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

  std::cout << (gnuplot << unset ("multiplot"));
  return 0;
}

GENERATE_TEST ()
