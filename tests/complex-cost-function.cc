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

#include <roboptim/core/io.hh>
#include <roboptim/core/finite-difference-gradient.hh>
#include <roboptim/core/solver-factory.hh>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>

#include <roboptim/trajectory/free-time-trajectory.hh>
#include <roboptim/trajectory/freeze.hh>
#include <roboptim/trajectory/frontal-speed.hh>
#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/limit-speed.hh>
#include <roboptim/trajectory/orthogonal-speed.hh>
#include <roboptim/trajectory/spline-length.hh>
#include <roboptim/trajectory/spline.hh>
#include <roboptim/trajectory/trajectory-cost.hh>

#include <roboptim/trajectory/visualization/limit-speed.hh>
#include <roboptim/trajectory/visualization/trajectory.hh>

#include <roboptim/core/plugin/cfsqp.hh>

#include "common.hh"

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;

typedef CFSQPSolver::problem_t::constraints_t constraint_t;
typedef CFSQPSolver solver_t;
typedef FreeTimeTrajectory<Spline::derivabilityOrder> freeTime_t;

// Problem parameters.
const unsigned configurationSpaceSize = 3;
const unsigned nControlPoints = 11;
const unsigned nDiscretizationPoints = 10;
const unsigned nConstraintsPerCtrlPts = 0;

const double finalX = 200.;
const double finalY = 200.;
const double finalTheta = 3.14;

const double vMax = 200.;

namespace roboptim
{
  // Cost function from ``An optimal control model unifying holonomic
  // and nonholonomic walking'' Katja Mombaur, Jean-Paul Laumond,
  // Eiichi Yoshida
  // (2008 8th IEEE-RAS Interational Conference on Humanoid Robots).
  template <typename T>
  class CostFunction : public DerivableFunction
  {
  public:
    CostFunction (const T& trajectory, const vector_t alpha) throw ()
      : DerivableFunction (trajectory.parameters ().size (), 1,
			   "cost function"),
	trajectory_ (trajectory),
	alpha_ (alpha),
	alpha3_ ()
    {
      //FIXME: +1 as we're doing free time trajectory.
      //should be generic.
      const vector_t& p = trajectory.parameters ();
      const value_type dx = p[1 + 0] - p[p.size () - 3];
      const value_type dy = p[1 + 1] - p[p.size () - 2];
      const value_type dtheta = p[p.size () - 1] - p[1 + 2];

      alpha3_ = this->alpha3 (dtheta, std::sqrt (dx * dx + dy * dy));
    }

    ~CostFunction () throw ()
    {}

  protected:
    void impl_compute (result_t& res, const argument_t& p) const throw ()
    {
      res.clear ();

      boost::scoped_ptr<T> updatedTrajectory (trajectory_.clone ());
      updatedTrajectory->setParameters (p);

      const value_type delta = 1. / nDiscretizationPoints;

      for (double i = 0.; i < 1.; i += delta)
	{
	  FrontalSpeed<T> frontalSpeed (i * tMax, *updatedTrajectory);
	  OrthogonalSpeed<T> orthogonalSpeed (i * tMax, *updatedTrajectory);

	  const value_type u1 = frontalSpeed.gradient (p)[0];
	  const value_type u2 = updatedTrajectory->derivative (i * tMax, 2)[2];
	  const value_type u3 = orthogonalSpeed.gradient (p)[0];
	  res[0] +=
	    alpha_[0]
	    + alpha_[1] * u1 * u1
	    + alpha_[2] * u2 * u2
	    + alpha3_ * u3 * u3;
	}
    }

    void impl_gradient (gradient_t& grad, const argument_t& p, size_type i)
      const throw ()
    {
      FiniteDifferenceGradient fdfunction (*this);
      fdfunction.gradient (grad, p, 0);
    }

  private:

    value_type alpha3 (value_type deltaTheta, value_type d) const throw ()
    {
      const value_type ksi1 = 0.174532925;
      const value_type ksi2 = 0.5;

      return alpha_[3] * (1. + (deltaTheta / ksi1)) * (1 + ((d * d) / ksi2));
    }

    const T& trajectory_;
    const vector_t alpha_;
    value_type alpha3_;
  };
} // end of namespace roboptim.

int run_test ()
{
  using namespace boost;
  using namespace boost::assign;

  Spline::vector_t params (nControlPoints * configurationSpaceSize);

  for (unsigned i = 0; i < nControlPoints; ++i)
    {
      params[i * configurationSpaceSize] = finalX / (nControlPoints - 1) * i;
      params[i * configurationSpaceSize + 1] = finalY / (nControlPoints - 1) * i;
      params[i * configurationSpaceSize + 2] = finalTheta / (nControlPoints - 1) * i;
    }

  // Make trajectories.
  Spline::interval_t timeRange = Spline::makeInterval (0., 4.);
  Spline spline (timeRange, configurationSpaceSize, params, "before");
  freeTime_t freeTimeTraj (spline, 1.);

  // Define cost.
  Function::vector_t alpha (4);
  alpha[0] = 1.;
  alpha[1] = 1.;
  alpha[2] = .5;
  alpha[3] = 2.;
  CostFunction<freeTime_t> cost (freeTimeTraj, alpha);

  // Create problem.
  solver_t::problem_t problem (cost);
  problem.startingPoint () = freeTimeTraj.parameters ();

  // Scale has to remain positive.
  problem.argumentBounds ()[0] = Function::makeLowerInterval (0.);

  const freeTime_t::vector_t freeTimeParams = freeTimeTraj.parameters ();

  std::vector<Function::size_type> indices;
  indices.push_back (1);
  indices.push_back (2);
  indices.push_back (3);
  indices.push_back (freeTimeParams.size () - 3);
  indices.push_back (freeTimeParams.size () - 2);
  indices.push_back (freeTimeParams.size () - 1);
  makeFreeze (problem) (indices, freeTimeParams);

  Function::interval_t vRange = Function::makeUpperInterval (.5 * vMax * vMax);
  LimitSpeed<FreeTimeTrajectory<Spline::derivabilityOrder> >::addToProblem
    (freeTimeTraj, problem, vRange, nControlPoints * nConstraintsPerCtrlPts);

  std::ofstream trajectoryStream ("trajectory.gp");
  std::ofstream limitSpeedStream ("limit-speed.gp");
  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();

  gnuplot
    << set ("multiplot layout 1,2 title "
	    "'variation of speed before and after optimization'")
    << set ("grid");
  gnuplot << plot_limitSpeed (freeTimeTraj, vMax);

  SolverFactory<solver_t> factory ("cfsqp", problem);
  solver_t& solver = factory ();

  std::cerr << solver << std::endl;

  Gnuplot gnuplotTraj = Gnuplot::make_interactive_gnuplot ();
  gnuplotTraj
    << set ("multiplot layout 1, 2 title 'trajectory'")
    << set ("grid")
    << plot_xy (spline);

  solver_t::result_t res = solver.minimum ();
  std::cerr << res << std::endl;

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

  Spline updatedSpline
    (optimizedTrajectory.timeRange (),
     configurationSpaceSize,
     removeScaleFromParameters (optimizedTrajectory.parameters ()),
     "after");

  gnuplotTraj << plot_xy (updatedSpline);
  trajectoryStream << (gnuplotTraj << unset ("multiplot"));

  gnuplot << plot_limitSpeed (optimizedTrajectory, vMax);
  limitSpeedStream << (gnuplot << unset ("multiplot"));

  return 0;
}

GENERATE_TEST ()

