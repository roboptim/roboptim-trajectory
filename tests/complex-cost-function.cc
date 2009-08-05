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
#include <roboptim/trajectory/orthogonal-speed.hh>
#include <roboptim/trajectory/spline.hh>
#include <roboptim/trajectory/trajectory-cost.hh>

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

const double initialX = 0.;
const double initialY = 0.;
const double initialTheta = M_PI / 2.;

const double finalX = 5.;
const double finalY = 0.;
const double finalTheta = M_PI / 2.;

namespace roboptim
{
  void normalizeAngle (Function::vector_t& parameters,
		       unsigned configurationSize,
		       unsigned thetaIdx,
		       bool isFreeTime);

  void normalizeAngle (Function::vector_t& parameters,
		       unsigned configurationSize,
		       unsigned thetaIdx,
		       bool isFreeTime)
  {
    double thetaPrev = 0.;
    double offset = isFreeTime ? 1. : 0.;
    for (unsigned i = offset; i < parameters.size () / configurationSize; ++i)
      {
	double& theta = parameters[offset + i * configurationSize + thetaIdx];
	if (theta - thetaPrev > M_PI)
	  theta -= M_PI * 2;
	else if (theta - thetaPrev < -M_PI)
	  theta += M_PI * 2;
	thetaPrev = theta;
      }
  }

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
      vector_t p = trajectory.parameters ();
      //normalizeAngle (p, 3, 2, true);

      //FIXME: +1 as we're doing free time trajectory.
      //should be generic.
      const value_type dx = p[p.size () - 3] - p[1 + 0];
      const value_type dy = p[p.size () - 2] - p[1 + 1];
      const value_type dtheta = p[p.size () - 1] - p[1 + 2];

      alpha3_ = this->alpha3 (std::fabs (dtheta), dx * dx + dy * dy);
    }

    ~CostFunction () throw ()
    {}

  protected:
    struct ComputeIntegral
    {
      ComputeIntegral (const T& traj,
		       const vector_t& alpha,
		       const double& alpha3,
		       double& res)
	: traj_ (traj),
	  alpha_ (alpha),
	  alpha3_ (alpha3),
	  res_ (res)
      {}

      void operator () (const double& t)
      {
	FrontalSpeed<T> frontalSpeed (traj_);
	OrthogonalSpeed<T> orthogonalSpeed (traj_);

	vector_t t_ (1);
	t_[0] = t;
	const value_type u1 = frontalSpeed.gradient (t_)[0];
	const value_type u2 = traj_.derivative (t, 2)[2];
	const value_type u3 = orthogonalSpeed.gradient (t_)[0];
	res_ +=
	  alpha_[0]
	  + alpha_[1] * u1 * u1
	  + alpha_[2] * u2 * u2
	  + alpha3_ * u3 * u3;
      }

    private:
      const T& traj_;
      const vector_t& alpha_;
      const double& alpha3_;
      double& res_;
    };

    void impl_compute (result_t& res, const argument_t& p) const throw ()
    {
      res.clear ();

      vector_t params = p;
      //normalizeAngle (params, 3, 2, true);

      boost::scoped_ptr<T> updatedTrajectory (trajectory_.clone ());
      updatedTrajectory->setParameters (params);

      discreteInterval_t interval =
	makeDiscreteInterval (updatedTrajectory->timeRange (), .1);

      ComputeIntegral ci (*updatedTrajectory, alpha_, alpha3_, res[0]);
      foreach (interval, ci);
    }

    void impl_gradient (gradient_t& grad, const argument_t& p, size_type i)
      const throw ()
    {
      FiniteDifferenceGradient fdfunction (*this);
      fdfunction.gradient (grad, p, 0);
    }

  private:

    value_type alpha3 (value_type deltaTheta, value_type dsquare) const throw ()
    {
      const value_type ksi1 = M_PI / 18.;
      const value_type ksi2 = .5;

      return alpha_[3] * (1. + (deltaTheta / ksi1)) * (1 + (dsquare / ksi2));
    }

    const T& trajectory_;
    const vector_t alpha_;
    value_type alpha3_;
  };
} // end of namespace roboptim.


template <typename T>
Command displayComponents (const T& traj,
			   typename T::value_type step = .01);

template <typename T>
Command displaySpeeds (const T& traj,
		       typename T::value_type step = .01);

template <typename T>
Command displayComponents (const T& traj,
			   typename T::value_type step)
{
  using boost::format;
  assert (traj.outputSize () == 3);
  Function::value_type min = Function::getLowerBound (traj.timeRange ());
  Function::value_type max = Function::getUpperBound (traj.timeRange ());
  Function::discreteInterval_t interval (min, max, step);

  if (min + step > max)
    throw std::string ("bad interval");

  std::string str = (boost::format ("plot '-' title '%1% (0)' with lines")
		     % traj.getName ()).str ();
  for (unsigned i = 1; i < traj.outputSize (); ++i)
    {
      str += (format (", '-' title '%1% (%2%)' with lines")
	      % traj.getName ()
	      % i).str ();
    }
  str += "\n";
  for (unsigned component = 0; component < traj.outputSize (); ++component)
    {
      for (double i = step; i < 1. - step; i += step)
	{
	  StableTimePoint timePoint = i * tMax;
	  Function::vector_t res = traj (timePoint);
	  str += (format ("%1f %2f\n")
		  % timePoint.getTime (traj.timeRange ())
		  % res [component]).str ();
	}
      str += "e\n";
    }
  return Command (str);
}

template <typename T>
Command displaySpeeds (const T& traj,
		       typename T::value_type step)
{
  using boost::format;
  assert (traj.outputSize () == 3);
  Function::value_type min = Function::getLowerBound (traj.timeRange ());
  Function::value_type max = Function::getUpperBound (traj.timeRange ());
  Function::discreteInterval_t interval (min, max, step);

  if (min + step > max)
    throw std::string ("bad interval");

  std::string str =
    (boost::format ("plot '-' title '%1% (frontal)' with lines, "
		    "'-' title '%1% (orthogonal)' with lines\n")
     % traj.getName ()).str ();

  {
    FrontalSpeed<T> frontalSpeed (traj);
    for (double i = step; i < 1. - step; i += step)
      {
	StableTimePoint timePoint = i * tMax;
	Function::vector_t t_ (1);
	t_[0] = timePoint.getTime (traj.timeRange ());

	str += (format ("%1f %2f\n")
		% t_[0]
		% frontalSpeed (t_)[0]).str ();
      }
    str += "e\n";
  }

  {
    OrthogonalSpeed<T> orthogonalSpeed (traj);
    for (double i = step; i < 1. - step; i += step)
      {
	StableTimePoint timePoint = i * tMax;
	Function::vector_t t_ (1);
	t_[0] = timePoint.getTime (traj.timeRange ());

	str += (format ("%1f %2f\n")
		% t_[0]
		% orthogonalSpeed (t_)[0]).str ();
      }
    str += "e\n";
  }
  return Command (str);
}


int run_test ()
{
  using namespace boost;
  using namespace boost::assign;

  Spline::vector_t params (nControlPoints * configurationSpaceSize);

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
  Spline::interval_t timeRange = Spline::makeInterval (0., 16.);
  Spline spline (timeRange, configurationSpaceSize, params, "before");
  freeTime_t freeTimeTraj (spline, 1.);

  // Define cost.
  Function::vector_t alpha (4);
  alpha[0] = 1.;
  alpha[1] = 10.;
  alpha[2] = 10.;
  alpha[3] = 5.;
  CostFunction<freeTime_t> cost (freeTimeTraj, alpha);

  // Create problem.
  solver_t::problem_t problem (cost);
  //problem.startingPoint () = freeTimeTraj.parameters ();

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

  SolverFactory<solver_t> factory ("cfsqp", problem);
  solver_t& solver = factory ();

  std::cerr << solver << std::endl;

  std::ofstream trajectoryStream ("trajectory.gp");
  Gnuplot gnuplotTraj = Gnuplot::make_interactive_gnuplot ();
  gnuplotTraj
    << set ("multiplot layout 2, 3 title 'trajectory'")
    << set ("grid")
    << plot_xy (spline)
    << displayComponents (spline)
    << displaySpeeds (spline);

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

  std::cerr
    << "Final time range: " << updatedSpline.timeRange () << std::endl
    << "Final length: " << updatedSpline.length () << std::endl;

  gnuplotTraj << plot_xy (updatedSpline)
	      << displayComponents (updatedSpline)
	      << displaySpeeds (updatedSpline)
	      << unset ("multiplot");
  trajectoryStream << gnuplotTraj;

  return 0;
}

GENERATE_TEST ()

