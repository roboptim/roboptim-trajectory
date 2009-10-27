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

#include <boost/format.hpp>
#include <boost/scoped_ptr.hpp>

#include <roboptim/core/finite-difference-gradient.hh>
#include <roboptim/core/io.hh>

#include <roboptim/trajectory/free-time-trajectory.hh>
#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/spline.hh>
#include <roboptim/trajectory/trajectory-cost.hh>

#include "common.hh"

using boost::format;
using boost::io::group;

using namespace roboptim;


typedef FreeTimeTrajectory<Spline::derivabilityOrder> freeTime_t;

struct ConfigWrtParam : public DerivableFunction
{
  ConfigWrtParam (const freeTime_t& traj, int order, double t) throw ()
    : DerivableFunction (traj.parameters ().size (),
			 (order < 0)
			 ? traj.outputSize () : traj.parameters ().size (),
			 "config wrt param"),
      traj_ (traj),
      order_ (order),
      t_ (t)
  {
  }

  ~ConfigWrtParam () throw ()
  {
  }

  void compute (result_t& res, const argument_t& p, int order) const throw ()
  {
    if (order < 0)
      {
	boost::scoped_ptr<freeTime_t> updatedTrajectory (traj_.clone ());
	updatedTrajectory->setParameters (p);
	res = (*updatedTrajectory) (t_);
      }
    else
      {
	boost::scoped_ptr<freeTime_t> updatedTrajectory (traj_.clone ());
	updatedTrajectory->setParameters (p);
	matrix_t tmp = updatedTrajectory->variationDerivWrtParam (t_, order);
	res = row (tmp, 0);
      }
  }

  void
  impl_compute (result_t& res, const argument_t& p) const throw ()
  {
    compute (res, p, order_);
  }

  void
  impl_gradient (gradient_t& grad, const argument_t& p, size_type i)
    const throw ()
  {
    compute (grad, p, order_ + 1);
  }

  const freeTime_t& traj_;
  const int order_;
  double t_;
};


void printTable (const Spline& spline, const freeTime_t& freeTimeTraj);

void printTable (const Spline& spline, const freeTime_t& freeTimeTraj)
{
  double tmin = Function::getLowerBound (spline.timeRange ());
  double tmax = Function::getUpperBound (spline.timeRange ());
  double fttTmin = Function::getLowerBound (freeTimeTraj.timeRange ());
  double fttTmax = Function::getUpperBound (freeTimeTraj.timeRange ());

  std::cout << format ("Spline range: [%1%, %2%]") % tmin % tmax << std::endl
	    << format ("FTT range: [%1%, %2%]") % fttTmin % fttTmax << std::endl;

  format fmter ("| %1% %|6t||| %2% %|20t|| %3% %|35t||| %4% %|55t|| %5% %|79t||");
  std::cout << "/---------------------------------------"
	    << "---------------------------------------\\" << std::endl
	    << fmter
    % "T" % "Value" % "Value (ftt)"
    % "Derivative" % "Derivative (ftt)"
	    << std::endl
	    << "----------------------------------------"
	    << "----------------------------------------"
	    << std::endl;

  for (double t = fttTmin; t <= fttTmax + 1e-3; t += .1)
    {
      if (t > fttTmax)
	t = fttTmax;

      fmter % t;
      if (tmin <= t && t <= tmax)
	fmter % spline (t)[0];
      else
	fmter % "N/A";
      fmter % freeTimeTraj (t)[0];
      if (tmin <= t && t <= tmax)
	fmter % spline.derivative (t, 1)[0];
      else
	fmter % "N/A";
      fmter % freeTimeTraj.derivative (t, 1)[0];

      std:: cout << fmter << std::endl;
    }
  std::cout << "\\---------------------------------------"
	    << "---------------------------------------/"
	    << std::endl << std::endl;

  std::cout << "Gradient check." << std::endl;

  for (double t = fttTmin; t <= fttTmax + 1e-3; t += .1)
    {
      if (t > fttTmax)
	t = fttTmax;

      Spline::vector_t x (1);
      x[0] = t;

      try
	{
	  std::cout << "Spline gradient." << std::endl;
	  if (tmin <= t && t <= tmax)
	    checkGradientAndThrow (spline, 0, x);
	  std::cout << "Free time trajectory gradient." << std::endl;
	  checkGradientAndThrow (freeTimeTraj, 0, x);
	}
      catch (BadGradient& bg)
	{
	  std::cout << bg << std::endl;
	}
    }
  std::cout << std::endl << std::endl;


  std::cout << "Variation of the configuration w.r.t to parameters:" << std::endl;
  format fmterConfig ("%1% %|50t|%2%");
  for (double t = fttTmin; t <= fttTmax + 1e-3; t += .1)
    {
      if (t > fttTmax)
	t = fttTmax;

      if (tmin <= t && t <= tmax)
	{
	  Spline::jacobian_t splineVarConfig =
	    spline.variationConfigWrtParam (t);
	  fmterConfig % splineVarConfig;
	}
      else
	fmterConfig % "N/A";

      freeTime_t::jacobian_t fttVarConfig =
	freeTimeTraj.variationConfigWrtParam (t);
      fmterConfig % fttVarConfig;

      try
	{
	  ConfigWrtParam configWrtParam (freeTimeTraj, -1, t);
	  checkGradientAndThrow (configWrtParam, 0,
				 freeTimeTraj.parameters ());
	}
      catch (BadGradient& bg)
	{
	  std::cout << bg << std::endl;
	}

      std::cout << fmterConfig << std::endl;
    }
  std::cout << std::endl << std::endl;

  std::cout << "Variation of the derivative w.r.t to parameters:" << std::endl;
  format fmterDeriv ("%1% %|50t|%2%");
  for (double t = fttTmin; t <= fttTmax + 1e-3; t += .1)
    {
      if (t > fttTmax)
	t = fttTmax;

      if (tmin <= t && t <= tmax)
	{
	  Spline::jacobian_t splineVarDeriv =
	    spline.variationDerivWrtParam (t, 1);
	  fmterDeriv % splineVarDeriv;
	}
      else
	fmterDeriv % "N/A";

      freeTime_t::jacobian_t fttVarDeriv =
	freeTimeTraj.variationDerivWrtParam (t, 1);
      fmterDeriv % fttVarDeriv;

      try
	{
	  ConfigWrtParam derivWrtParam (freeTimeTraj, 0, t);
	  checkGradientAndThrow (derivWrtParam, 0,
				 freeTimeTraj.parameters ());
	}
      catch (BadGradient& bg)
	{
	  std::cout << bg << std::endl;
	}

      std::cout << fmterDeriv << std::endl;
    }

  std::cout << std::endl;
}

int run_test ()
{
  Spline::vector_t params (5);

  // Scale.
  params[0] = 1.;
  // Initial position.
  params[1] = 0.;
  // Control point 1.
  params[2] = 25.;
  // Control point 2.
  params[3] = 75.;
  // Final position.
  params[4] = 100.;

  // Make trajectories.
  Spline::interval_t timeRange = Spline::makeInterval (0., 4.);
  Spline spline (timeRange, 1, removeScaleFromParameters (params), "before");
  FreeTimeTrajectory<Spline::derivabilityOrder> freeTimeTraj (spline, 1.);

  assert (freeTimeTraj.inputSize () == 1);
  assert (freeTimeTraj.outputSize () == 1);

  printTable (spline, freeTimeTraj);

  params[0] = .5;
  freeTimeTraj.setParameters (params);
  printTable (spline, freeTimeTraj);

  params[0] = 2.;
  freeTimeTraj.setParameters (params);
  printTable (spline, freeTimeTraj);

  return 0;
}

GENERATE_TEST ()
