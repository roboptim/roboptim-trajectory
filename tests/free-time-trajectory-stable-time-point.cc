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

#include <boost/format.hpp>
#include <boost/scoped_ptr.hpp>

#include <roboptim/core/finite-difference-gradient.hh>
#include <roboptim/core/io.hh>
#include <roboptim/core/util.hh>
#include <roboptim/core/visualization/gnuplot.hh>

#include <roboptim/trajectory/free-time-trajectory.hh>
#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/cubic-b-spline.hh>
#include <roboptim/trajectory/trajectory-cost.hh>

#include "shared-tests/common.hh"

using boost::format;
using boost::io::group;

using namespace roboptim;
using namespace roboptim::visualization;


typedef FreeTimeTrajectory<CubicBSpline> freeTime_t;

template <typename T>
bool isAlmostEqual (const T& x, const T& y, const T& epsilon = 1e10-8)
{
  return (x - y) * (x - y) < epsilon;
}

struct ConfigWrtParam : public DerivableFunction
{
  ConfigWrtParam (const freeTime_t& traj, StableTimePoint stp) throw ()
    : DerivableFunction (traj.parameters ().size (),
			 traj.outputSize (),
			 "config wrt param"),
      traj_ (traj),
      stp_ (stp)
  {
  }

  ~ConfigWrtParam () throw ()
  {
  }

  void
  impl_compute (result_t& res, const argument_t& p) const throw ()
  {
    boost::scoped_ptr<freeTime_t> updatedTrajectory (traj_.clone ());
    updatedTrajectory->setParameters (p);
    res = (*updatedTrajectory) (stp_);
  }

  void
  impl_gradient (gradient_t& grad, const argument_t& p, size_type)
    const throw ()
  {
    boost::scoped_ptr<freeTime_t> updatedTrajectory (traj_.clone ());
    updatedTrajectory->setParameters (p);
    matrix_t tmp = updatedTrajectory->variationDerivWrtParam (stp_, 0);
    grad = row (tmp, 0);
  }

  const freeTime_t& traj_;
  StableTimePoint stp_;
};

struct DerivWrtParam : public DerivableFunction
{
  DerivWrtParam (const freeTime_t& traj) throw ()
    : DerivableFunction (traj.inputSize (),
			 traj.parameters ().size (),
			 "config wrt param"),
      traj_ (traj.clone ())
  {}

  ~DerivWrtParam () throw ()
  {}

  void
  impl_compute (result_t& res, const argument_t& stp) const throw ()
  {
    // Make sure one never evaluates under zero.
    // This case happens because of the finite difference gradient
    // checking!
    value_type alpha = stp[0];
    if (alpha <= 0.)
      alpha = 0.;
    if (alpha >= 1.)
      alpha = 1.;


    matrix_t tmp = traj_->variationDerivWrtParam (alpha * tMax, 0);
    res = row (tmp, 0);
  }

  void
  impl_gradient (gradient_t& grad, const argument_t& stp, size_type i)
    const throw ()
  {
    matrix_t tmp = traj_->variationDerivWrtParam (stp[0] * tMax, 1);
    grad[0] = row (tmp, 0)[i];
  }

  const boost::scoped_ptr<const freeTime_t> traj_;
};


void printTable (const CubicBSpline& spline, const freeTime_t& freeTimeTraj);

void printTable (const CubicBSpline& spline, const freeTime_t& freeTimeTraj)
{
  double tmin = Function::getLowerBound (spline.timeRange ());
  double tmax = Function::getUpperBound (spline.timeRange ());
  double fttTmin = Function::getLowerBound (freeTimeTraj.timeRange ());
  double fttTmax = Function::getUpperBound (freeTimeTraj.timeRange ());

  std::cout << format ("Spline range: [%1%, %2%]") % tmin % tmax
	    << std::endl
	    << format ("FTT range: [%1%, %2%]") % fttTmin % fttTmax
	    << std::endl
	    << format ("FTT scale: %1%") % freeTimeTraj.timeScale ()
	    << std::endl;

  format fmterTitle ("| %1% %|6t||| %2% %|20t|| %3% %|35t|||"
		     " %4% %|55t|| %5% %|79t||");

  format fmter ("| %1.2f %|6t||| %2.2f %|20t|| %3.2f %|35t|||"
		" %4.2f %|55t|| %5.2f %|79t||");

  std::cout << "/---------------------------------------"
	    << "---------------------------------------\\" << std::endl
	    << fmterTitle
    % "T" % "Value" % "Value (ftt)"
    % "Derivative" % "Derivative (ftt)"
	    << std::endl
	    << "----------------------------------------"
	    << "----------------------------------------"
	    << std::endl;

  for (double alpha = 0.; alpha <= 1.; alpha += .1)
    {
      if (alpha >= 1.)
	alpha = 1.;

      fmter % alpha;
      fmter % normalize (spline (alpha * tMax)[0]);
      fmter % normalize (freeTimeTraj (alpha * tMax)[0]);
      fmter % normalize (spline.derivative (alpha * tMax, 1)[0]);
      fmter % normalize (freeTimeTraj.derivative (alpha * tMax, 1)[0]);
      std:: cout << fmter << std::endl;
    }
  std::cout << "\\---------------------------------------"
	    << "---------------------------------------/"
	    << std::endl << std::endl;


  std::cout << "Gradient check." << std::endl;
  for (double alpha = 0.; alpha <= 1.; alpha += .1)
    {
      if (alpha >= 1.)
	alpha = 1.;

      CubicBSpline::vector_t x (1);
      x[0] = alpha;

      try
	{
	  std::cout << "Spline gradient." << std::endl;
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
  format fmterConfig ("%1.2f %|50t|%2.2f");
  for (double alpha = 0.; alpha <= 1.; alpha += .1)
    {
      if (alpha >= 1.)
	alpha = 1.;

      CubicBSpline::jacobian_t splineVarConfig =
	spline.variationConfigWrtParam (alpha * tMax);
      fmterConfig % normalize (splineVarConfig);

      freeTime_t::jacobian_t fttVarConfig =
	freeTimeTraj.variationConfigWrtParam (alpha * tMax);
      fmterConfig % normalize (fttVarConfig);

      try
	{
	  ConfigWrtParam configWrtParam (freeTimeTraj, alpha * tMax);
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
  format fmterDeriv ("%1.2f %|50t|%2.2f");
  for (double alpha = 0. + .1; alpha <= 1. - .1; alpha += .1)
    {
      if (alpha >= 1.)
	alpha = 1.;
      if (alpha <= 0.)
	alpha = 0.;

      CubicBSpline::jacobian_t splineVarDeriv =
	spline.variationDerivWrtParam (alpha * tMax, 1);
      fmterDeriv % normalize (splineVarDeriv);

      freeTime_t::jacobian_t fttVarDeriv =
	freeTimeTraj.variationDerivWrtParam (alpha * tMax, 1);
      fmterDeriv % normalize (fttVarDeriv);

      for (unsigned gradientId = 0;
	   gradientId < freeTimeTraj.parameters ().size (); ++gradientId)
	try
	  {
	    CubicBSpline::vector_t t_ (1);
	    t_[0] = alpha;
	    DerivWrtParam derivWrtParam (freeTimeTraj);
	    checkGradientAndThrow (derivWrtParam, gradientId, t_);
	  }
	catch (BadGradient& bg)
	  {
	    std::cout << "Component "
		      << gradientId
		      << " of the gradient is wrong:" << std::endl;
	    std::cout << bg << std::endl;
	  }

      std::cout << fmterDeriv << std::endl;
    }

  std::cout << std::endl;
}

int run_test ()
{
  typedef CubicBSpline::value_type value_type;
  CubicBSpline::vector_t params (9);

  // Scale.
  params[0] = 1.;
  // Initial position.
  params[1] = 0.;
  params[2] = 0.;
  params[3] = 0.;
  // Control point 1.
  params[4] = 25.;
  // Control point 2.
  params[5] = 75.;
  // Final position.
  params[6] = 100.;
  params[7] = 100.;
  params[8] = 100.;

  // Make trajectories.
  CubicBSpline::interval_t timeRange = CubicBSpline::makeInterval (0., 4.);
  CubicBSpline spline (timeRange, 1, removeScaleFromParameters (params),
		       "before");
  FreeTimeTrajectory<CubicBSpline> freeTimeTraj (spline, 1.);

  assert (freeTimeTraj.inputSize () == 1);
  assert (freeTimeTraj.outputSize () == 1);

  params[0] = 1.;
  freeTimeTraj.setParameters (params);
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
