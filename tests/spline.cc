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

#include "shared-tests/common.hh"

#include <boost/numeric/ublas/io.hpp>

#include <roboptim/core/finite-difference-gradient.hh>

#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>

#include <roboptim/trajectory/fwd.hh>
#include <roboptim/trajectory/cubic-b-spline.hh>

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::gnuplot;

struct SplineDerivWrtParameters : public DifferentiableFunction
{
  SplineDerivWrtParameters (const CubicBSpline& spline, value_type t)
    : DifferentiableFunction
      (spline.parameters ().size (), spline.outputSize (),
       "spline differentiable w.r.t parameters"),
      spline_ (spline),
      t_ (t)
  {}

  virtual void
  impl_compute (result_t& result, const argument_t& x)
    const throw ()
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    CubicBSpline spline (spline_);
    spline.setParameters (x);
    result = spline (t_);
  }

  virtual void
  impl_gradient (gradient_t& gradient,
		 const argument_t& x,
		 size_type functionId = 0)
    const throw ()
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    CubicBSpline spline (spline_);
    spline.setParameters (x);
    gradient = spline.variationConfigWrtParam (t_).row (functionId);
  }

private:
  const CubicBSpline& spline_;
  value_type t_;
};

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE (trajectory_spline)
{
  CubicBSpline::vector_t params (16);

  // Initial position.
  params[0] = 0.,  params[1] = 0.;
  params[2] = 0.,  params[3] = 0.;
  params[4] = 0.,  params[5] = 0.;
  // Control point 3.
  params[6] = 25.,  params[7] = 50.;
  // Control point 4.
  params[8] = 50.,  params[9] = 25.;
  // Final position.
  params[10] = 100., params[11] = 100.;
  params[12] = 100., params[13] = 100.;
  params[14] = 100., params[15] = 100.;

  //FIXME: change interval.
  CubicBSpline spline (std::make_pair (0., 5.), 2, params, "spline");

  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();
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
    << (gnuplot << plot_xy (spline, interval));

  // Reset parameters and recompute some values.
  params[6] = 50.,  params[7] = 25.;
  spline.setParameters (params);

  std::cout
    << "# 1st derivative:" << std::endl
    << "# " << spline.derivative (0., 1) << std::endl
    << "# " << spline.derivative (2.5, 1) << std::endl
    << "# " << spline.derivative (5., 1) << std::endl

    << "# 2nd derivative:" << std::endl
    << "# " << spline.derivative (0., 2) << std::endl
    << "# " << spline.derivative (2.5, 2) << std::endl
    << "# " << spline.derivative (5., 2) << std::endl;

  for (double t = 0.5; t < 5.; t += 0.5)
    {
      try
	{
	  Function::vector_t x (1);
	  x[0] = t;
	  checkGradientAndThrow (spline, 0, x);

	  SplineDerivWrtParameters splineDerivWrtParams (spline, t);
	  checkGradientAndThrow
	    (splineDerivWrtParams, 0, spline.parameters ());
	}
      catch (BadGradient<EigenMatrixDense>& bg)
	{
	  std::cerr << bg << std::endl;
	  BOOST_CHECK(false);
	}
    }

  for (double x = 0.; x < 10.; x += 0.25)
    {
      Function::vector_t params = spline.parameters ();
      params[0 * 2] = 321;
      params[1 * 2] = 123;
      params[2 * 2] =
	6. * (x - 1./6. * params[0 * 2] - 2. / 3. * params[1 * 2]);

      assert (1. / 6. * params[0 * 2]
	      + 2. / 3. * params[1 * 2]
	      + 1. / 6. * params[2 * 2] - x < 1e-8);

      spline.setParameters (params);
      if (std::fabs (spline (0.)[0] - x) >= 1e-8)
	std::cout << "# " << spline (0.)[0] << " != " << x << std::endl;
    }
}

BOOST_AUTO_TEST_SUITE_END ()
