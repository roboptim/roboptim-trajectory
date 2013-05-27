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

#include <roboptim/core/io.hh>
#include <roboptim/core/finite-difference-gradient.hh>
#include <roboptim/core/visualization/gnuplot.hh>

#include <roboptim/trajectory/frontal-speed.hh>
#include <roboptim/trajectory/orthogonal-speed.hh>
#include <roboptim/trajectory/cubic-b-spline.hh>
#include <roboptim/trajectory/state-function.hh>

using namespace roboptim;
using namespace roboptim::visualization;

int run_test ()
{
  const unsigned orderMax = 1;
  CubicBSpline::vector_t params (24);

  // Initial position.
  params[0] = 0.,  params[1] = 0., params[2] = 0.;
  params[3] = 0.,  params[4] = 0., params[5] = 0.;
  params[6] = 0.,  params[7] = 0., params[8] = 0.;
  // Control point 3.
  params[9] = 25.,  params[10] = 100., params[11] = 0.;
  // Control point 4.
  params[12] = 75.,  params[13] = 0., params[14] = 0.;
  // Final position.
  params[15] = 100., params[16] = 100., params[17] = 0.;
  params[18] = 100., params[19] = 100., params[20] = 0.;
  params[21] = 100., params[22] = 100., params[23] = 0.;

  CubicBSpline::interval_t timeRange = CubicBSpline::makeInterval (0., 4.);
  CubicBSpline spline (timeRange, 3, params, "before");

  for (unsigned i = 0; i < 10; ++i)
    {
      const StableTimePoint timePoint = i / 10. * tMax;
      const double t = timePoint.getTime (spline.timeRange ());

      boost::shared_ptr<DerivableFunction> frontalSpeed (new FrontalSpeed ());
      StateFunction<CubicBSpline> stateFunction
	(spline, frontalSpeed, timePoint, orderMax);

      boost::shared_ptr<DerivableFunction> orthogonalSpeed
	(new OrthogonalSpeed ());
      StateFunction<CubicBSpline> orthoStateFunction
	(spline, orthogonalSpeed, timePoint, orderMax);

      std::cout << "State cost evaluation:" << std::endl
		<< normalize (stateFunction (params)) << std::endl
		<< "State cost gradient:" << std::endl
		<< normalize (stateFunction.gradient (params)) << std::endl
		<< "Trajectory state (splitted):" << std::endl;
      for (unsigned o = 0; o <= orderMax; ++o)
	std::cout << normalize (spline.derivative (t, o)) << std::endl;
      std::cout << "Trajectory state (one call):" << std::endl
		<< normalize (spline.state (t, orderMax)) << std::endl
		<< "Trajectory state variation (splitted):" << std::endl;
      for (unsigned o = 0; o <= orderMax; ++o)
	std::cout << normalize (spline.variationDerivWrtParam (t, o))
		  << std::endl;
      std::cout << "Trajectory state (one call):" << std::endl
		<< normalize (spline.variationStateWrtParam (t, orderMax))
		<< std::endl;

      try
	{
	  std::cout << "Check frontal speed gradient." << std::endl;
	  checkGradientAndThrow (*frontalSpeed, 0, spline.state (t, orderMax));

	  std::cout << "Check orthogonal speed gradient." << std::endl;
	  checkGradientAndThrow (*orthogonalSpeed, 0,
				 spline.state (t, orderMax));

	  std::cout << "Check state cost gradient." << std::endl;
	  checkGradientAndThrow (stateFunction, 0, params);
	}
      catch (BadGradient<EigenMatrixDense>& bg)
	{
	  std::cout << bg << std::endl;
	}
      std::cout << "\n\n" << std::endl;
    }

  return 0;
}

GENERATE_TEST ()
