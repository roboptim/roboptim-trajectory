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

#include <roboptim-core/cfsqp.hh>

#include <roboptim-trajectory/fwd.hh>
#include <roboptim-trajectory/spline.hh>
#include <roboptim-trajectory/sum-cost.hh>
#include <roboptim-trajectory/state-cost.hh>

#include "common.hh"

using namespace roboptim;

struct MyStateCost : public StateCost<Spline>
{
  MyStateCost (size_type m)
    : StateCost<Spline> (m)
  {
  }

  virtual vector_t operator () (const vector_t&) const throw ()
  {
    vector_t res (1);
    res[0] = 0.;
    return res;
  }

  virtual gradient_t gradient (const vector_t&, int) const throw ()
  {
    gradient_t grad (n);
    grad.clear ();
    return grad;
  }
};

int run_test ()
{
  Spline::vector_t params (8);

  // Initial position.
  params[0] = 0.,  params[1] = 0.;
  // Control point 1.
  params[2] = 25.,  params[3] = 50.;
  // Control point 2.
  params[4] = 50.,  params[5] = 25.;
  // Final position.
  params[6] = 100., params[7] = 100.;

  Spline spline (std::make_pair (0., 5.), 2, params);

  std::cout << "# Start generating GNU plot information...." << std::endl
	    << "set term x11 enhanced persist" << std::endl
	    << "plot '-' with line, '-' with line" << std::endl;

  const double step = 0.01;
  for (double t = 0.; t < 5.; t += step)
    {
      Spline::vector_t res = spline (t);
      std::cout << res[0] << " " << res[1] << std::endl;
    }
  std::cout << "e\n" << std::endl;

  // Optimize.
  SumCost<Spline>::vector_t pts (10);
  pts.clear ();

  MyStateCost statecost (4);
  SumCost<Spline> cost (spline, statecost, pts);

  CFSQPSolver::problem_t problem (cost);
  CFSQPSolver solver (problem);

  CFSQPSolver::result_t res = solver.minimum ();
  Result& result = boost::get<Result> (res);

  spline.parameters () = result.x;
  spline.updateParameters ();

  for (double t = 0.; t < 5.; t += step)
    {
      Spline::vector_t res = spline (t);
      std::cout << res[0] << " " << res[1] << std::endl;
    }
  std::cout << "e\n" << std::endl;

  return 0;
}

GENERATE_TEST ()
