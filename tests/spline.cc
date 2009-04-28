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

#include "common.hh"

#include <roboptim-trajectory/fwd.hh>
#include <roboptim-trajectory/spline.hh>

using namespace roboptim;

int run_test ()
{
  Spline::vector_t params (8);

  // Initial position.
  params[0] = 1.,  params[1] = 2.;
  // Control point 1.
  params[2] = 0.,  params[3] = 6.;
  // Control point 2.
  params[4] = 4.,  params[5] = 5.;
  // Final position.
  params[6] = 3., params[7] = 4.;

  Spline spline (std::make_pair (0., 5.), 2, params);

  std::cout << spline << std::endl;

  std::cout << "Values:" << std::endl;
  std::cout << spline (0.) << std::endl;
  std::cout << spline (2.5) << std::endl;
  std::cout << spline (5.) << std::endl;

  std::cout << "1st derivative:" << std::endl;
  std::cout << spline.derivative (0., 1) << std::endl;
  std::cout << spline.derivative (2.5, 1) << std::endl;
  std::cout << spline.derivative (5., 1) << std::endl;

  std::cout << "2nd derivative:" << std::endl;
  std::cout << spline.derivative (0., 2) << std::endl;
  std::cout << spline.derivative (2.5, 2) << std::endl;
  std::cout << spline.derivative (5., 2) << std::endl;

  return 0;
}

GENERATE_TEST ()
