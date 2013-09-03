// Copyright (C) 2013 by Thomas Moulard, AIST, CNRS, INRIA.
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

#include <boost/mpl/list.hpp>

#include "shared-tests/fixture.hh"

#include <boost/test/test_case_template.hpp>

#include <iostream>

#include <roboptim/core/io.hh>
#include <roboptim/core/filter/map.hh>
#include <roboptim/core/function/cos.hh>
#include <roboptim/trajectory/vector-interpolation.hh>

using namespace roboptim;


typedef boost::mpl::list< ::roboptim::EigenMatrixDense> functionTypes_t;

BOOST_FIXTURE_TEST_SUITE (core, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE_TEMPLATE (filter_vector_interpolation_test, T, functionTypes_t)
{
  boost::shared_ptr<Cos<T> > cosinus = boost::make_shared<Cos<T> > ();
  boost::shared_ptr<GenericDifferentiableFunction<T> >
    mappedCosinus = map (cosinus, 30);

  typename GenericDifferentiableFunction<T>::argument_t x (30);
  for (unsigned i = 0; i < 30; ++i)
    x.coeffRef (i) = i;
  typename GenericDifferentiableFunction<T>::result_t cosinusVector
    = (*mappedCosinus) (x);

  typename VectorInterpolation::vector_t params
    = typename VectorInterpolation::vector_t (cosinusVector);
  typename VectorInterpolation::size_type outputSize = 1;

  boost::shared_ptr<VectorInterpolation >
    interpolation = vectorInterpolation (params, outputSize);

  BOOST_CHECK (interpolation->outputSize () == 1);

  for (unsigned i = 0; i < 30; ++i)
    {
      std::cout
	<< (*interpolation) (i) << "\n"
	<< interpolation->derivative (i);
    }
}

BOOST_AUTO_TEST_SUITE_END ()
