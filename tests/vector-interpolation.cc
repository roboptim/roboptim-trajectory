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
#include <fstream>

#include <boost/mpl/list.hpp>

#include "shared-tests/fixture.hh"

#include <boost/test/test_case_template.hpp>

#include <iostream>

#include <roboptim/core/io.hh>
#include <roboptim/core/alloc.hh>
#include <roboptim/core/operator/derivative.hh>
#include <roboptim/core/operator/map.hh>
#include <roboptim/core/operator/selection.hh>
#include <roboptim/core/function/cos.hh>
#include <roboptim/core/visualization/gnuplot.hh>
#include <roboptim/core/visualization/gnuplot-commands.hh>
#include <roboptim/core/visualization/gnuplot-function.hh>
#include <roboptim/trajectory/vector-interpolation.hh>

using namespace roboptim;
using namespace roboptim::trajectory;
using namespace roboptim::visualization;


struct VectorInterpolationDerivWrtParameters : public DifferentiableFunction
{
  VectorInterpolationDerivWrtParameters
  (const VectorInterpolation& vectorInterpolation, value_type t, size_type variableId)
    : DifferentiableFunction
      (1,
       vectorInterpolation.outputSize (),
       "vectorInterpolation differentiable w.r.t parameters"),
      vectorInterpolation_ (vectorInterpolation),
      t_ (t),
      variableId_ (variableId)
  {}

  virtual void
  impl_compute (result_ref result, const_argument_ref x)
    const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    bool cur_malloc_allowed = is_malloc_allowed ();
    set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    VectorInterpolation vectorInterpolation (vectorInterpolation_);
    vector_t params = vectorInterpolation.parameters ();
    params[variableId_] = x[0];
    vectorInterpolation.setParameters (params);
    result = vectorInterpolation (t_);

#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    set_is_malloc_allowed (cur_malloc_allowed);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION
  }

  virtual void
  impl_gradient (gradient_ref gradient,
		 const_argument_ref x,
		 size_type functionId = 0)
    const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    bool cur_malloc_allowed = is_malloc_allowed ();
    set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    VectorInterpolation vectorInterpolation (vectorInterpolation_);
    vector_t params = vectorInterpolation.parameters ();
    params[variableId_] = x[0];
    vectorInterpolation.setParameters (params);

    gradient = vectorInterpolation.variationConfigWrtParam (t_).block
      (functionId, variableId_, 1, 1);

#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
    set_is_malloc_allowed (cur_malloc_allowed);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION
  }

private:
  const VectorInterpolation& vectorInterpolation_;
  value_type t_;
  size_type variableId_;
};


typedef boost::mpl::list< ::roboptim::EigenMatrixDense> functionTypes_t;

BOOST_FIXTURE_TEST_SUITE (core, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE_TEMPLATE (filter_vector_interpolation_test, T,
			       functionTypes_t)
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

  // Display initial and final trajectory.
  VectorInterpolation::discreteInterval_t intervalS (0., 10., 0.01);
  using namespace roboptim::visualization::gnuplot;
  // Fix clash with std::set
  using roboptim::visualization::gnuplot::set;

  std::ofstream f ("/tmp/vector-interpolation-0.gp");
  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();
  f  << (gnuplot
	 << set ("multiplot layout 4, 1")
	 << plot (*cosinus, intervalS)
	 << plot (*derivative(cosinus, 0), intervalS)
	 << plot (*interpolation, intervalS)
	 << plot (*derivative(interpolation, 0), intervalS)
	 << unset ("multiplot")
	 );

}

BOOST_AUTO_TEST_CASE_TEMPLATE (filter_vector_interpolation_nonscalar_test,
			       T, functionTypes_t)
{
  typename VectorInterpolation::value_type dt = 1.;
  typename VectorInterpolation::size_type outputSize = 2;
  typename VectorInterpolation::vector_t params (6 * outputSize);

  params <<
    0. , 0., // t = 0
    1. , 1., // t = 1
    2. , 2., // t = 2
    3. , 1., // t = 3
    4. , 0., // t = 4
    5. , 1.; // t = 5

  boost::shared_ptr<VectorInterpolation >
    interpolation = vectorInterpolation (params, outputSize, dt);

  VectorInterpolation::discreteInterval_t intervalS (0., 10., 0.01);

  BOOST_CHECK (interpolation->outputSize () == outputSize);

  for (VectorInterpolation::vector_t::Index i = 0; i < 6; ++i)
    {
      VectorInterpolation::value_type pt =
	static_cast<VectorInterpolation::value_type> (i) * dt;
      BOOST_CHECK (params[i * outputSize + 0] == (*interpolation) (pt)[0]);
      BOOST_CHECK (params[i * outputSize + 1] == (*interpolation) (pt)[1]);
    }

  // Display initial and final trajectory.
  using namespace roboptim::visualization::gnuplot;
  // Fix clash with std::set
  using roboptim::visualization::gnuplot::set;

  std::ofstream f ("/tmp/vector-interpolation-1.gp");
  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();

  f << (gnuplot
	<< set ("multiplot layout 3, 1")
	<< plot (*interpolation, intervalS)
	<< plot (*derivative(interpolation, 0), intervalS)
	<< unset ("multiplot")
	);
}

BOOST_AUTO_TEST_CASE_TEMPLATE (filter_vector_derivWrtParams_test, T,
			       functionTypes_t)
{
  typename VectorInterpolation::value_type dt = 1.;
  typename VectorInterpolation::size_type outputSize = 2;
  typename VectorInterpolation::vector_t params (6 * outputSize);

  params <<
    0. , 0., // t = 0
    1. , 1., // t = 1
    2. , 2., // t = 2
    3. , 1., // t = 3
    4. , 0., // t = 4
    5. , 1.; // t = 5

  boost::shared_ptr<VectorInterpolation >
    interpolation = vectorInterpolation (params, outputSize, dt);

  boost::shared_ptr<VectorInterpolation >
    interpolationCopy =
    boost::shared_ptr<VectorInterpolation> (interpolation->clone ());
  for (typename VectorInterpolation::vector_t::Index i = 0;
       i < params.size (); ++i)
    BOOST_CHECK_SMALL
      (interpolation->parameters ()[i] -
       interpolationCopy->parameters ()[i],
       1e-8);

  boost::shared_ptr<VectorInterpolationDerivWrtParameters> derivWrtParams
    = boost::make_shared<VectorInterpolationDerivWrtParameters>
    (*interpolation, 1.5, 2);


  VectorInterpolation::discreteInterval_t intervalS (0., 10., 0.01);

  // Display initial and final trajectory.
  using namespace roboptim::visualization::gnuplot;
  // Fix clash with std::set
  using roboptim::visualization::gnuplot::set;

  std::ofstream f ("/tmp/vector-interpolation-2.gp");
  Gnuplot gnuplot = Gnuplot::make_interactive_gnuplot ();

  f << (gnuplot
	<< set ("multiplot layout 3, 1")
	<< plot (*derivWrtParams, intervalS)
	<< plot (*derivative(derivWrtParams, 0), intervalS)
	<< unset ("multiplot")
	);

}

BOOST_AUTO_TEST_SUITE_END ()
