// Copyright (C) 2015 by Benjamin Chrétien, CNRS-LIRMM.
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

#include "shared-tests/fixture.hh"

#include <boost/make_shared.hpp>
#include <boost/filesystem/fstream.hpp>

#include <roboptim/core/io.hh>
#include <roboptim/core/util.hh>
#include <roboptim/core/decorator/finite-difference-gradient.hh>
#include <roboptim/core/visualization/matplotlib-function.hh>

#include <roboptim/trajectory/cubic-b-spline.hh>
#include <roboptim/trajectory/b-spline.hh>
#include <roboptim/trajectory/constraints-over-splines.hh>

#include <roboptim/trajectory/visualization/cubic-b-spline-matplotlib.hh>
#include <roboptim/trajectory/visualization/b-spline-matplotlib.hh>

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::matplotlib;
using namespace roboptim::trajectory;
using namespace roboptim::trajectory::visualization::matplotlib;

typedef boost::mpl::list<CubicBSpline, BSpline<3> > splinesType_t;

template <typename T>
std::string splineName();

template <>
std::string splineName<roboptim::trajectory::CubicBSpline> ()
{
  return "cubic-b-spline";
}

template <>
std::string splineName<roboptim::trajectory::BSpline<3> > ()
{
  return "b-spline";
}

template <typename T>
std::string getOutputFilename ()
{
  std::string s = "constraints-over-splines-";
  s += splineName<T> ();
  return s;
}

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE_TEMPLATE (trajectory_constraints_over_splines,
                               spline_t, splinesType_t)
{
  boost::shared_ptr<boost::test_tools::output_test_stream>
    output = retrievePattern ("constraints-over-splines");

  typedef ConstraintsOverSplines<EigenMatrixDense, spline_t>
    constraintsOverSplines_t;
  typedef typename constraintsOverSplines_t::splines_t splines_t;
  typedef typename spline_t::size_type size_type;

  const size_type n = 10;

  typename spline_t::vector_t x (n);
  x << 10., 5., -3., 5., 4., -2., 8., 4., 1., 0.;

  typename spline_t::interval_t timeRange = CubicBSpline::makeInterval (0., 1.);
  boost::shared_ptr<spline_t> spline = boost::make_shared<spline_t>
                                         (timeRange, 1, x, "spline");
  splines_t splines;
  splines.push_back (spline);

  for (double t = 1e-6; t < 1.; t += 1./7.)
  {
    boost::shared_ptr<constraintsOverSplines_t>
      f = boost::make_shared<constraintsOverSplines_t> (splines, 0, 0, t, n);

    typedef typename finiteDifferenceGradientPolicies::Simple<EigenMatrixDense>
      fdRule_t;
    GenericFiniteDifferenceGradient<EigenMatrixDense, fdRule_t> fd (f);

    (*output)
      << x << "\n"
      << *f << "\n"
      << (*f) (x) << "\n"
      << normalize (toDense (f->gradient (x, 0))) << "\n"
      << normalize (toDense (f->gradient (x, 1))) << "\n"
      << normalize (toDense (f->jacobian (x))) << std::endl;

    // Numerical errors between spline implementations prevent basic string
    // comparisons with FD.
    std::cout
      << fd << "\n"
      << fd (x) << "\n"
      << normalize (toDense (fd.gradient (x, 0)), 1e-6) << "\n"
      << normalize (toDense (fd.gradient (x, 1)), 1e-6) << "\n"
      << normalize (toDense (fd.jacobian (x)), 1e-6) << std::endl;

    BOOST_CHECK (allclose (toDense (f->jacobian (x)),
          toDense (fd.jacobian (x)),
          1e-4, 1e-4));
  }

  Matplotlib matplotlib = Matplotlib::make_matplotlib ();
  matplotlib << plot_spline (*splines[0], 0.005);

  boost::filesystem::ofstream pythonPlot
    (getOutputFilename<spline_t> () + ".py");
  pythonPlot << matplotlib << std::endl;

  std::cout << output->str () << std::endl;
  BOOST_CHECK (output->match_pattern ());
}

BOOST_AUTO_TEST_SUITE_END ()
