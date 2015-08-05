// Copyright (C) 2015 by Félix Darricau, AIST, CNRS, EPITA.
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

#include <iostream>
#include <boost/filesystem/fstream.hpp>

#include <roboptim/core/io.hh>

#include "shared-tests/fixture.hh"

#include <roboptim/core/function/constant.hh>

#include <roboptim/trajectory/cubic-b-spline.hh>
#include <roboptim/trajectory/b-spline.hh>
#include <roboptim/trajectory/jerk-over-splines-factory.hh>
#include <roboptim/trajectory/problem-over-splines-factory.hh>

#include <roboptim/trajectory/visualization/cubic-b-spline-matplotlib.hh>
#include <roboptim/trajectory/visualization/b-spline-matplotlib.hh>

#include <roboptim/trajectory/fwd.hh>

#include <roboptim/core/visualization/matplotlib-function.hh>

using namespace roboptim;
using namespace roboptim::visualization;
using namespace roboptim::visualization::matplotlib;
using namespace roboptim::trajectory;
using namespace roboptim::trajectory::visualization::matplotlib;

typedef boost::mpl::list< ::roboptim::trajectory::CubicBSpline,
			  ::roboptim::trajectory::BSpline<3> > splinesType_t;

template <typename T>
void splineName(std::string& str);

template <>
void splineName<roboptim::trajectory::CubicBSpline>(std::string& str)
{
  str = "cubic-b-spline";
}

template <>
void splineName<roboptim::trajectory::BSpline<3> >(std::string& str)
{
  str = "b-spline";
}

BOOST_FIXTURE_TEST_SUITE (trajectory, TestSuiteConfiguration)

BOOST_AUTO_TEST_CASE_TEMPLATE (problem_over_splines, spline_t, splinesType_t)
{
  std::cout.precision(10);
  typedef EigenMatrixSparse T;
  typedef Solver<GenericDifferentiableFunction<T>,
		 boost::mpl::vector<GenericLinearFunction<T>,
				    GenericDifferentiableFunction<T> > > solver_t;
  typedef typename spline_t::vector_t param_t;
  typedef GenericConstantFunction<T>::value_type value_type;
  std::string spline_name;
  splineName<spline_t>(spline_name);
  boost::filesystem::ofstream file("problem-over-splines-" + spline_name + ".py");
  value_type step = 0.005;
  param_t params(13);
  params << 40,25,5,0,8,15,0,10,3,-5,15,22,35;
  boost::shared_ptr<spline_t> spline = boost::make_shared<spline_t>(std::make_pair(0,1), 1, params, "First spline", false);
  std::vector<boost::shared_ptr<spline_t> > splines;
  splines.push_back(spline);
  param_t params2(13);
  params2 << 10,12,9,8,25,15,24,20,1,5,12,-12,3;
  boost::shared_ptr<spline_t> spline2 = boost::make_shared<spline_t>(std::make_pair(0,1), 1, params2, "Second spline", true);
  splines.push_back(spline2);
  solver_t::problem_t::intervals_t range;
  std::vector<value_type> range2;
  JerkOverSplinesFactory<spline_t, T> jerkFactory(splines, GenericConstantFunction<T>::makeInterval(0, 1));
  solver_t::problem_t pb (*jerkFactory.getJerk());
  ProblemOverSplinesFactory<T, spline_t> constraints(splines, pb);
  BOOST_CHECK(pb.constraints().size() == 0);
  constraints.updateStartingPoint(0.02);
  range.clear();
  range.push_back(std::make_pair<value_type, value_type>(0, 5));
  range.push_back(std::make_pair<value_type, value_type>(0, 5));
  constraints.addConstraint(0.02, 0, range);
  range.clear();
  range.push_back(std::make_pair<value_type, value_type>(1, 10));
  range.push_back(std::make_pair<value_type, value_type>(3, 8));
  constraints.addConstraint(0.62, 1, range);
  BOOST_CHECK(constraints.getProblem().constraints().size() == 4);
  solver_t::problem_t problem (constraints.getProblem());

  param_t startingPoint(26);
  startingPoint << params, params2;
  problem.startingPoint() = startingPoint;
  SolverFactory<solver_t> factory ("ipopt-sparse", problem);

  solver_t& solver = factory ();

  std::cout << solver << std::endl;

  solver.parameters()["ipopt.tol"].value = 1e-5;
  solver.parameters()["ipopt.output_file"].value = std::string("test.log");
  solver.parameters()["ipopt.print_level"].value = 12;
  typename solver_t::result_t res = solver.minimum();
  Matplotlib matplotlib = Matplotlib::make_matplotlib (std::make_pair(3, 2));
  (file)
    << (matplotlib
        << plot_spline (*spline, step)
        << title ("initial B-spline 1")
        << plot_spline (*spline2, step)
        << title ("initial B-spline 2")
        ) << std::endl;
  switch (res.which ())
    {
    case solver_t::SOLVER_VALUE:
      {
        // Get the result.
        roboptim::Result& result = boost::get<roboptim::Result> (res);

        spline->setParameters(result.x.segment(0, 13));
        spline2->setParameters(result.x.segment(13, 13));
	file
	  << (matplotlib
	      << plot_spline (*spline, step)
	      << title ("initial B-spline 1 after first optimization pass")
	      << plot_spline (*spline2, step)
	      << title ("initial B-spline 2 after first optimization pass")
	      ) << std::endl;

        std::cout << result << std::endl;
        break;
      }

    case solver_t::SOLVER_VALUE_WARNINGS:
      {
        // Get the result.
        roboptim::ResultWithWarnings& result = boost::get<roboptim::ResultWithWarnings> (res);

        spline->setParameters(result.x.segment(0, 13));
        spline2->setParameters(result.x.segment(13, 13));
	file
	  << (matplotlib
	      << plot_spline (*spline, step)
	      << title ("initial B-spline 1 after first optimization pass")
	      << plot_spline (*spline2, step)
	      << title ("initial B-spline 2 after first optimization pass")
	      ) << std::endl;

        std::cout << result << std::endl;
        break;
      }

    case solver_t::SOLVER_NO_SOLUTION:
    case solver_t::SOLVER_ERROR:
      {
        roboptim::SolverError err = boost::get<roboptim::SolverError> (res);
        std::cout << "A solution should have been found. Failing..."
		  << std::endl
		  << err.what ()
		  << std::endl;
        roboptim::Result result = boost::get<roboptim::Result> (err.lastState());
        std::cout << result << std::endl;

        BOOST_CHECK(false);
      }
    }
  constraints.updateStartingPoint(0.32);
  range2.clear();
  range2.push_back((*spline)(0.32)[0]);
  range2.push_back((*spline2)(0.32)[0]);
  constraints.addConstraint(0.32, 0, range2);
  range2.clear();
  range2.push_back(spline->derivative(0.32, 1)[0]);
  range2.push_back(spline2->derivative(0.32, 1)[0]);
  constraints.addConstraint(0.32, 1, range2);
  range2.clear();
  range2.push_back(spline->derivative(0.32, 2)[0]);
  range2.push_back(spline2->derivative(0.32, 2)[0]);
  constraints.addConstraint(0.32, 2, range2);
  solver_t::problem_t newproblem (constraints.getProblem());
  BOOST_CHECK(newproblem.constraints().size() == 8);
  startingPoint << spline->parameters(), spline2->parameters();
  newproblem.startingPoint() = startingPoint;
  SolverFactory<solver_t> newfactory ("ipopt-sparse", newproblem);
  solver_t& newsolver = newfactory ();

  std::cout << newsolver << std::endl;

  newsolver.parameters()["ipopt.tol"].value = 1e-5;
  newsolver.parameters()["ipopt.output_file"].value = std::string("test2.log");
  res = newsolver.minimum();
  switch (res.which ())
    {
    case solver_t::SOLVER_VALUE:
      {
        // Get the result.
        roboptim::Result& result = boost::get<roboptim::Result> (res);

        spline->setParameters(result.x.segment(0, 13));
        spline2->setParameters(result.x.segment(13, 13));
	file
	  << (matplotlib
	      << plot_spline (*spline, step)
	      << title ("initial B-spline 1 after second optimization pass")
	      << plot_spline (*spline2, step)
	      << title ("initial B-spline 2 after second optimization pass")
	      ) << std::endl;

        std::cout << result << std::endl;
        break;
      }

    case solver_t::SOLVER_VALUE_WARNINGS:
      {
        // Get the result.
        roboptim::ResultWithWarnings& result = boost::get<roboptim::ResultWithWarnings> (res);

        spline->setParameters(result.x.segment(0, 13));
        spline2->setParameters(result.x.segment(13, 13));
	file
	  << (matplotlib
	      << plot_spline (*spline, step)
	      << title ("initial B-spline 1 after first optimization pass")
	      << plot_spline (*spline2, step)
	      << title ("initial B-spline 2 after first optimization pass")
	      ) << std::endl;

        std::cout << result << std::endl;
        break;
      }

    case solver_t::SOLVER_NO_SOLUTION:
    case solver_t::SOLVER_ERROR:
      {
        roboptim::SolverError err = boost::get<roboptim::SolverError> (res);
        std::cout << "A solution should have been found. Failing..."
		  << std::endl
		  << err.what ()
		  << std::endl;
        roboptim::Result result = boost::get<roboptim::Result> (err.lastState());
        std::cout << result << std::endl;

        BOOST_CHECK(false);
      }
    }
}
BOOST_AUTO_TEST_SUITE_END ()
