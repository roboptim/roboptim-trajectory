# Copyright 2010, Thomas Moulard, LAAS-CNRS
#
# This file is part of roboptim-core.
# roboptim-core is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# roboptim-core is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Lesser Public License for more details.
# You should have received a copy of the GNU Lesser General Public License
# along with roboptim-core.  If not, see <http://www.gnu.org/licenses/>.

# Requires at least CMake 2.6 to configure the package.
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/cpack.cmake)

SET(PROJECT_NAME roboptim-trajectory)
SET(PROJECT_DESCRIPTION "RobOptim trajectory")
SET(PROJECT_URL "http://github.com/roboptim/roboptim-trajectory")

# Use MathJax for Doxygen formulae
SET(DOXYGEN_USE_MATHJAX "YES")

SET(CSD ${CMAKE_SOURCE_DIR})
SET(HEADERS
  ${CSD}/include/roboptim/trajectory/anthropomorphic-cost-function.hh
  ${CSD}/include/roboptim/trajectory/anthropomorphic-cost-function.hxx
  ${CSD}/include/roboptim/trajectory/b-spline.hh
  ${CSD}/include/roboptim/trajectory/b-spline.hxx
  ${CSD}/include/roboptim/trajectory/constrained-b-spline.hh
  ${CSD}/include/roboptim/trajectory/constrained-b-spline.hxx
  ${CSD}/include/roboptim/trajectory/constraints-over-splines.hh
  ${CSD}/include/roboptim/trajectory/constraints-over-splines.hxx
  ${CSD}/include/roboptim/trajectory/cubic-b-spline.hh
  ${CSD}/include/roboptim/trajectory/cubic-b-spline.hxx
  ${CSD}/include/roboptim/trajectory/free-time-trajectory.hh
  ${CSD}/include/roboptim/trajectory/free-time-trajectory.hxx
  ${CSD}/include/roboptim/trajectory/freeze.hh
  ${CSD}/include/roboptim/trajectory/freeze.hxx
  ${CSD}/include/roboptim/trajectory/frontal-speed.hh
  ${CSD}/include/roboptim/trajectory/fwd.hh
  ${CSD}/include/roboptim/trajectory/jerk-over-splines-factory.hh
  ${CSD}/include/roboptim/trajectory/jerk-over-splines-factory.hxx
  ${CSD}/include/roboptim/trajectory/limit-omega.hh
  ${CSD}/include/roboptim/trajectory/limit-omega.hxx
  ${CSD}/include/roboptim/trajectory/limit-speed.hh
  ${CSD}/include/roboptim/trajectory/limit-speed.hxx
  ${CSD}/include/roboptim/trajectory/orthogonal-speed.hh
  ${CSD}/include/roboptim/trajectory/polynomial.hh
  ${CSD}/include/roboptim/trajectory/polynomial.hxx
  ${CSD}/include/roboptim/trajectory/polynomial-3.hh
  ${CSD}/include/roboptim/trajectory/problem-over-splines-factory.hh
  ${CSD}/include/roboptim/trajectory/problem-over-splines-factory.hxx
  ${CSD}/include/roboptim/trajectory/spline-length.hh
  ${CSD}/include/roboptim/trajectory/stable-point-state-function.hh
  ${CSD}/include/roboptim/trajectory/stable-point-state-function.hxx
  ${CSD}/include/roboptim/trajectory/stable-time-point.hh
  ${CSD}/include/roboptim/trajectory/state-function.hh
  ${CSD}/include/roboptim/trajectory/state-function.hxx
  ${CSD}/include/roboptim/trajectory/sys.hh
  ${CSD}/include/roboptim/trajectory/trajectory-cost.hh
  ${CSD}/include/roboptim/trajectory/trajectory-cost.hxx
  ${CSD}/include/roboptim/trajectory/trajectory-sum-cost.hh
  ${CSD}/include/roboptim/trajectory/trajectory-sum-cost.hxx
  ${CSD}/include/roboptim/trajectory/trajectory.hh
  ${CSD}/include/roboptim/trajectory/trajectory.hxx
  ${CSD}/include/roboptim/trajectory/vector-interpolation.hh
  ${CSD}/include/roboptim/trajectory/vector-interpolation.hxx
  ${CSD}/include/roboptim/trajectory/visualization/limit-speed.hh
  ${CSD}/include/roboptim/trajectory/visualization/speed.hh
  ${CSD}/include/roboptim/trajectory/visualization/trajectory.hh
  ${CSD}/include/roboptim/trajectory/visualization/matplotlib.hh
  ${CSD}/include/roboptim/trajectory/visualization/b-spline-matplotlib.hh
  ${CSD}/include/roboptim/trajectory/visualization/cubic-b-spline-matplotlib.hh
  ${CSD}/include/roboptim/trajectory.hh
)

SETUP_PROJECT()

# Add main library to pkg-config file.
PKG_CONFIG_APPEND_LIBS(roboptim-trajectory)

# Search for dependencies.
SET(BOOST_COMPONENTS
  filesystem system thread program_options unit_test_framework serialization)
SEARCH_FOR_BOOST()

# Libtool dynamic loading
# This project does not use Libtool directly but still uses ltdl for
# plug-in loading.
INCLUDE(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(ltdl.h LTDL_H_FOUND)
IF (NOT LTDL_H_FOUND)
  MESSAGE(FATAL_ERROR
    "Failed to find ltdl.h, check that Libtool ltdl is installed.")
ENDIF()
#FIXME: check that libltdl.so is available.

ADD_REQUIRED_DEPENDENCY("roboptim-core >= 3.1")
ADD_REQUIRED_DEPENDENCY("eigen3 >= 3.2.0")

HEADER_INSTALL("${HEADERS}")

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tests)

SETUP_PROJECT_FINALIZE()
SETUP_PROJECT_CPACK()
