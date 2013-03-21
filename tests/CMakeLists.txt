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

# Generate local-libdir.hh.
CONFIG_FILES(shared-tests/local-libdir.hh)

# Make sure local headers are found
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/tests)


# Add Boost path to include directories.
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

# add roboptim core to include directories
INCLUDE_DIRECTORIES(${ROBOPTIM_CORE_INCLUDE_DIRS})

# Make Boost.Test generates the main function in test cases.
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

SET(COMMON_HEADERS
  ${HEADERS}
  shared-tests/common.hh
  shared-tests/hs071.hh
  )

# ROBOPTIM_TRAJECTORY_TEST(NAME)
# ------------------------
#
# Define a test named `NAME'.
#
# This macro will create a binary from `NAME.cc', link it against
# Boost and add it to the test suite.
#
MACRO(ROBOPTIM_TRAJECTORY_TEST NAME)
  ADD_EXECUTABLE(${NAME} ${NAME}.cc ${COMMON_HEADERS})
  ADD_TEST(${NAME} ${RUNTIME_OUTPUT_DIRECTORY}/${NAME})

  TARGET_LINK_LIBRARIES(${NAME} roboptim-trajectory)

  # Link against Boost.
  TARGET_LINK_LIBRARIES(${NAME} ${Boost_LIBRARIES})

  PKG_CONFIG_USE_DEPENDENCY(${NAME} roboptim-core)
ENDMACRO(ROBOPTIM_TRAJECTORY_TEST)

# ROBOPTIM_TRAJECTORY_ANTHROPOMORPHIC_TEST(NAME)
# ----------------------------------------------
#
# Define a test named `NAME' and link it wit the
# anthropomorphic cost function.
#
# This macro will create a binary from `NAME.cc', link it against
# Boost and add it to the test suite.
#
MACRO(ROBOPTIM_TRAJECTORY_ANTHROPOMORPHIC_TEST NAME)
  ADD_EXECUTABLE(${NAME} ${NAME}.cc
    ${COMMON_HEADERS} anthropomorphic-cost-function.cc)
  ADD_TEST(${NAME} ${RUNTIME_OUTPUT_DIRECTORY}/${NAME})

  TARGET_LINK_LIBRARIES(${NAME} roboptim-trajectory)

  # Link against Boost.
  TARGET_LINK_LIBRARIES(${NAME} ${Boost_LIBRARIES})

  PKG_CONFIG_USE_DEPENDENCY(${NAME} roboptim-core)
ENDMACRO(ROBOPTIM_TRAJECTORY_ANTHROPOMORPHIC_TEST)

# Tools.
ROBOPTIM_TRAJECTORY_TEST(simple)
ROBOPTIM_TRAJECTORY_TEST(free-time-trajectory)
ROBOPTIM_TRAJECTORY_TEST(free-time-trajectory-stable-time-point)
ROBOPTIM_TRAJECTORY_TEST(stable-point-state-function)
ROBOPTIM_TRAJECTORY_TEST(state-function)

# Spline.
ROBOPTIM_TRAJECTORY_TEST(spline-bicycle-optim)
ROBOPTIM_TRAJECTORY_TEST(spline)
ROBOPTIM_TRAJECTORY_TEST(spline-gradient)
ROBOPTIM_TRAJECTORY_TEST(spline-optimization)
ROBOPTIM_TRAJECTORY_TEST(spline-time-optimization)
ROBOPTIM_TRAJECTORY_TEST(polynomial)
ROBOPTIM_TRAJECTORY_TEST(b-spline)
ROBOPTIM_TRAJECTORY_TEST(constraint-b-spline)


# Anthropomorphic cost function.
ROBOPTIM_TRAJECTORY_TEST(anthropomorphic-cost-function-case-1)
ROBOPTIM_TRAJECTORY_TEST(anthropomorphic-cost-function-case-2)
ROBOPTIM_TRAJECTORY_TEST(anthropomorphic-cost-function-case-3)
ROBOPTIM_TRAJECTORY_TEST(anthropomorphic-cost-function-case-4)
ROBOPTIM_TRAJECTORY_TEST(anthropomorphic-cost-function-case-5)
ROBOPTIM_TRAJECTORY_TEST(anthropomorphic-cost-function-case-6)
