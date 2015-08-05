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

# Which solver should be used?
SET(TESTSUITE_SOLVER ""
  CACHE STRING "Which solver should be used during tests?")
ADD_DEFINITIONS(-DTESTSUITE_SOLVER="${TESTSUITE_SOLVER}")

# Which linear solver should Ipopt use?
SET(IPOPT_LINEAR_SOLVER "mumps"
  CACHE STRING "Which linear solver should be used by Ipopt during tests?")
ADD_DEFINITIONS(-DIPOPT_LINEAR_SOLVER="${IPOPT_LINEAR_SOLVER}")

# Whether to run all the (potentially slow) tests or not
SET(FULL_TEST "OFF"
  CACHE STRING "Whether to run all the (potentially slow) tests or not")

# Include shared test suite.
SET(PLUGIN_PATH "${CMAKE_BINARY_DIR}/src")
SET(SOLVER_NAME ${TESTSUITE_SOLVER})
SET(FUNCTION_TYPE ::roboptim::EigenMatrixDense)
SET(PROGRAM_SUFFIX "")
SET(COST_FUNCTION_TYPE ::roboptim::GenericDifferentiableFunction)
SET(CONSTRAINT_TYPE_1 ::roboptim::GenericLinearFunction)
SET(CONSTRAINT_TYPE_2 ::roboptim::GenericDifferentiableFunction)

# Make sure local headers are found
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/tests)

# Add Boost path to include directories.
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

# Make Boost.Test generates the main function in test cases.
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

SET(COMMON_HEADERS
  ${HEADERS}
  shared-tests/common.hh
  shared-tests/fixture.hh
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
  SET_TARGET_PROPERTIES(${NAME}
    PROPERTIES COMPILE_FLAGS
    "-DSOLVER_NAME='\"${SOLVER_NAME}\"' -DPLUGIN_PATH='\"${PLUGIN_PATH}\"' -DFUNCTION_TYPE=${FUNCTION_TYPE} -DCOST_FUNCTION_TYPE=${COST_FUNCTION_TYPE} ${CONSTRAINT_TYPE_1_STR} ${CONSTRAINT_TYPE_2_STR} -DTESTS_DATA_DIR='\"${CMAKE_CURRENT_SOURCE_DIR}\"' -DLOG_FILENAME='\"${NAME}.log\"' ")

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
MACRO(ROBOPTIM_TRAJECTORY_TEST NAME)
  ADD_EXECUTABLE(${NAME} ${NAME}.cc ${COMMON_HEADERS})
  SET_TARGET_PROPERTIES(${NAME}
    PROPERTIES COMPILE_FLAGS
    "-DSOLVER_NAME='\"${SOLVER_NAME}\"' -DPLUGIN_PATH='\"${PLUGIN_PATH}\"' -DFUNCTION_TYPE=${FUNCTION_TYPE} -DCOST_FUNCTION_TYPE=${COST_FUNCTION_TYPE} ${CONSTRAINT_TYPE_1_STR} ${CONSTRAINT_TYPE_2_STR} -DTESTS_DATA_DIR='\"${CMAKE_CURRENT_SOURCE_DIR}\"' -DLOG_FILENAME='\"${NAME}.log\"' ")

  ADD_TEST(${NAME} ${RUNTIME_OUTPUT_DIRECTORY}/${NAME})

  TARGET_LINK_LIBRARIES(${NAME} roboptim-trajectory)

  # Link against Boost.
  TARGET_LINK_LIBRARIES(${NAME} ${Boost_LIBRARIES} boost_serialization)

  PKG_CONFIG_USE_DEPENDENCY(${NAME} roboptim-core)
ENDMACRO(ROBOPTIM_TRAJECTORY_TEST)

# Tools.
ROBOPTIM_TRAJECTORY_TEST(simple)
ROBOPTIM_TRAJECTORY_TEST(polynomial-3)
ROBOPTIM_TRAJECTORY_TEST(polynomial)
ROBOPTIM_TRAJECTORY_TEST(free-time-trajectory)
ROBOPTIM_TRAJECTORY_TEST(free-time-trajectory-stable-time-point)
ROBOPTIM_TRAJECTORY_TEST(stable-point-state-function)
ROBOPTIM_TRAJECTORY_TEST(state-function)

# Spline.
ROBOPTIM_TRAJECTORY_TEST(cubic-b-spline)
CONFIGURE_FILE ("cubic-b-spline-gradient.dat" "cubic-b-spline-gradient.dat")
ROBOPTIM_TRAJECTORY_TEST(cubic-b-spline-2)
CONFIGURE_FILE ("cubic-b-spline-2.dat" "cubic-b-spline-2.dat")
ROBOPTIM_TRAJECTORY_TEST(cubic-b-spline-gradient)
ROBOPTIM_TRAJECTORY_TEST(b-spline)
ROBOPTIM_TRAJECTORY_TEST(constrained-b-spline)
ROBOPTIM_TRAJECTORY_TEST(spline-length)
ROBOPTIM_TRAJECTORY_TEST(b-spline-matplotlib)
ROBOPTIM_TRAJECTORY_TEST(cubic-b-spline-matplotlib)

# Other
ROBOPTIM_TRAJECTORY_TEST(vector-interpolation)

IF(TESTSUITE_SOLVER)
  # TODO: fix these tests and re-add them
  ROBOPTIM_TRAJECTORY_TEST(problem-over-splines)
  #ROBOPTIM_TRAJECTORY_TEST(spline-bicycle-optim)
  #ROBOPTIM_TRAJECTORY_TEST(spline-optimization)
  ROBOPTIM_TRAJECTORY_TEST(spline-time-optimization)
  CONFIGURE_FILE ("spline-time-optimization.dat" "spline-time-optimization.dat")
ENDIF()

# Anthropomorphic cost function.
IF(TESTSUITE_SOLVER AND FULL_TEST)
  ROBOPTIM_TRAJECTORY_TEST(anthropomorphic-cost-function-case-1)
  ROBOPTIM_TRAJECTORY_TEST(anthropomorphic-cost-function-case-2)
  ROBOPTIM_TRAJECTORY_TEST(anthropomorphic-cost-function-case-3)
  ROBOPTIM_TRAJECTORY_TEST(anthropomorphic-cost-function-case-4)
  ROBOPTIM_TRAJECTORY_TEST(anthropomorphic-cost-function-case-5)
  ROBOPTIM_TRAJECTORY_TEST(anthropomorphic-cost-function-case-6)
ENDIF()
