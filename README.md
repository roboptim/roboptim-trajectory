roboptim-core
=============

[![Build Status](https://travis-ci.org/roboptim/roboptim-trajectory.png?branch=master)](https://travis-ci.org/roboptim/roboptim-trajectory)

For general information about the project, please refer to its
homepage: https://github.com/laas/roboptim-trajectory

Optimization problems in test suite
-----------------------------------

To run all tests, you must to have at least *one* solver plug-in
supporting non-linear problems (DifferentiableFunction as cost
function and (LinearFunction, DifferentiableFunction) as constraints).

Valid examples are: cfsqp, ipopt.

You can then choose which solver will be used through CMake:

    $ cmake <usual cmake flags> -DTESTSUITE_SOLVER=ipopt

By default, no solver is chosen and these tests are not compiled nor
run.

Before sending a pull request, please make sure that these tests are
succeeding.
