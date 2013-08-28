roboptim-trajectory
===================

[![Build Status](https://travis-ci.org/roboptim/roboptim-trajectory.png?branch=master)](https://travis-ci.org/roboptim/roboptim-trajectory)
[![Coverage Status](https://coveralls.io/repos/roboptim/roboptim-core/badge.png)](https://coveralls.io/r/roboptim/roboptim-trajectory)

For general information about the project, please refer to its
homepage: https://github.com/roboptim/roboptim-trajectory


How can I install roboptim-core?
--------------------------------

RobOptim uses [CMake](http://www.cmake.org/) to generate build files. For
instance, if you want to build roboptim-core in release with debug info, and
install it in `/my/prefix`, go to the root of the project folder and type:

    $ mkdir -p build && cd build
    $ cmake -DCMAKE_INSTALL_PREFIX=/my/prefix -DCMAKE_BUILD_TYPE=RELWITHDEBINFO ..
    $ make && make test
    $ sudo make install

Please note that roboptim-trajectory requires roboptim-core to be
installed first. If you want to solve optimization problems (or run
the test suite), you will also need at least one solver
package. Recommended solver package is roboptim-core-plugin-ipopt.


Where is the library documentation?
-----------------------------------

This **README** only covers configure/building issues. For more information
regarding this library usage, please refer to the Doxygen documentation.

If you have configured the package as explained in the first section, go
into your `build` directory and type:

    $ make doc

To view the HTML documentation: go in the `doc/doxygen-html` directory
and open `index.html` with your favorite internet browser.


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


How to use Valgrind with the test suite?
----------------------------------------

All the tests launched by the test suite can be prefixed
with the environment variable `CHECK_PREFIX`.

    $ cmake -DCHECK_PREFIX='valgrind --log-file=valgrind.log' ..
    $ make && make check


Available packages
------------------

 * Fedora (Release 0.5):
   https://apps.fedoraproject.org/packages/roboptim-trajectory
