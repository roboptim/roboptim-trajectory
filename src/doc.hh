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

/**

   \mainpage User manual

   \section intro Introduction

   RobOptim trajectory is a toolbox for trajectory optimization.

   It defines cost functions and constraints to simplify problem
   definition. Abstraction for trajectory definition are provided and
   an implementation of cubic splines can be used.


   \section reporting Reporting bugs

   As this package is still in its early development steps, bugs report
   and enhancement proposals are highly welcomed.

   To report a bug, please go this package's
   <a href="http://roboptim.sourceforge.net/">web page</a>.

   To discuss about this package, please use <a
href="https://lists.sourceforge.net/mailman/listinfo/roboptim-user">roboptim-user@lists.sourceforge.net</a>.
*/

/// \namespace roboptim
/// \brief Meta-functions, functions and solvers related classes.


/// \namespace roboptim::visualization
/// \brief Graphic visualization
///
/// Visualization related code. Only Gnuplot is supported currently.


/// \namespace roboptim::visualization::gnuplot
/// \brief Gnuplot rendering
///
/// Gnuplot display classes.


/// \defgroup roboptim_meta_function Mathematical abstract functions
/// \defgroup roboptim_function Mathematical functions
/// \defgroup roboptim_problem Optimization problems
/// \defgroup roboptim_solver Optimization solvers
/// \defgroup roboptim_visualization Visualization
