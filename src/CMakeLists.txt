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

# Define the directory where plug-ins will be installed.
SET(PLUGINDIR lib/${PROJECT_NAME})

# Main library.
ADD_LIBRARY(roboptim-trajectory SHARED
  ${HEADERS}
  doc.hh
  orthogonal-speed.cc
  CMakeLists.txt
  frontal-speed.cc
  spline-length.cc
  debug.hh
  cubic-b-spline.cc
  )
PKG_CONFIG_USE_DEPENDENCY(roboptim-trajectory roboptim-core)
SET_TARGET_PROPERTIES(roboptim-trajectory PROPERTIES SOVERSION 3.0.0)
INSTALL(TARGETS roboptim-trajectory DESTINATION ${CMAKE_INSTALL_LIBDIR})
