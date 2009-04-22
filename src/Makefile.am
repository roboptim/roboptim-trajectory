include $(top_srcdir)/build-aux/init.mk

# Add available warnings flags.
AM_CXXFLAGS += $(WARNING_CXXFLAGS)

lib_LTLIBRARIES = libroboptim-trajectory.la

libroboptim_trajectory_la_SOURCES = 	\
	doc.hh				\
	trajectory.cc

libroboptim_trajectory_la_SOURCES += 	\
	spline/basis_function.cpp	\
	spline/basis_function.h		\
	spline/bspline.cpp		\
	spline/bspline.h


libroboptim_trajectory_la_LIBADD = $(PTHREAD_LIBS) -lroboptim-core
libroboptim_trajectory_la_CXXFLAGS = $(PTHREAD_CFLAGS)
