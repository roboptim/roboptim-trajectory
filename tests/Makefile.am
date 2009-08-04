include $(top_srcdir)/build-aux/init.mk
include $(top_srcdir)/build-aux/autotest.mk

# Add RobOptim flags.
AM_CPPFLAGS += $(ROBOPTIMCORE_CFLAGS)

# Add available warnings flags.
AM_CXXFLAGS += $(WARNING_CXXFLAGS)

# Package the following files.
EXTRA_DIST += package.m4.in

# Clean atconfig (do not distribute it).
CLEANFILES += atconfig

# Source which are shared by all tests.
COMMON_SOURCES = common.hh

# simple
check_PROGRAMS = simple
simple_SOURCES = simple.cc $(COMMON_SOURCES)
simple_LDADD = $(top_builddir)/src/libroboptim-trajectory.la

# spline
check_PROGRAMS += spline
spline_SOURCES = spline.cc $(COMMON_SOURCES)
spline_LDADD = $(top_builddir)/src/libroboptim-trajectory.la

# spline-gradient
check_PROGRAMS += spline-gradient
spline_gradient_SOURCES = spline-gradient.cc $(COMMON_SOURCES)
spline_gradient_LDADD = $(top_builddir)/src/libroboptim-trajectory.la

# spline-optimization
check_PROGRAMS += spline-optimization
spline_optimization_SOURCES = spline-optimization.cc $(COMMON_SOURCES)
spline_optimization_LDADD = $(top_builddir)/src/libroboptim-trajectory.la

# free-time-trajectory
check_PROGRAMS += free-time-trajectory
free_time_trajectory_SOURCES = free-time-trajectory.cc $(COMMON_SOURCES)
free_time_trajectory_LDADD = $(top_builddir)/src/libroboptim-trajectory.la

# spline-time-optimization
check_PROGRAMS += spline-time-optimization
spline_time_optimization_SOURCES = spline-time-optimization.cc $(COMMON_SOURCES)
spline_time_optimization_LDADD = $(top_builddir)/src/libroboptim-trajectory.la

# complex-cost-function
check_PROGRAMS += complex-cost-function
complex_cost_function_SOURCES = complex-cost-function.cc $(COMMON_SOURCES)
complex_cost_function_LDADD = $(top_builddir)/src/libroboptim-trajectory.la


.PHONY: generate-reference

# Re-generate reference outputs.
generate-reference:
	@echo "Generating references...";	\
	for i in $(check_PROGRAMS); do		\
	 IN="$(srcdir)/$$i.stdin";		\
	 OUT="$(srcdir)/$$i.stdout";		\
	 if test -f "$$IN"; then 		\
	  "./$$i" < $$IN > $$OUT || true;	\
	 else					\
	  "./$$i" > $$OUT || true;		\
	 fi					\
	done


# ------------ #
# Test suite.  #
# ------------ #

# Distribute reference files.
EXTRA_DIST +=					\
	complex-cost-function.stdout		\
	free-time-trajectory.stdout		\
	simple.stdout				\
	spline-gradient.stdout			\
	spline-optimization.stdout		\
	spline-time-optimization.stdout		\
	spline.stdout
