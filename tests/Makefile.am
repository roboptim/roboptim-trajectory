include $(top_srcdir)/build-aux/init.mk
include $(top_srcdir)/build-aux/autotest.mk

# Add RobOptim flags.
AM_CPPFLAGS += $(ROBOPTIMCORE_CFLAGS)

# Add shared headers to include path.
AM_CPPFLAGS += -I$(srcdir)/shared-tests

# Add available warnings flags.
AM_CXXFLAGS += $(WARNING_CXXFLAGS)

# Package the following files.
EXTRA_DIST += package.m4.in

# Clean atconfig (do not distribute it).
CLEANFILES += atconfig

# Source which are shared by all tests.
COMMON_SOURCES = shared-tests/common.hh

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

# free-time-trajectory
check_PROGRAMS += free-time-trajectory
free_time_trajectory_SOURCES = free-time-trajectory.cc $(COMMON_SOURCES)
free_time_trajectory_LDADD = $(top_builddir)/src/libroboptim-trajectory.la

# free-time-trajectory-stable-time-point
check_PROGRAMS += free-time-trajectory-stable-time-point
free_time_trajectory_stable_time_point_SOURCES =	\
	free-time-trajectory-stable-time-point.cc 	\
	$(COMMON_SOURCES)
free_time_trajectory_stable_time_point_LDADD = 		\
	$(top_builddir)/src/libroboptim-trajectory.la


# state-function
check_PROGRAMS += state-function
state_function_SOURCES = state-function.cc $(COMMON_SOURCES)
state_function_LDADD = $(top_builddir)/src/libroboptim-trajectory.la

# stable-point-state-function
check_PROGRAMS += stable-point-state-function
stable_point_state_function_SOURCES = stable-point-state-function.cc \
				     $(COMMON_SOURCES)
stable_point_state_function_LDADD = $(top_builddir)/src/libroboptim-trajectory.la


#########################################
# Tests which require the CFSQP solver. #
#########################################
if HAVE_ROBOPTIM_CORE_CFSQP_PLUGIN

# spline-optimization
check_PROGRAMS += spline-optimization
spline_optimization_SOURCES = spline-optimization.cc $(COMMON_SOURCES)
spline_optimization_LDADD = $(top_builddir)/src/libroboptim-trajectory.la

# spline-time-optimization
check_PROGRAMS += spline-time-optimization
spline_time_optimization_SOURCES = spline-time-optimization.cc $(COMMON_SOURCES)
spline_time_optimization_LDADD = $(top_builddir)/src/libroboptim-trajectory.la

# anthropomorphic-cost-function-case-1
check_PROGRAMS += anthropomorphic-cost-function-case-1
anthropomorphic_cost_function_case_1_SOURCES = 	\
	anthropomorphic-cost-function-case-1.cc \
	anthropomorphic-cost-function.cc 	\
	anthropomorphic-cost-function.hh 	\
	$(COMMON_SOURCES)
anthropomorphic_cost_function_case_1_LDADD = 	\
	$(top_builddir)/src/libroboptim-trajectory.la

# anthropomorphic-cost-function-case-2
check_PROGRAMS += anthropomorphic-cost-function-case-2
anthropomorphic_cost_function_case_2_SOURCES = 	\
	anthropomorphic-cost-function-case-2.cc \
	anthropomorphic-cost-function.cc 	\
	anthropomorphic-cost-function.hh 	\
	$(COMMON_SOURCES)
anthropomorphic_cost_function_case_2_LDADD = 	\
	$(top_builddir)/src/libroboptim-trajectory.la

# anthropomorphic-cost-function-case-3
check_PROGRAMS += anthropomorphic-cost-function-case-3
anthropomorphic_cost_function_case_3_SOURCES = 	\
	anthropomorphic-cost-function-case-3.cc \
	anthropomorphic-cost-function.cc 	\
	anthropomorphic-cost-function.hh 	\
	$(COMMON_SOURCES)
anthropomorphic_cost_function_case_3_LDADD = 	\
	$(top_builddir)/src/libroboptim-trajectory.la

# anthropomorphic-cost-function-case-4
check_PROGRAMS += anthropomorphic-cost-function-case-4
anthropomorphic_cost_function_case_4_SOURCES = 	\
	anthropomorphic-cost-function-case-4.cc \
	anthropomorphic-cost-function.cc 	\
	anthropomorphic-cost-function.hh 	\
	$(COMMON_SOURCES)
anthropomorphic_cost_function_case_4_LDADD = 	\
	$(top_builddir)/src/libroboptim-trajectory.la

# anthropomorphic-cost-function-case-5
check_PROGRAMS += anthropomorphic-cost-function-case-5
anthropomorphic_cost_function_case_5_SOURCES = 	\
	anthropomorphic-cost-function-case-5.cc \
	anthropomorphic-cost-function.cc 	\
	anthropomorphic-cost-function.hh 	\
	$(COMMON_SOURCES)
anthropomorphic_cost_function_case_5_LDADD = 	\
	$(top_builddir)/src/libroboptim-trajectory.la

# anthropomorphic-cost-function-case-6
check_PROGRAMS += anthropomorphic-cost-function-case-6
anthropomorphic_cost_function_case_6_SOURCES = 	\
	anthropomorphic-cost-function-case-6.cc \
	anthropomorphic-cost-function.cc 	\
	anthropomorphic-cost-function.hh 	\
	$(COMMON_SOURCES)
anthropomorphic_cost_function_case_6_LDADD = 	\
	$(top_builddir)/src/libroboptim-trajectory.la

endif


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
EXTRA_DIST +=						\
	anthropomorphic-cost-function-case-1.stdout	\
	anthropomorphic-cost-function-case-2.stdout	\
	anthropomorphic-cost-function-case-3.stdout	\
	anthropomorphic-cost-function-case-4.stdout	\
	anthropomorphic-cost-function-case-5.stdout	\
	anthropomorphic-cost-function-case-6.stdout	\
	free-time-trajectory.stdout			\
	free-time-trajectory-stable-time-point.stdout	\
	simple.stdout					\
	spline-gradient.stdout				\
	spline-optimization.stdout			\
	spline-bicycle-optim.stdout			\
	spline-time-optimization.stdout			\
	spline.stdout					\
	stable-point-state-function.stdout		\
	state-function.stdout
