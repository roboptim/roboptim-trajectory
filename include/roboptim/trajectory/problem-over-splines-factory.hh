// Copyright (C) 2015 by Félix Darricau, EPITA, AIST, CNRS.
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

#ifndef ROBOPTIM_TRAJECTORY_PROBLEM_OVER_SPLINES_FACTORY_HH
# define ROBOPTIM_TRAJECTORY_PROBLEM_OVER_SPLINES_FACTORY_HH

# include <vector>

# include <boost/tuple/tuple.hpp>
# include <boost/shared_ptr.hpp>
# include <boost/variant.hpp>

# include <roboptim/core/problem.hh>
# include <roboptim/core/numeric-linear-function.hh>
# include <roboptim/trajectory/constraints-over-splines.hh>
# include <roboptim/trajectory/jerk-over-splines-factory.hh>

namespace roboptim
{
  namespace trajectory
  {
    /// \brief Factory of problems over splines.
    ///
    /// \tparam T Matrix type.
    /// \tparam S Spline type.
    template <typename T, typename S>
    class ProblemOverSplinesFactory
    {
      typedef typename S::value_type value_type;
      typedef typename S::interval_t interval_t;

      typedef Problem<T> problem_t;
      typedef typename problem_t::scaling_t scaling_t;
      typedef typename problem_t::scalingVect_t scalingVect_t;
      typedef typename problem_t::intervals_t intervals_t;
      typedef typename problem_t::function_t function_t;
      typedef typename problem_t::constraint_t constraint_t;
      typedef std::vector<constraint_t> constraints_t;
      typedef typename function_t::vector_t vector_t;
      typedef typename function_t::matrix_t matrix_t;

      typedef ConstraintsOverSplines<T, S> splinesConstraint_t;
      typedef boost::shared_ptr<splinesConstraint_t> splinesConstraintPtr_t;
      typedef boost::tuple<splinesConstraintPtr_t, intervals_t, scaling_t>
      globalConstraint_t;

      typedef GenericNumericLinearFunction<T> numericLinearConstraint_t;
      typedef boost::shared_ptr<numericLinearConstraint_t> numericLinearConstraintPtr_t;
      typedef boost::tuple<numericLinearConstraintPtr_t, interval_t, value_type>
      freeze_t;

      typedef boost::variant<globalConstraint_t, freeze_t> supportedConstraint_t;

      typedef std::vector<boost::shared_ptr<S> > splines_t;

      /// The constraints are stored with their given startingPoint/time.
      typedef std::vector<std::pair<value_type, std::vector<supportedConstraint_t> > >
      problemConstraints_t;

    public:
      ProblemOverSplinesFactory (const splines_t& splines, problem_t& problem);

      value_type t0 () const;
      value_type& t0 ();

      value_type tmax () const;
      value_type& tmax ();

      /// \brief Updates the range of the optimization problem
      ///
      /// \param newRange desired timeRange
      /// \param buildCostFunction whether we want to build a new Jerk
      void updateRange(interval_t newRange, bool buildCostFunction = true);

      /// \brief Updates the startingPoint of the problem, using updateRange
      void updateStartingPoint(value_type startingPoint, bool buildCostFunction = true);

      /// \brief Updates the endingPoint of the problem, using updateRange
      void updateEndingPoint(value_type endingPoint, bool buildCostFunction = true);

      /// \brief Adds a spline to the problem
      ///
      /// Warning, the size of the existing cost function will not change,
      /// it can therefore break your optimization problem.
      ///
      /// It is best practice to give all the splines to the constructor.
      void addSpline(S& spline);

      /// \brief Adds freezing contraints on every spline at a given time
      ///
      /// These constraints are equality constraints on the value, derivative
      /// or second derivative (depending on the given order) of the spline 
      /// at t=time. It constraints every spline.
      ///
      /// \param time Time of the constraint
      /// \param order Derivation order (0 for no derivation)
      /// \param value Goals of the constraints for each spline
      /// \param scaling Scalings of the constraints for each spline
      void addConstraint(value_type time, int order, std::vector<value_type> value, scaling_t scaling);

      /// \brief Adds contraints on every spline starting at a given time
      ///
      /// These constraints are bounding constraints on the value, derivative
      /// or second derivative (depending on the given order) of the spline, on
      /// the time range between the given startingPoint and the end of the
      /// interval. It constraints every spline.
      ///
      /// \param startingPoint Time when the constraint interval starts
      /// \param order Derivation order (0 for no derivation)
      /// \param range Bounds of the constraint for each spline
      /// \param scaling Scalings of the constraints for each spline
      void addConstraint(value_type startingPoint, int order, intervals_t range, scalingVect_t scaling);

      /// \brief Adds freezing contraints on every spline at a given time
      ///
      /// Calls the corresponding addConstraint with a scaling set to 1 for each
      /// spline.
      void addConstraint(value_type time, int order, std::vector<value_type> value);

      /// \brief Adds contraints on every spline starting at a given time
      ///
      /// Calls the corresponding addConstraint with a scaling set to 1 for each
      /// spline.
      void addConstraint(value_type startingPoint, int order, intervals_t range);

      problem_t& getProblem() const
      {
        return *problem_;
      }

    private:
      /// \brief Generate a new problem, given the actual range and constraints
      ///
      /// This function is called when the time range is updated.
      /// It creates a new problem with, either on a new Jerk cost function or
      /// on the previous one (determined by buildCostFunction) and adds to it
      /// the constraints corresponding to the new range.
      ///
      /// \param buildCostFunction whether we want to build a new Jerk
      void updateProblem(bool buildCostFunction);

      /// \brief Creates and retrieves a new equality constraint on a spline
      ///
      /// \returns shared pointer to the constraint
      /// \param time Time of the constraint
      /// \param order Derivation order (0 for no derivation)
      /// \param value Goal of the constraint
      /// \param spline Index of the spline in the vector
      const boost::shared_ptr<GenericNumericLinearFunction<T> > localConstraint(value_type time, int order, value_type value, unsigned long spline);

      /// \brief Shared pointers to the splines used in the problem
      splines_t splines_;

      /// \brief Dimension of the problem
      unsigned long dimension_;

      /// \brief Shared pointer to the built problem
      boost::shared_ptr<problem_t> problem_;

      /// \brief Starting point of the problem
      value_type t0_;

      /// \brief Ending point of the problem
      value_type tmax_;

      /// \brief Intervals used to get the true range for bounds constraints
      intervals_t range_;

      /// \brief Input size of the cost function
      long inputsize_;

      /// \brief Jerk cost function factory
      /// Default behaviour is to use it, since when the range is updated, the
      /// factory also updates the jerk. But if the user wants to provide its
      /// own cost function, it will be ignored
      boost::shared_ptr<JerkOverSplinesFactory<S, T> > factory_;

      /// \brief Constraints of the problem
      problemConstraints_t constraints_;
    };
  }
}

#include <roboptim/trajectory/problem-over-splines-factory.hxx>

#endif //! ROBOPTIM_TRAJECTORY_PROBLEM_OVER_SPLINES_FACTORY_HH
