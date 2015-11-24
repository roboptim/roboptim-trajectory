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

#ifndef ROBOPTIM_TRAJECTORY_PROBLEM_OVER_SPLINES_FACTORY_HXX
# define ROBOPTIM_TRAJECTORY_PROBLEM_OVER_SPLINES_FACTORY_HXX

# include <stdexcept>

# include <boost/format.hpp>
# include <boost/make_shared.hpp>

namespace roboptim
{
  namespace trajectory
  {
    template <typename T, typename S>
    ProblemOverSplinesFactory<T,S>::ProblemOverSplinesFactory
    (const splines_t& splines, const problem_t& problem, CostType cost,
     value_type scaling)
      : splines_ (splines),
	dimension_ (),
	problem_ (),
	t0_ (),
	tmax_ (),
	epsilon_ (0.),
	jerkFactory_ (),
	constraints_ (),
	objScaling_ (scaling)
    {
      int dimension = splines[0]->dimension();
      std::pair<value_type, value_type> timerange = splines[0]->timeRange();

      for (size_t i = 0; i < splines.size(); ++i)
	{
	  if (splines[i]->dimension () != dimension
              || splines[i]->timeRange () != timerange)
	    throw std::runtime_error ("splines are not comparable");
	}

      dimension_ = static_cast<size_type> (dimension);

      t0_ = timerange.first;
      tmax_ = timerange.second;

      if (cost == COST_JERK)
      {
        initializeJerkFactory ();
        problem_ = boost::make_shared<problem_t> (jerkFactory_->getJerk ());
      }
      else
      {
        problem_ = boost::make_shared<problem_t> (problem);
      }
    }

    template <typename T, typename S>
    typename ProblemOverSplinesFactory<T, S>::value_type
    ProblemOverSplinesFactory<T, S>::t0 () const
    {
      return t0_;
    }

    template <typename T, typename S>
    typename ProblemOverSplinesFactory<T, S>::value_type&
    ProblemOverSplinesFactory<T, S>::t0 ()
    {
      return t0_;
    }

    template <typename T, typename S>
    typename ProblemOverSplinesFactory<T, S>::value_type
    ProblemOverSplinesFactory<T, S>::tmax () const
    {
      return tmax_;
    }

    template <typename T, typename S>
    typename ProblemOverSplinesFactory<T, S>::value_type&
    ProblemOverSplinesFactory<T, S>::tmax ()
    {
      return tmax_;
    }

    template <typename T, typename S>
    typename ProblemOverSplinesFactory<T, S>::value_type
    ProblemOverSplinesFactory<T, S>::epsilon () const
    {
      return epsilon_;
    }

    template <typename T, typename S>
    typename ProblemOverSplinesFactory<T, S>::value_type&
    ProblemOverSplinesFactory<T, S>::epsilon ()
    {
      return epsilon_;
    }

    template <typename T, typename S>
    const typename ProblemOverSplinesFactory<T, S>::problem_t&
    ProblemOverSplinesFactory<T, S>::problem () const
    {
      return *problem_;
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::addSpline (const S& spline)
    {
      if (static_cast<size_type> (spline.dimension ()) == dimension_
          && spline.timeRange().first == t0_
          && spline.timeRange().second == tmax_)
	{
	  splines_.push_back(boost::make_shared<S>(spline));
	}
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::updateRange (const interval_t& newRange,
						       CostType cost)
    {
      assert(newRange.first >= t0_);
      assert(newRange.first < newRange.second);

      t0_ = newRange.first;
      tmax_ = newRange.second;

      updateProblem (cost);
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::updateStartingPoint
    (value_type startingPoint, CostType cost)
    {
      updateRange(S::makeInterval(startingPoint, tmax_), cost);
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::updateEndingPoint
    (value_type endingPoint, CostType cost)
    {
      updateRange(S::makeInterval(t0_, endingPoint), cost);
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::addConstraint
    (value_type time, int order, const vector_t& values,
     const scaling_t& scaling, value_type eps)
    {
      if (values.size () != static_cast<size_type> (splines_.size ()))
        throw std::range_error ("invalid number of spline values");

      if (eps < 0.) eps = std::abs (epsilon_);

      constraints_.resize(constraints_.size() + 1);
      constraints_.back().first = time;

      for (size_t i = 0; i < splines_.size(); ++i)
	{
          size_type ii = static_cast<size_type> (i);
	  constraints_.back().second.push_back
            (boost::make_tuple (localConstraint (time, order, values[ii], i),
                                S::makeInterval (-eps, eps),
                                scaling[i]));

	  if (time <= tmax_ && time >= t0_)
	    problem_->addConstraint
              (boost::get<0> (boost::get<freeze_t>
                              (constraints_.back ().second.back ())),
               S::makeInterval (-eps, eps), scaling[i]);
	}
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::addConstraint
    (value_type startingPoint, int order, const intervals_t& range,
     const scalingVect_t& scaling)
    {
      size_type inputSize = problem_->function ().inputSize ();

      constraints_.resize (constraints_.size () + 1);
      constraints_.back ().first = startingPoint;

      intervals_t ranges (2);

      for (size_t i = 0; i < splines_.size (); ++i)
	{
	  ranges[0] = S::makeLowerInterval (range[i].first);
	  ranges[1] = S::makeUpperInterval (range[i].second);

	  constraints_.back ().second.push_back
            (boost::make_tuple (boost::make_shared<splinesConstraint_t>
                                (splines_, i, order, startingPoint, inputSize),
                                ranges, scaling[i]));

	  if (startingPoint < tmax_ && startingPoint >= t0_)
	    problem_->addConstraint
              (boost::get<0> (boost::get<globalConstraint_t>
                              (constraints_.back ().second.back ())),
               ranges, scaling[i]);
	}
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::addConstraint
    (value_type time, int order, const vector_t& values)
    {
      scaling_t scaling (static_cast<size_t> (values.size ()), 1.);
      addConstraint (time, order, values, scaling);
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::addConstraint
    (value_type startingPoint, int order, const intervals_t& range)
    {
      scalingVect_t scaling (range.size (), scaling_t (2, 1.));
      addConstraint (startingPoint, order, range, scaling);
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::updateProblem (CostType cost)
    {
      boost::shared_ptr<problem_t> pb;

      if (cost == COST_JERK)
	{
	  if (!jerkFactory_) initializeJerkFactory ();
	  jerkFactory_->updateRange (S::makeInterval (t0_, tmax_));
	  problem_ = boost::make_shared<problem_t> (jerkFactory_->getJerk ());
	}
      else
	{
	  // Reinitialize the constraints.
	  problem_->clearConstraints ();
	}

      globalConstraint_t g;
      freeze_t f;
      for (typename problemConstraints_t::iterator
	     c = constraints_.begin (); c != constraints_.end (); ++c)
        if (c->first >= t0_ && c->first <= tmax_)
          for (typename std::vector<supportedConstraint_t>::iterator
		 ci  = c->second.begin (); ci != c->second.end (); ++ci)
	    {
	      switch (ci->which())
		{
		case 0 :
		  g = boost::get<globalConstraint_t> (*ci);
		  problem_->addConstraint (boost::get<0> (g), boost::get<1> (g),
					   boost::get<2> (g));
		  break;
		case 1 :
		  f = boost::get<freeze_t> (*ci);
		  problem_->addConstraint (boost::get<0>(f), boost::get<1>(f),
					   boost::get<2>(f));
		  break;
		}
	    }
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::initializeJerkFactory ()
    {
      jerkFactory_ = boost::make_shared<JerkOverSplinesFactory<S, T> >
        (splines_, S::makeInterval (t0_, tmax_), objScaling_);
    }

    template <typename T, typename S>
    const typename ProblemOverSplinesFactory<T, S>::numericLinearConstraintPtr_t
    ProblemOverSplinesFactory<T, S>::localConstraint (value_type t,
						      int order,
						      value_type value,
                                                      size_t spline_idx)
    {
      size_type inputSize = problem_->function ().inputSize ();
      matrix_t A (1, static_cast<typename matrix_t::Index> (inputSize));
      vector_t B (1);
      A.setZero();
      B.coeffRef(0) = -value;
      int splineOffset = 0;

      const splinePtr_t& spline = splines_[spline_idx];

      for (size_t i = 0; i < spline_idx; ++i)
        splineOffset += static_cast<int> (splines_[i]->getNumberControlPoints ());

      size_type int_idx = static_cast<size_type> (spline->interval (t));
      for (int j = 0; j < spline->dimension() + 1; ++j)
	{
	  int col = splineOffset + static_cast<int> (int_idx) - j;
	  size_t k = static_cast<size_t> (int_idx - j);
	  A.coeffRef (0, col) = spline->
	    basisPolynomials()[k][static_cast<size_t> (j)].derivative (t, order);
	}

      std::string name;
      if (order == 0)
        name = (boost::format ("f(%.3f) = %.3f (%s)")
		% t % value % spline->getName ()).str ();
      else
        name = (boost::format ("f^(%i) (%.3f) = %.3f (%s)")
		% order % t % value % spline->getName ()).str ();

      return boost::make_shared<numericLinearConstraint_t> (A, B, name);
    }
  }
}

#endif //! ROBOPTIM_TRAJECTORY_PROBLEM_OVER_SPLINES_FACTORY_HXX
