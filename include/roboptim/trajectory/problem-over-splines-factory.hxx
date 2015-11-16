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

namespace roboptim
{
  namespace trajectory
  {
    template <typename T, typename S>
    ProblemOverSplinesFactory<T,S>::ProblemOverSplinesFactory
    (const splines_t& splines, const problem_t& problem)
      : splines_ (splines),
	dimension_ (),
	problem_ (boost::make_shared<problem_t>(problem)),
	t0_ (),
	range_ (2),
	inputsize_ (0),
	epsilon_ (0.),
	factory_ (),
	constraints_ ()
    {
      int dimension = splines[0]->dimension();
      std::pair<value_type, value_type> timerange = splines[0]->timeRange();
      for (unsigned long i = 0; i < splines.size(); ++i)
	{
	  if (splines[i]->dimension() != dimension || splines[i]->timeRange() != timerange)
	    throw std::runtime_error ("Splines are not comparable");
	  else
	    inputsize_ += static_cast<size_type> (splines[i]->getNumberControlPoints ());
	}
      dimension_ = static_cast<unsigned long>(dimension);
      t0_ = timerange.first;
      tmax_ = timerange.second;
      factory_ = boost::make_shared<JerkOverSplinesFactory<S, T> >(splines_, S::makeInterval(t0_, tmax_));
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
    void ProblemOverSplinesFactory<T, S>::addSpline (S& spline)
    {
      if (static_cast<unsigned long>(spline.dimension()) == dimension_
          && spline.timeRange().first == t0_
          && spline.timeRange().second == tmax_)
	{
	  splines_.push_back(boost::make_shared<S>(spline));
	  inputsize_ += static_cast<size_type> (spline.getNumberControlPoints());
	}
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::updateRange (interval_t newRange,
						       bool buildCostFunction)
    {
      assert(newRange.first >= t0_);
      assert(newRange.first < newRange.second);

      t0_ = newRange.first;
      tmax_ = newRange.second;

      updateProblem(buildCostFunction);
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::updateStartingPoint (value_type startingPoint, bool buildCostFunction)
    {
      updateRange(S::makeInterval(startingPoint, tmax_), buildCostFunction);
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::updateEndingPoint (value_type endingPoint, bool buildCostFunction)
    {
      updateRange(S::makeInterval(t0_, endingPoint), buildCostFunction);
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
      assert (problem_->function ().inputSize () == inputsize_);
      constraints_.resize (constraints_.size () + 1);
      constraints_.back ().first = startingPoint;

      for (unsigned long i = 0; i < splines_.size (); ++i)
	{
	  range_[0] = S::makeLowerInterval (range[i].first);
	  range_[1] = S::makeUpperInterval (range[i].second);

	  constraints_.back ().second.push_back
            (boost::make_tuple (boost::make_shared<splinesConstraint_t>
                                (splines_, i, order, startingPoint, inputsize_),
                                range_, scaling[i]));

	  if (startingPoint < tmax_ && startingPoint >= t0_)
	    problem_->addConstraint
              (boost::get<0> (boost::get<globalConstraint_t>
                              (constraints_.back ().second.back ())),
               range_, scaling[i]);
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
    void ProblemOverSplinesFactory<T, S>::updateProblem(bool buildCostFunction)
    {
      boost::shared_ptr<problem_t> pb;
      if (buildCostFunction)
	{
	  factory_->updateRange(S::makeInterval(t0_, tmax_));
	  problem_ = boost::make_shared<problem_t>(*(factory_->getJerk()));
	}
      else
	{
	  pb = boost::make_shared<problem_t>(problem_->function());
	  problem_ = pb;
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
    const typename ProblemOverSplinesFactory<T, S>::numericLinearConstraintPtr_t
    ProblemOverSplinesFactory<T, S>::localConstraint (value_type time,
						      int order,
						      value_type value,
                                                      unsigned long spline_idx)
    {
      matrix_t A (1, static_cast<typename matrix_t::Index> (inputsize_));
      vector_t B (1);
      A.setZero();
      B.coeffRef(0) = -value;
      int splineOffset = 0;

      const splinePtr_t& spline = splines_[spline_idx];

      for (unsigned long i = 0; i < spline_idx; ++i)
        splineOffset += static_cast<int> (splines_[i]->getNumberControlPoints ());

      for (int j = 0; j < spline->dimension() + 1; ++j)
	{
	  int col = splineOffset + static_cast<int> (spline->interval(time)) - j;
	  unsigned long k = static_cast<unsigned long> (spline->interval(time) - j);
	  A.coeffRef (0, col) = spline->
	    basisPolynomials()[k][static_cast<unsigned long>(j)].derivative(time, order);
	}

      std::string name;
      if (order == 0)
        name = (boost::format ("f(%.3f) = %.3f (%s)")
		% time % value % spline->getName ()).str ();
      else
        name = (boost::format ("f^(%i) (%.3f) = %.3f (%s)")
		% order % time % value % spline->getName ()).str ();
      return boost::make_shared<numericLinearConstraint_t> (A, B, name);
    }
  }
}

#endif //! ROBOPTIM_TRAJECTORY_PROBLEM_OVER_SPLINES_FACTORY_HXX
