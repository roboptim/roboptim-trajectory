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

namespace roboptim
{
  namespace trajectory
  {
    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::addSpline (S& spline)
    {
      if (static_cast<unsigned long>(spline.dimension()) == dimension_
          && spline.timeRange().first == t0_
          && spline.timeRange().second == tmax_)
	{
	  splines_.push_back(boost::make_shared<S>(spline));
	  inputsize_ += spline.getNumberControlPoints();
	}
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::updateRange (interval_t newRange, bool buildCostFunction)
    {
      assert(newRange.first >= t0_);
      assert(newRange.first < newRange.second);
      t0 (newRange.first);
      tmax (newRange.second);
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
    void ProblemOverSplinesFactory<T, S>::addConstraint(value_type time, int order, std::vector<value_type> value, scaling_t scaling)
    {
      assert (problem_->function().inputSize() == inputsize_);
      constraints_.resize(constraints_.size() + 1);
      constraints_.back().first = time;
      for (unsigned long i = 0; i < splines_.size(); ++i)
	{
	  constraints_.back().second.push_back(static_cast<freeze_t>(boost::make_tuple(localConstraint(time, order, value[i], i), S::makeInterval(0, 0), scaling[i])));
	  if (time <= tmax_ && time >= t0_)
	    problem_->addConstraint(boost::get<0>(boost::get<freeze_t>(constraints_.back().second.back())), S::makeInterval(0, 0), scaling[i]);
	}
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::addConstraint(value_type startingPoint, int order, intervals_t range, scalingVect_t scaling)
    {
      assert (problem_->function().inputSize() == inputsize_);
      constraints_.resize(constraints_.size() + 1);
      constraints_.back().first = startingPoint;
      for (unsigned long i = 0; i < splines_.size(); ++i)
	{
	  range_[0] = S::makeLowerInterval(range[i].first);
	  range_[1] = S::makeUpperInterval(range[i].second);
	  constraints_.back().second.push_back(boost::make_tuple(boost::make_shared<ConstraintsOverSplines<T, S> >(splines_, i, order, startingPoint, inputsize_), range_, scaling[i]));
	  if (startingPoint < tmax_ && startingPoint >= t0_)
	    problem_->addConstraint(boost::get<0>(boost::get<constraint_t>(constraints_.back().second.back())), range_, scaling[i]);
	}
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::addConstraint(value_type time, int order, std::vector<value_type> value)
    {
      scaling_t scaling (value.size());
      for (unsigned long i = 0; i < scaling.size(); ++i)
        scaling[i] = 1;
      addConstraint(time, order, value, scaling);
    }

    template <typename T, typename S>
    void ProblemOverSplinesFactory<T, S>::addConstraint(value_type startingPoint, int order, intervals_t range)
    {
      scalingVect_t scaling (range.size());
      for (unsigned long i = 0; i < scaling.size(); ++i)
	{
	  scaling[i] = scaling_t(2);
	  scaling[i][0] = 1;
	  scaling[i][1] = 1;
	}
      addConstraint(startingPoint, order, range, scaling);
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
      constraint_t c;
      freeze_t f;
      for (unsigned long i = 0; i < constraints_.size(); ++i)
        if (constraints_[i].first >= t0_ && constraints_[i].first <= tmax_)
          for (unsigned long j = 0; j < constraints_[i].second.size(); ++j)
	    {
	      switch (constraints_[i].second[j].which())
		{
		case 0 :
		  c = boost::get<constraint_t>(constraints_[i].second[j]);
		  problem_->addConstraint(boost::get<0>(c), boost::get<1>(c), boost::get<2>(c));
		  break;
		case 1 :
		  f = boost::get<freeze_t>(constraints_[i].second[j]);
		  problem_->addConstraint(boost::get<0>(f), boost::get<1>(f), boost::get<2>(f));
		  break;
		}
	    }
    }

    template <typename T, typename S>
    const boost::shared_ptr<GenericNumericLinearFunction<T> > ProblemOverSplinesFactory<T, S>::localConstraint(value_type time, int order, value_type value, unsigned long spline)
    {
      matrix_t A (1, static_cast<int>(inputsize_));
      vector_t B (1);
      A.setZero();
      B.coeffRef(0) = -value;
      int splineOffset = 0;
      for (unsigned long i = 0; i < spline; ++i)
        splineOffset += static_cast<int>(splines_[i]->getNumberControlPoints());
      for (int j = 0; j < splines_[spline]->dimension() + 1; ++j)
	A.coeffRef(0, splineOffset + static_cast<int>(splines_[spline]->interval(time)) - j) =
	  splines_[spline]->basisPolynomials()[static_cast<unsigned long>(splines_[spline]->interval(time) - j)][static_cast<unsigned long>(j)].derivative(time, order);
      return boost::make_shared<GenericNumericLinearFunction<T> >(A, B);
    }
  }
}

#endif //! ROBOPTIM_TRAJECTORY_PROBLEM_OVER_SPLINES_FACTORY_HXX
