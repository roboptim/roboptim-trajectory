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

#ifndef ROBOPTIM_TRAJECTORY_JERK_OVER_SPLINES_FACTORY_HXX
# define ROBOPTIM_TRAJECTORY_JERK_OVER_SPLINES_FACTORY_HXX

# include <boost/make_shared.hpp>

namespace roboptim
{
  namespace trajectory
  {
    template <typename S, typename T>
    JerkOverSplinesFactory<S, T>::JerkOverSplinesFactory
    (const std::vector<splinePtr_t>& splines, const interval_t& range,
     value_type scaling)
      : splines_ (splines),
        range_ (range),
        inputSize_ (static_cast<size_type> (splines.size ())
                    * static_cast<size_type> (splines[0]->getNumberControlPoints ())),
        A_ (inputSize_, inputSize_),
        B_ (inputSize_),
        scaling_ (scaling > 0 ? scaling : chooseScaling ()),
        f_ ()
    {
      assert (splines[0]->order () == 3);
      ROBOPTIM_DEBUG_ONLY
        (typename spline_t::size_type nbp_ = splines[0]->getNumberControlPoints();
         for (size_t i = 1; i < splines.size (); ++i)
           assert (nbp_ == splines[i]->getNumberControlPoints());
        );

      A_.setZero ();
      B_.setZero ();

      assert (range_.first < range_.second);
      size_t start = static_cast<size_t> (splines[0]->interval (range_.first));
      size_t end = static_cast<size_t> (splines[0]->interval (range_.second));

      size_t N = static_cast<size_t> (splines[0]->order ());
      size_t nbParameters = static_cast<size_t> (splines[0]->getNumberControlPoints());

      for (size_t n = 0; n < splines.size (); ++n)
        for (size_t i = start; i < end + 1; ++i)
          for (size_t j = 0; j < N+1; ++j)
            for (size_t k = 0; k < N+1; ++k)
            {
              size_type ii = static_cast<size_type> (i);
              A_.coeffRef (static_cast<size_type> (i + n*nbParameters - j),
                           static_cast<size_type> (i + n*nbParameters - k))
                += splines_[n]->basisPolynomials ()[i-j][j].derivative (0,3) *
                   splines_[n]->basisPolynomials ()[i-k][k].derivative (0,3) *
                   (std::min (splines_[n]->knotVector ()[ii+1], range_.second)
                    - std::max (splines_[n]->knotVector ()[ii], range_.first));
            }

      A_ *= scaling_;
      f_ = boost::make_shared<numericQuadraticFunction_t> (A_, B_, "Jerk");
    }

    template <typename S, typename T>
    void JerkOverSplinesFactory<S, T>::updateRange (const interval_t& range)
    {
      assert (range.first < range.second);
      range_ = range;

      size_t newStart = static_cast<size_t> (splines_[0]->interval (range_.first));
      size_t newEnd = static_cast<size_t> (splines_[0]->interval (range_.second));
      size_t N = static_cast<size_t> (splines_[0]->order ());
      size_t nbParameters = static_cast<size_t> (splines_[0]->getNumberControlPoints ());
      A_.setZero ();

      for (size_t n = 0; n < splines_.size (); ++n)
        for (size_t i = newStart; i < newEnd + 1; ++i)
          for (size_t j = 0; j < N+1; ++j)
            for (size_t k = 0; k < N+1; ++k)
            {
              size_type ii = static_cast<size_type> (i);
              A_.coeffRef (static_cast<size_type> (i + n*nbParameters - j),
                           static_cast<size_type> (i + n*nbParameters - k))
                += splines_[n]->basisPolynomials ()[i-j][j].derivative (0,3) *
                   splines_[n]->basisPolynomials ()[i-k][k].derivative (0,3) *
                   (std::min (splines_[n]->knotVector ()[ii+1], range_.second)
                    - std::max (splines_[n]->knotVector ()[ii], range_.first));
            }

      A_ *= scaling_;
      f_ = boost::make_shared<numericQuadraticFunction_t> (A_, B_, "Jerk");
    }

    template <typename S, typename T>
    typename JerkOverSplinesFactory<S, T>::numericQuadraticFunctionPtr_t
    JerkOverSplinesFactory<S, T>::getJerk ()
    {
      return f_;
    }

    template <typename S, typename T>
    typename JerkOverSplinesFactory<S, T>::value_type
    JerkOverSplinesFactory<S, T>::chooseScaling () const
    {
      value_type length = range_.second - range_.first;
      size_t n = splines_.size ();
      return 1. / (static_cast<value_type> (n) * length * 1e7);
    }
  }
}

#endif //! ROBOPTIM_TRAJECTORY_JERK_OVER_SPLINES_FACTORY_HXX
