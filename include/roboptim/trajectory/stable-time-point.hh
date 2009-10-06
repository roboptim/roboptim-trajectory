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

#ifndef ROBOPTIM_TRAJECTORY_STABLE_TIME_POINT_HH
# define ROBOPTIM_TRAJECTORY_STABLE_TIME_POINT_HH
# include <roboptim/trajectory/sys.hh>

# include <roboptim/core/function.hh>

namespace roboptim
{
  class TMax {};
  static const TMax tMax = TMax ();

  class StableTimePoint
  {
  public:
    typedef Function::value_type value_type;

    explicit StableTimePoint (value_type alpha)
      : alpha_ (alpha)
    {
      assert (alpha_ >= 0. && alpha <= 1.);
    }

    const value_type& getAlpha () const
    {
      return alpha_;
    }

    value_type getTime (Function::interval_t timeRange) const
    {
      value_type tmin = Function::getLowerBound (timeRange);
      value_type tmax = Function::getUpperBound (timeRange);

      value_type res = tmin + alpha_ * (tmax - tmin);
      res = (res > tmax) ? tmax : ((res < tmin) ? tmin : res);
      return res;
    }
  private:
    value_type alpha_;
  };

  inline StableTimePoint operator* (Function::value_type alpha, TMax)
  {
    return StableTimePoint (alpha);
  }

  inline StableTimePoint operator* (TMax, Function::value_type alpha)
  {
    return StableTimePoint (alpha);
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_STABLE_TIME_POINT_HH
