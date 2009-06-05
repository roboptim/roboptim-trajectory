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

#ifndef ROBOPTIM_TRAJECTORY_SPLINE_LENGTH_HH
# define ROBOPTIM_TRAJECTORY_SPLINE_LENGTH_HH
# include <roboptim/trajectory/fwd.hh>
# include <roboptim/trajectory/trajectory-cost.hh>

namespace roboptim
{
  /// Approximate the length of a Spline.
  ///
  /// The length is computed using:
  ///
  /// \f[\frac{1}{2} \int_{tmin}^{tmax} ||\ddot{\Gamma_p(t)}||^2 dt\f]
  ///
  /// \f$tmin\f$ and \f$tmax\f$ are given when the object is instantiated.
  class SplineLength : public TrajectoryCost<Spline>
  {
  public:
    /// Construct the function from a Spline and a definition interval.
    ///
    /// The interval allows to only compute the length of the Spline
    /// on a specific interval. The step associated to the interval controls
    /// the approximation precision.
    /// \param spline spline used for length computation
    /// \param interval discrete interval on which the length is computed.
    SplineLength (const Spline& spline, discreteInterval_t interval) throw ();

    virtual ~SplineLength () throw ();

  protected:
    void impl_compute (result_t& res, const argument_t& p) const throw ();
    void impl_gradient (gradient_t& grad, const argument_t& p, int i) const throw ();

  private:
    /// Interval on which the length is computed.
    discreteInterval_t interval_;
  };

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HH
