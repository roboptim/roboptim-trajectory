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
# include <roboptim/trajectory/sys.hh>

# include <boost/optional/optional_fwd.hpp>

# include <roboptim/trajectory/fwd.hh>
# include <roboptim/trajectory/cubic-b-spline.hh>
# include <roboptim/trajectory/trajectory-cost.hh>


namespace roboptim
{
  /// \addtogroup roboptim_function
  /// @{

  /// \brief Approximate the length of a Spline.
  ///
  /// The length is computed using:
  ///
  /// \f[\frac{1}{2} \int_{t_{min}}^{t_{max}} ||\ddot{\Gamma_p(t)}||^2 dt\f]
  ///
  /// \f$t_{min}\f$ and \f$t_{max}\f$ are given when the object is instantiated.
  class SplineLength : public TrajectoryCost<CubicBSpline>
  {
  public:
    /// Construct the function from a Spline and a definition interval.
    ///
    /// The interval allows to only compute the length of the Spline
    /// on a specific interval. The step associated to the interval controls
    /// the approximation precision.
    /// \param spline spline used for length computation
    /// \param interval interval on which the length is computed.
    /// \param nDiscretizationPoints number of discretization points
    SplineLength (const CubicBSpline& spline,
		  const size_type nDiscretizationPoints = 100,
		  boost::optional<interval_t> interval = boost::none_t ()) throw ();

    virtual ~SplineLength () throw ();

  protected:
    void impl_compute (result_t&, const argument_t&) const throw ();
    void impl_gradient (gradient_t&, const argument_t&, size_type) const throw ();

  private:
    /// \brief Interval on which the length is computed.
    interval_t interval_;
    /// \brief Number of discretization points.
    size_type nDiscretizationPoints_;
  };

  /// \@

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HH
