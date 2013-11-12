// Copyright (C) 2013 by Thomas Moulard, AIST, CNRS, INRIA.
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

#ifndef ROBOPTIM_TRAJECTORY_FILTER_VECTOR_INTERPOLATION_HH
# define ROBOPTIM_TRAJECTORY_FILTER_VECTOR_INTERPOLATION_HH
# include <vector>
# include <stdexcept>
# include <boost/make_shared.hpp>
# include <boost/shared_ptr.hpp>
# include <roboptim/trajectory/trajectory.hh>

namespace roboptim
{
  /// \brief Takes a vector or argument and differentiate it.
  ///
  /// It is common to have an optimization vector where a pattern is
  /// repeated:
  ///
  /// [x_0^0 ... x_N^0 ... x_0^M x_N^M]
  ///
  /// The pattern is here of length N and repeated M times.
  ///
  /// By using numerical interpolation, it is possible to
  /// differentiate this pattern.
  ///
  /// A classical example are trajectories.
  class VectorInterpolation : public Trajectory<3>
  {
  public:
    ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS
    (Trajectory<3>);

    ROBOPTIM_IMPLEMENT_CLONE (VectorInterpolation);

    typedef boost::shared_ptr<VectorInterpolation> VectorInterpolationShPtr_t;

    explicit VectorInterpolation
    (const vector_t& x, size_type outputSize, value_type dt)
      throw (std::runtime_error);
    ~VectorInterpolation () throw ();

    /// \brief Store parameters and update coefficients.
    void setParameters (const vector_t&) throw (std::runtime_error);

    jacobian_t variationConfigWrtParam (double t) const throw ();
    jacobian_t variationDerivWrtParam (double t, size_type order)
      const throw ();
    value_type singularPointAtRank (size_type rank) const;
    vector_t derivBeforeSingularPoint (size_type rank, size_type order) const;
    vector_t derivAfterSingularPoint (size_type rank, size_type order) const;

    jacobian_t variationConfigWrtParam (StableTimePoint tp)
      const throw ();
    jacobian_t
    variationDerivWrtParam (StableTimePoint tp, size_type order)
      const throw ();

  protected:
    void impl_compute (result_t& result, double t) const throw ();
    void impl_derivative (gradient_t& derivative,
			  double argument,
			  size_type order = 1) const throw ();
    void impl_derivative (gradient_t& g, StableTimePoint, size_type order)
      const throw ();
    Trajectory<3>* resize (interval_t timeRange)
      const throw ();
  private:
    vector_t dx_;
    value_type dt_;
  };

  boost::shared_ptr<VectorInterpolation>
  vectorInterpolation (VectorInterpolation::vector_t x,
		       VectorInterpolation::size_type outputSize,
		       VectorInterpolation::value_type dt = 1.)
  {
    return boost::make_shared<VectorInterpolation> (x, outputSize, dt);
  }

} // end of namespace roboptim.

# include <roboptim/trajectory/vector-interpolation.hxx>
#endif //! ROBOPTIM_TRAJECTORY_FILTER_VECTOR_INTERPOLATION_HH
