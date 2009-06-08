// Copyright (C) 2009 by Florent Lamiraux, Thomas Moulard, AIST, CNRS, INRIA.
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

#ifndef ROBOPTIM_TRAJECTORY_FREETIMETRAJECTORY_HH
# define ROBOPTIM_TRAJECTORY_FREETIMETRAJECTORY_HH
# include <roboptim/trajectory/trajectory.hh>

namespace roboptim
{
  /// \addtogroup roboptim_function
  /// @{

  /// \brief Decorate a trajectory to make time scalable.
  ///
  /// Build a trajectory from an input trajectory and a time scale
  /// factor.
  template <unsigned DerivabilityOrder>
  class FreeTimeTrajectory : public Trajectory<DerivabilityOrder>
  {
  public:
    /// Constructor with fixed definition interval trajectory
    ///
    /// \param traj trajectory defining this one by reparameterization
    /// \param s time scale
    FreeTimeTrajectory (const Trajectory<DerivabilityOrder>& traj, double s) throw ();

    virtual ~FreeTimeTrajectory () throw ();

    virtual vector_t operator () (double) const throw ();

    /// \brief Compute the derivative of the function.
    ///
    /// Derivative is computed for a certain order, at a given point.
    /// \param x point at which the derivative will be computed
    /// \param order derivative order (if 0 then function is evaluated)
    /// \return derivative vector
    virtual vector_t derivative (double x, size_type order) const throw ();

    virtual jacobian_t variationConfigWrtParam (double t) const throw ();
    virtual jacobian_t variationDerivWrtParam (double t, size_type order)
      const throw ();
    virtual value_type singularPointAtRank (size_type rank) const;
    virtual vector_t derivBeforeSingularPoint (size_type rank, size_type order)
      const;
    virtual vector_t derivAfterSingularPoint (size_type rank, size_type order)
      const;

    /// \brief Get time scale factor.
    /// \return time scale factor.
    double timeScale () const throw ();

    /// \brief Display the function on the specified output stream.
    ///
    /// \param o output stream used for display
    /// \return output stream
    virtual std::ostream& print (std::ostream& o) const throw ();
  private:
    /// \brief Scale input time argument.
    ///
    /// Scale input argument with the same factor that the input
    /// trajectory:
    /// \f[t' = t_{min} + \lambda * (t - t_{min})\f]
    /// where \f$[t_{min}, t_{max}]\f$ is the input trajectory time range and
    /// \f[\lambda\f] the scale factor.
    ///
    /// \param t input time
    /// \return new scaled time
    double scaleTime (double t) const throw ();

    /// \brief Input fixed time trajectory.
    const Trajectory<DerivabilityOrder>& trajectory_;
    /// \brief Store time scaling parameter \f$\textbf{p}_{m+1}\f$.
    double timeScale_;
  };

  /// @}

} // end of namespace roboptim.

# include <roboptim/trajectory/free-time-trajectory.hxx>
#endif //! ROBOPTIM_TRAJECTORY_FREETIMETRAJECTORY_HH
