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

    virtual vector_t derivative (double x, size_type order) const throw ();

    virtual jacobian_t variationConfigWrtParam (double t) const throw ();
    virtual jacobian_t variationDerivWrtParam (double t, size_type order)
      const throw ();
    virtual value_type singularPointAtRank (size_type rank) const;
    virtual vector_t derivBeforeSingularPoint (size_type rank, size_type order)
      const;
    virtual vector_t derivAfterSingularPoint (size_type rank, size_type order)
      const;

    /// Get index of time scaling parameter in vector
    double timeScale () const throw ();

    virtual std::ostream& print (std::ostream&) const throw ();
  private:
    double scaleTime (double) const throw ();

    /// Fixed time trajectory
    const Trajectory<DerivabilityOrder>& trajectory_;
    /// Store time scaling parameter \f$\textbf{p}_{m+1}\f$.
    double timeScale_;
  };
} // end of namespace roboptim.

# include <roboptim/trajectory/free-time-trajectory.hxx>
#endif //! ROBOPTIM_TRAJECTORY_FREETIMETRAJECTORY_HH
