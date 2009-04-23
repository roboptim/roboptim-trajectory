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

/**
 * \brief Class Spline declaration.
 */

#ifndef ROBOPTIM_TRAJECTORY_SPLINE_HH
# define ROBOPTIM_TRAJECTORY_SPLINE_HH
# include <roboptim-trajectory/trajectory.hh>

# include <roboptim-trajectory/fwd.hh>

class bspline;

namespace roboptim
{
  class Spline : public Trajectory<4>
  {
  public:
    Spline (size_type, const vector_t&, int, int, int) throw ();
    virtual ~Spline () throw ();

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

  private:
    bspline* spline_;
  };
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HH
