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

#ifndef ROBOPTIM_TRAJECTORY_SPLINE_HH
# define ROBOPTIM_TRAJECTORY_SPLINE_HH
# include <roboptim/trajectory/sys.hh>

# include <roboptim/trajectory/trajectory.hh>
# include <roboptim/trajectory/fwd.hh>

/// \internal
class bspline;

namespace roboptim
{
  /// \addtogroup roboptim_function
  /// @{

  /// \brief Spline trajectory.
  ///
  /// Implement B-Spline.
  class Spline : public Trajectory<4>
  {
  public:
    /// \brief Instantiate a Spline from its definition.
    ///
    /// Instantiate a Spline from its triplet definition:
    /// - time range (\f$t \in [f_min, t_max]\f$): interval on
    ///   which the spline is defined,
    /// - dimension: output space dimension,
    /// - parameters: vector containing concatenated control points.
    ///
    /// The vector parameters has to satisfy the following conditions:
    /// - size should be at least 4,
    /// - there should be at least two control points,
    /// - the vector has to be valid (it has to contain complete control
    ///   points so its size has to be a multiple of the dimension size).
    ///
    /// \param timeRange spline time range
    /// \param dimension spline dimension
    /// \param parameters vector of parameters defining control points
    /// \param name function title
    Spline (interval_t timeRange, size_type dimension,
	    const vector_t& parameters, std::string name = "spline") throw ();

    /// \brief Copy constructor.
    /// \param spline spline that will be copied
    Spline (const Spline& spline) throw ();

    virtual ~Spline () throw ();

    /// \brief Modify spline parameters.
    virtual void setParameters (const vector_t&) throw ();

    virtual jacobian_t variationConfigWrtParam (double t) const throw ();
    virtual jacobian_t variationDerivWrtParam (double t, size_type order)
      const throw ();
    virtual value_type singularPointAtRank (size_type rank) const;
    virtual vector_t derivBeforeSingularPoint (size_type rank, size_type order)
      const;
    virtual vector_t derivAfterSingularPoint (size_type rank, size_type order)
      const;

    ROBOPTIM_IMPLEMENT_CLONE(Spline)

    virtual Trajectory<derivabilityOrder>* resize (interval_t timeRange)
      const throw ()
    {
      return new Spline (timeRange, this->outputSize (), this->parameters ());
    }

    /// \brief Display the function on the specified output stream.
    ///
    /// \param o output stream used for display
    /// \return output stream
    virtual std::ostream& print (std::ostream& o) const throw ();

  protected:
    void impl_compute (result_t&, double) const throw ();
    void impl_derivative (gradient_t& g, double x, size_type order)
      const throw ();

  private:
    /// \brief Convert parameters to internal representation.
    /// \return internal parameter representation
    vector_t makeBSplineVector ();

    /// \brief Number of control points.
    int nbp_;
    /// \brief Pointer to internal spline implementation.
    bspline* spline_;
  };

  /// Example shows Spline use.
  /// \example spline-gradient.cc

  /// @}

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_TRAJECTORY_HH
