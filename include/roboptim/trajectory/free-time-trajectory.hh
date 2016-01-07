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
# include <roboptim/trajectory/sys.hh>

# include <boost/static_assert.hpp>
# include <boost/type_traits/is_base_of.hpp>

# include <roboptim/trajectory/trajectory.hh>

namespace roboptim
{
namespace trajectory
{
  /// \addtogroup roboptim_function
  /// @{

  /// \brief Decorate a trajectory to make time scalable.
  ///
  /// Build a trajectory from an input trajectory and a time scale
  /// factor.
  template <typename T>
  class FreeTimeTrajectory : public Trajectory<T::derivabilityOrder>
  {
    /// Check that T is a trajectory type.
    BOOST_STATIC_ASSERT((boost::is_base_of
			 <Trajectory<T::derivabilityOrder>, T>::value));

  public:
    /// \brief Parent type and imports.
    ROBOPTIM_NTIMES_DERIVABLE_FUNCTION_FWD_TYPEDEFS_
      (Trajectory<T::derivabilityOrder>);
    /// \brief Fixed point trajectory type.
    typedef T fixedTimeTrajectory_t;
    /// \brief Self type.
    typedef FreeTimeTrajectory<T> self_t;


    /// \brief Import interval type.
    typedef typename parent_t::interval_t interval_t;

    using parent_t::normalizeAngles;
    using parent_t::variationConfigWrtParam;
    using parent_t::variationDerivWrtParam;

    /// Constructor with fixed definition interval trajectory
    ///
    /// \param traj trajectory defining this one by reparameterization
    /// \param s time scale
    FreeTimeTrajectory (const fixedTimeTrajectory_t& traj, value_type s);

    FreeTimeTrajectory (const self_t& traj);

    virtual ~FreeTimeTrajectory ();


    virtual jacobian_t variationConfigWrtParam (double t) const;
    virtual jacobian_t variationDerivWrtParam (double t, size_type order)
      const;

    virtual value_type singularPointAtRank (size_type rank) const;
    virtual vector_t derivBeforeSingularPoint (size_type rank, size_type order)
      const;
    virtual vector_t derivAfterSingularPoint (size_type rank, size_type order)
      const;

    virtual void setParameters (const_vector_ref);

    /// \brief Get time scale factor.
    /// \return time scale factor.
    value_type timeScale () const;

    size_type getTimeScalingIndex () const
    {
      return 0;
    }

    ROBOPTIM_IMPLEMENT_CLONE (self_t)

    /// \brief Display the function on the specified output stream.
    ///
    /// \param o output stream used for display
    /// \return output stream
    virtual std::ostream& print (std::ostream& o) const;

    /// \brief Normalize angles in parameters array.
    ///
    /// Make sure angles are continuous.
    /// \param index Angles index in parameter array.
    virtual void normalizeAngles (size_type index);

    const fixedTimeTrajectory_t&
    getFixedTimeTrajectory () const
    {
      assert (trajectory_);
      return *trajectory_;
    }

    self_t*
    resize (interval_t) const
    {
      assert (trajectory_);
      value_type tMin = this->getLowerBound (this->timeRange ());
      value_type tMax = this->getUpperBound (this->timeRange ());
      value_type tmin = this->getLowerBound (this->trajectory_->timeRange ());
      value_type tmax = this->getUpperBound (this->trajectory_->timeRange ());
      assert (tMax != tMin);

      value_type scale = (tmax - tmin) / (tMax - tMin);

      assert (this->scaleTime (tMin) == tmin);
      assert (this->scaleTime (tMax) == tmax);

      self_t* res =
	new self_t (*trajectory_, scale);
      assert (res->timeRange () == this->timeRange ());
      return res;
    }

    fixedTimeTrajectory_t*
    makeFixedTimeTrajectory () const
    {
      assert (trajectory_);
      return trajectory_->resize (this->timeRange ());
    }

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
    double scaleTime (double t) const;
    double unscaleTime (double t) const;

    jacobian_t variationConfigWrtParam (StableTimePoint tp) const;
    jacobian_t variationDerivWrtParam (StableTimePoint tp, size_type order)
      const;

  protected:
    void impl_compute (result_ref, double) const;
    void impl_derivative (derivative_ref g, double x, size_type order)
      const;
    void impl_derivative (derivative_ref g, StableTimePoint, size_type order)
      const;
  private:
    /// \brief Input fixed time trajectory.
    fixedTimeTrajectory_t* trajectory_;
  };

  /// Example shows FreeTimeTrajectory use.
  /// \example spline-time-optimization.cc

  /// @}

  Function::vector_t
  addScaleToParameters (Function::const_vector_ref p,
			Function::value_type t = 1.);
  Function::vector_t
  removeScaleFromParameters (Function::const_vector_ref v);

  inline Function::vector_t
  addScaleToParameters (Function::const_vector_ref p,
			Function::value_type t)
  {

    Function::vector_t res (p.size () + 1);
    res[0] = t;
    res.segment( 1, p.size () ) = p;
    return res;
  }

  inline Function::vector_t
  removeScaleFromParameters (Function::const_vector_ref p)
  {
    Function::vector_t res (p.size () - 1);
    res = p.segment( 1, p.size () - 1);
    return res;
  }

} // end of namespace trajectory.
} // end of namespace roboptim.

# include <roboptim/trajectory/free-time-trajectory.hxx>
#endif //! ROBOPTIM_TRAJECTORY_FREETIMETRAJECTORY_HH
