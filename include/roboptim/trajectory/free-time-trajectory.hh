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
    /// \brief Parent type.
    typedef Trajectory<T::derivabilityOrder> parent_t;
    /// \brief Fixed point trajectory type.
    typedef T fixedTimeTrajectory_t;
    /// \brief Self type.
    typedef FreeTimeTrajectory<T> self_t;


    /// \brief Import value type.
    typedef typename parent_t::value_type value_type;
    /// \brief Import size type.
    typedef typename parent_t::size_type size_type;
    /// \brief Import result type.
    typedef typename parent_t::result_t result_t;
    /// \brief Import gradient type.
    typedef typename parent_t::gradient_t gradient_t;
    /// \brief Import vector type.
    typedef typename parent_t::vector_t vector_t;
    /// \brief Import jacobian type.
    typedef typename parent_t::jacobian_t jacobian_t;
    /// \brief Import interval type.
    typedef typename parent_t::interval_t interval_t;

    using parent_t::normalizeAngles;
    using parent_t::variationConfigWrtParam;
    using parent_t::variationDerivWrtParam;

    /// Constructor with fixed definition interval trajectory
    ///
    /// \param traj trajectory defining this one by reparameterization
    /// \param s time scale
    FreeTimeTrajectory (const fixedTimeTrajectory_t& traj, value_type s)
      throw ();

    FreeTimeTrajectory (const self_t& traj) throw ();

    virtual ~FreeTimeTrajectory () throw ();


    virtual jacobian_t variationConfigWrtParam (double t) const throw ();
    virtual jacobian_t variationDerivWrtParam (double t, size_type order)
      const throw ();

    virtual value_type singularPointAtRank (size_type rank) const;
    virtual vector_t derivBeforeSingularPoint (size_type rank, size_type order)
      const;
    virtual vector_t derivAfterSingularPoint (size_type rank, size_type order)
      const;

    virtual void setParameters (const vector_t&) throw ();

    /// \brief Get time scale factor.
    /// \return time scale factor.
    value_type timeScale () const throw ();

    size_type getTimeScalingIndex () const throw ()
    {
      return 0;
    }

    ROBOPTIM_IMPLEMENT_CLONE (self_t)

    /// \brief Display the function on the specified output stream.
    ///
    /// \param o output stream used for display
    /// \return output stream
    virtual std::ostream& print (std::ostream& o) const throw ();

    /// \brief Normalize angles in parameters array.
    ///
    /// Make sure angles are continuous.
    /// \param index Angles index in parameter array.
    virtual void normalizeAngles (size_type index) throw ();

    const fixedTimeTrajectory_t&
    getFixedTimeTrajectory () const throw ()
    {
      assert (trajectory_);
      return *trajectory_;
    }

    self_t*
    resize (interval_t) const throw ()
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
    makeFixedTimeTrajectory () const throw ()
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
    double scaleTime (double t) const throw ();
    double unscaleTime (double t) const throw ();

    jacobian_t variationConfigWrtParam (StableTimePoint tp) const throw ();
    jacobian_t variationDerivWrtParam (StableTimePoint tp, size_type order)
      const throw ();

  protected:
    void impl_compute (result_t&, double) const throw ();
    void impl_derivative (gradient_t& g, double x, size_type order)
      const throw ();
    void impl_derivative (gradient_t& g, StableTimePoint, size_type order)
      const throw ();
  private:
    /// \brief Input fixed time trajectory.
    fixedTimeTrajectory_t* trajectory_;
  };

  /// Example shows FreeTimeTrajectory use.
  /// \example spline-time-optimization.cc

  /// @}

  Function::vector_t
  addScaleToParameters (const Function::vector_t& p,
			Function::value_type t = 1.);
  Function::vector_t
  removeScaleFromParameters (const Function::vector_t& v);

  inline Function::vector_t
  addScaleToParameters (const Function::vector_t& p,
			Function::value_type t)
  {

    Function::vector_t res (p.size () + 1);
    res[0] = t;
    res.segment( 1, p.size () ) = p;
    return res;
  }

  inline Function::vector_t
  removeScaleFromParameters (const Function::vector_t& p)
  {
    Function::vector_t res (p.size () - 1);
    res = p.segment( 1, p.size () - 1);
    return res;
  }

} // end of namespace roboptim.

# include <roboptim/trajectory/free-time-trajectory.hxx>
#endif //! ROBOPTIM_TRAJECTORY_FREETIMETRAJECTORY_HH
