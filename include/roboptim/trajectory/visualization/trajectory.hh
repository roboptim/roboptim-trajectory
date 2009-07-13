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

#ifndef ROBOPTIM_TRAJECTORY_VISUALIZATION_TRAJECTORY_HH
# define ROBOPTIM_TRAJECTORY_VISUALIZATION_TRAJECTORY_HH
# include <boost/format.hpp>

# include <roboptim/core/visualization/gnuplot-commands.hh>
# include <roboptim/trajectory/trajectory.hh>

namespace roboptim
{
  namespace visualization
  {
    namespace gnuplot
    {
      /// \brief Plot a 2D trajectory with Gnuplot.
      ///
      /// Plot a 2D trajectory in Gnuplot.
      /// \param traj trajectory to be displayed
      /// \param step discretization step.
      /// \return Gnuplot command
      template <unsigned N>
      Command plot_xy (const Trajectory<N>& traj,
		       typename Trajectory<N>::value_type step = .01);

      namespace detail
      {
	template <typename T>
	struct PlotTrajectory
	{
	  PlotTrajectory (const T& traj, std::string& str)
	    : traj_ (traj),
	      str_ (str)
	  {}

	  void operator () (const typename T::value_type t) const
	  {
	    Function::vector_t res = traj_ (t);
	    assert (res.size () >= 2);
	    str_ += (boost::format ("%1% %2%\n") % res[0] % res [1]).str ();
	  }

	private:
	  const T& traj_;
	  std::string& str_;
	};
      } // end of namespace detail.

      template <unsigned N>
      Command plot_xy (const Trajectory<N>& traj,
		       typename Trajectory<N>::value_type step)
      {
	using namespace detail;
	assert (traj.outputSize () >= 2);
	Function::value_type min = Function::getLowerBound (traj.timeRange ());
	Function::value_type max = Function::getUpperBound (traj.timeRange ());
	Function::discreteInterval_t interval (min, max, step);

	std::string str = "plot '-' with line\n";
	traj.foreach (interval, PlotTrajectory<Trajectory<N> > (traj, str));
	str += "e\n";
	return Command (str);
      }
    } // end of namespace gnuplot.
  } // end of namespace visualization.
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_VISUALIZATION_TRAJECTORY_HH
