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
# include <roboptim/trajectory/sys.hh>

# include <boost/format.hpp>

# include <roboptim/core/visualization/gnuplot.hh>
# include <roboptim/core/visualization/gnuplot-commands.hh>
# include <roboptim/trajectory/trajectory.hh>
# include <roboptim/trajectory/stable-time-point.hh>

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

      template <typename T>
      Command plot_xytheta (const T& traj,
			    typename T::value_type step = .01);

      namespace detail
      {
	template <typename T>
	struct PlotTrajectory
	{
	  PlotTrajectory (const T& traj, std::string& str)
	    : traj_ (traj),
	      str_ (str)
	  {}

	  void operator () (const typename T::value_type) const
	  {}

	private:
	  const T& traj_;
	  std::string& str_;
	};
      } // end of namespace detail.

      template <unsigned N>
      Command plot_xy (const Trajectory<N>& traj,
		       typename Trajectory<N>::value_type step)
      {
	using boost::format;
	using namespace detail;
	assert (traj.outputSize () >= 2);
	Function::value_type min = Function::getLowerBound (traj.timeRange ());
        Function::value_type max = Function::getUpperBound (traj.timeRange ());

	if (min + step > max)
	  throw std::string ("bad interval");

        std::string str = (boost::format ("plot '-' title '%1%' with line\n")
			   % traj.getName ()).str ();

	for (double i = step; i < 1. - step; i += step)
	  {
	    Function::vector_t res = traj (i * tMax);
            str += (format ("%2.8f %2.8f\n")
		    % normalize (res[0])
                    % normalize (res[1])).str ();
	  }

	str += "e\n";
	return Command (str);
      }


      template <typename T>
      Command plot_xytheta (const T& traj, typename T::value_type step)
      {
	using boost::format;
	assert (traj.outputSize () == 3);
	Function::value_type min = Function::getLowerBound (traj.timeRange ());
	Function::value_type max = Function::getUpperBound (traj.timeRange ());
	Function::discreteInterval_t interval (min, max, step);

	if (min + step > max)
	  throw std::string ("bad interval");

	const char* color[] =
	  {
	    "red",
	    "blue",
	    "green"
	  };

	std::string str =
	  (boost::format ("plot '-' title '%1% (0)' with lines lc rgb '%2%'")
	   % traj.getName ()
	   % color[0]).str ();

	for (unsigned i = 1; i < traj.outputSize (); ++i)
	  {
	    str += (format (", '-' title '%1% (%2%)' with lines lc rgb '%3%'")
		    % traj.getName ()
		    % i
		    % color[i]).str ();
	  }
	str += "\n";
	for (unsigned component = 0; component < traj.outputSize (); ++component)
	  {
	    for (double i = step; i < 1. - step; i += step)
	      {
		StableTimePoint timePoint = i * tMax;
		Function::vector_t res = traj (timePoint);
                str += (format ("%2.8f %2.8f\n")
			% normalize (timePoint.getTime (traj.timeRange ()))
			% normalize (res [component])).str ();
	      }
	    str += "e\n";
	  }
	return Command (str);
      }

    } // end of namespace gnuplot.
  } // end of namespace visualization.
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_VISUALIZATION_TRAJECTORY_HH
