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

#ifndef ROBOPTIM_TRAJECTORY_VISUALIZATION_SPEED_HH
# define ROBOPTIM_TRAJECTORY_VISUALIZATION_SPEED_HH
# include <roboptim/trajectory/sys.hh>

# include <boost/format.hpp>

# include <roboptim/core/visualization/gnuplot-commands.hh>
# include <roboptim/trajectory/trajectory.hh>
# include <roboptim/trajectory/stable-time-point.hh>
# include <roboptim/trajectory/frontal-speed.hh>
# include <roboptim/trajectory/orthogonal-speed.hh>

namespace roboptim
{
  namespace visualization
  {
    namespace gnuplot
    {
      template <typename T>
      Command plot_speeds (const T& traj,
			   typename T::value_type step = .01);

      template <typename T>
      Command plot_speeds (const T& traj,
			   typename T::value_type step)
      {
	using boost::format;
	assert (traj.outputSize () == 3);
	Function::value_type min = Function::getLowerBound (traj.timeRange ());
	Function::value_type max = Function::getUpperBound (traj.timeRange ());
	Function::discreteInterval_t interval (min, max, step);

	if (min + step > max)
	  throw std::string ("bad interval");

	std::string str =
	  (boost::format
	   ("plot '-' title '%1% (frontal)' with lines lc rgb 'blue', "
	    "'-' title '%1% (orthogonal)' with lines lc rgb 'violet', "
	    "'-' title '%1% (omega)' with lines lc rgb 'yellow'\n")
	   % traj.getName ()).str ();

	{
	  FrontalSpeed frontalSpeed;
	  for (double i = step; i < 1. - step; i += step)
	    {
	      StableTimePoint timePoint = i * tMax;
	      double t = timePoint.getTime (traj.timeRange ());

	      str += (format ("%1f %2f\n")
		      % t
		      % frontalSpeed (traj.state (t, 1))[0]).str ();
	    }
	  str += "e\n";
	}

	{
	  OrthogonalSpeed orthogonalSpeed;
	  for (double i = step; i < 1. - step; i += step)
	    {
	      StableTimePoint timePoint = i * tMax;
	      double t = timePoint.getTime (traj.timeRange ());

	      str += (format ("%1f %2f\n")
		      % t
		      % orthogonalSpeed (traj.state(t, 1))[0]).str ();
	    }
	  str += "e\n";
	}

	{
	  for (double i = step; i < 1. - step; i += step)
	    {
	      StableTimePoint timePoint = i * tMax;
	      double t = timePoint.getTime (traj.timeRange ());

	      str += (format ("%1f %2f\n")
		      % t
		      % traj.derivative (timePoint, 1)[2]).str ();
	    }
	  str += "e\n";
	}
	return Command (str);
      }


    } // end of namespace gnuplot.
  } // end of namespace visualization.
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_VISUALIZATION_SPEED_HH
