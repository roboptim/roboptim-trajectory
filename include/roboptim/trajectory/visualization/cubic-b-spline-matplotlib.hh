// Copyright (C) 2015 by Félix Darricau, EPITA, AIST, CNRS.
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

#ifndef ROBOPTIM_TRAJECTORY_VISUALIZATION_CUBIC_B_SPLINE_MATPLOTLIB_HH
# define ROBOPTIM_TRAJECTORY_VISUALIZATION_CUBIC_B_SPLINE_MATPLOTLIB_HH
# include <roboptim/trajectory/sys.hh>
# include <roboptim/trajectory/cubic-b-spline.hh>

# include <roboptim/core/visualization/matplotlib.hh>
# include <roboptim/core/visualization/matplotlib-function.hh>
# include <roboptim/core/visualization/matplotlib-commands.hh>

namespace roboptim
{
  namespace trajectory
  {
    namespace visualization
    {
      namespace matplotlib
      {
        typedef roboptim::visualization::matplotlib::Command Command;

        Command plot_spline(const CubicBSpline& spline,
                            CubicBSpline::value_type step)
        {
          using boost::format;
          using namespace detail;
          using namespace ::roboptim::visualization;

          std::stringstream ss;
          std::string data_name = spline.getName();
          long dimension = spline.outputSize();
          double start = spline.getLowerBound(spline.timeRange ());
          double end = spline.getUpperBound(spline.timeRange ());
          double nbp = static_cast<double> (spline.getNumberControlPoints());
          double intervals = nbp - 3;
          CubicBSpline::polynomials3vector_t polynomials;
          double t;

          if (dimension != 2)
	    {
	      throw std::string("This tool is only designed to print 2D splines");
	    }

          double inc = (end-start)/intervals;
          spline.toPolynomials(polynomials);

          if (inc <= step)
	    {
	      throw std::string("Inadapted step, should be much smaller");
	    }
          ss << data_name << " = np.array ([";

          for (t = start; t < end; t += step)
	    {
	      Function::vector_t res = spline (t);

	      ss << (boost::format ("(%2.8f") % normalize (t)).str ();

	      for (Function::size_type i = 0; i < spline.outputSize(); ++i)
		{
		  ss << (boost::format (", %2.8f") % normalize (res[i])).str();
		}
	      ss << "), ";
	    }

          Function::vector_t res = spline (end);

          ss << (boost::format ("(%2.8f") % normalize (end)).str ();

          for (Function::size_type i = 0; i < spline.outputSize(); ++i)
	    {
	      ss << (boost::format (", %2.8f") % normalize (res[i])).str();
	    }
          ss << ")";

          ss << "])" << std::endl;

          ss << "CP = np.array([";

          for (int i = 0; i < nbp; ++i)
            ss << (boost::format("(%1%, %2%),")
                   % spline.parameters()[2*i]
                   % spline.parameters()[2*i+1]).str();

          ss << "])" << std::endl;

          ss << (boost::format("plt.plot(%1%[:,1], %1%[:,2], label=\"%1%\", color='b')")
                 % data_name).str() << std::endl;

          for (unsigned long i = 0; i < intervals; ++i, start += inc)
	    {
	      Function::interval_t interval = Function::makeInterval(start, start+inc);
	      double bound_min = polynomials[i].min(interval).second;
	      double bound_max = polynomials[i].max(interval).second;
	      ss << (boost::format("plt.hlines(%1%, %2%, %3%, lw=2, color='c')")
		     % bound_min
		     % spline(interval.first)[0]
		     % spline(interval.second)[0]).str() << std::endl;

	      ss << (boost::format("plt.hlines(%1%, %2%, %3%, lw=2, color='g')")
		     % bound_max
		     % spline(interval.first)[0]
		     % spline(interval.second)[0]).str() << std::endl;

	      ss << (boost::format("plt.axvline(%1%, ls='--', color='gray')")
		     % spline(interval.first)[0]).str() << std::endl;
	    }

          ss << "plt.plot(CP[:,0], CP[:,1], '*', ls=':', lw=3, color='m', label=\"Control Points\")"
	     << std::endl;

          return Command(ss.str(), true);
        }
      }
    }
  }
}

#endif //! ROBOPTIM_TRAJECTORY_VISUALIZATION_CUBIC_B_SPLINE_MATPLOTLIB_HH