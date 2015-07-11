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

#ifndef ROBOPTIM_TRAJECTORY_VISUALIZATION_B_SPLINE_MATPLOTLIB_HH
# define ROBOPTIM_TRAJECTORY_VISUALIZATION_B_SPLINE_MATPLOTLIB_HH

# include <exception>

# include <roboptim/trajectory/sys.hh>
# include <roboptim/trajectory/b-spline.hh>

# include <roboptim/core/visualization/matplotlib.hh>
# include <roboptim/core/visualization/matplotlib-function.hh>
# include <roboptim/core/visualization/matplotlib-commands.hh>

# include <roboptim/trajectory/visualization/matplotlib.hh>

namespace roboptim
{
  namespace trajectory
  {
    namespace visualization
    {
      namespace matplotlib
      {
        typedef roboptim::visualization::matplotlib::Command Command;

        template <int N>
        Command plot_spline(const BSpline<N>& spline,
                            typename BSpline<N>::value_type step)
        {
          using boost::format;
          using namespace ::roboptim::visualization;

          std::stringstream ss;
          std::string data_name = detail::formattedVarName (spline.getName());
          long dimension = spline.outputSize();
          double start = spline.getLowerBound(spline.timeRange ());
          double end = spline.getUpperBound(spline.timeRange ());
          double nbp = static_cast<double> (spline.getNumberControlPoints());
          double intervals = nbp - N;
          typename BSpline<N>::basisPolynomials_t polynomials;
          const typename BSpline<N>::vector_t& kv = spline.knotVector();
          double t;

          if (dimension > 2)
	    {
	      throw std::runtime_error ("This tool is not designed to print splines of degree > 2");
	    }

          spline.toPolynomials(polynomials);

          ss << data_name << " = np.array ([";

          for (t = start; t < end; t += step)
	    {
	      typename BSpline<N>::vector_t res = spline (t);

	      ss << (boost::format ("(%2.8f") % normalize (t)).str ();

	      for (typename BSpline<N>::size_type i = 0; i < spline.outputSize(); ++i)
		{
		  ss << (boost::format (", %2.8f") % normalize (res[i])).str();
		}
	      ss << "), ";
	    }

          typename BSpline<N>::vector_t res = spline (end);

          ss << (boost::format ("(%2.8f") % normalize (end)).str ();

          for (typename BSpline<N>::size_type i = 0; i < spline.outputSize(); ++i)
	    {
	      ss << (boost::format (", %2.8f") % normalize (res[i])).str();
	    }
          ss << ")";

          ss << "])" << std::endl;

          ss << "CP = np.array([";

          if (dimension == 2)
	    for (int i = 0; i < nbp; ++i)
	      ss << (boost::format("(%1%, %2%),")
		     % spline.parameters()[2*i]
		     % spline.parameters()[2*i+1]).str();

          ss << "])" << std::endl;

          if (dimension == 2)
	    {
	      ss << (boost::format("plt.plot(%1%[:,1], %1%[:,2], label=\"%2%\", color='b')")
		     % data_name % spline.getName ()).str() << std::endl;
	    }
          else
	    {
	      ss << (boost::format("plt.plot(%1%[:,0], %1%[:,1], label=\"%2%\", color='b')")
		     % data_name % spline.getName ()).str() << std::endl;
	    }

          for (unsigned long i = 0; i < intervals; ++i)
	    {
	      long n = static_cast<long>(i+N);
	      double inc = kv[n+1] - kv[n]; //Actual inc in time
	      assert(step < inc);

	      Function::interval_t interval = Function::makeInterval(start, start+inc);
	      double bound_min = polynomials[i].min(interval).second;
	      double bound_max = polynomials[i].max(interval).second;

	      double value_min = (dimension==2) ? spline(interval.first)[0] : interval.first;
	      double value_max = (dimension==2) ? spline(interval.second)[0] : interval.second;
	      ss << (boost::format("plt.hlines(%1%, %2%, %3%, lw=2, color='c')")
		     % bound_min
		     % value_min
		     % value_max).str() << std::endl;

	      ss << (boost::format("plt.hlines(%1%, %2%, %3%, lw=2, color='g')")
		     % bound_max
		     % value_min
		     % value_max).str() << std::endl;

	      ss << (boost::format("plt.axvline(%1%, ls='--', color='gray')")
		     % value_min).str() << std::endl;
	      start += inc;
	    }

          if (dimension == 2)
	    ss << "plt.plot(CP[:,0], CP[:,1], '*', ls=':', lw=3, color='m', label=\"Control Points\")" << std::endl;

          return Command(ss.str(), true);
        }
      }
    }
  }
}

#endif //! ROBOPTIM_TRAJECTORY_VISUALIZATION_B_SPLINE_MATPLOTLIB_HH
