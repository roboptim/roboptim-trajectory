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

#ifndef ROBOPTIM_TRAJECTORY_FREEZE_HH
# define ROBOPTIM_TRAJECTORY_FREEZE_HH
# include <roboptim/trajectory/sys.hh>

# include <vector>
# include <utility>

# include <roboptim/core/fwd.hh>
# include <roboptim/core/problem.hh>

namespace roboptim
{
namespace trajectory
{
  /// \brief Add constraints that freeze parameters.
  ///
  /// This constraint expects a vector of pairs (argument id, value) to
  /// build a function that will freeze the given arguments to their
  /// associated value.
  ///
  /// For instance, the vector: [(0, 5.), (3, -12.)] forces the first
  /// parameter to five and the fourth one to minus twelve.
  ///
  /// This adds to the problem one linear constraint per frozen parameter.
  template <typename P>
  class Freeze
  {
  public:
    /// \brief Pair representing an argument index and a value.
    typedef std::pair<Function::size_type,
		      Function::value_type> frozenArgument_t;

    /// \brief Problem type.
    typedef P problem_t;

    /// \brief Vector of pairs (argument index, value).
    ///
    /// This type define what are the frozen arguments and what their
    /// value.
    typedef std::vector<frozenArgument_t> frozenArguments_t;

    /// \brief Create the constraint from a vector of pairs.
    ///
    /// \param problem problem that will be modified.
    Freeze (problem_t& problem);

    virtual ~Freeze ();

    /// \brief Apply modification.
    /// \param fa Vector of pairs containing what to freeze and to what value.
    void operator () (const frozenArguments_t fa);

    /// \brief Apply modification.
    /// \param indices Vector of parameters index that will be frozen.
    /// \param values Vector of parameters values.
    void operator () (const std::vector<Function::size_type>& indices,
		      const Function::vector_t& values);

  private:
    /// \brief Reference to the problem that will be modified.
    problem_t& problem_;
  };

  template <typename P>
  Freeze<P> makeFreeze (P& problem)
  {
    return Freeze<P> (problem);
  }

  /// Example shows freeze use.
  /// \example spline-optimization.cc

} // end of namespace trajectory.
} // end of namespace roboptim.

# include <roboptim/trajectory/freeze.hxx>
#endif //! ROBOPTIM_TRAJECTORY_FREEZE_HH
