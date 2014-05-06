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

#ifndef ROBOPTIM_TRAJECTORY_FREEZE_HXX
# define ROBOPTIM_TRAJECTORY_FREEZE_HXX
# include <roboptim/core/numeric-linear-function.hh>

namespace roboptim
{
  template <typename P>
  Freeze<P>::Freeze (problem_t& problem)
    : problem_ (problem)
  {
  }

  template <typename P>
  Freeze<P>::~Freeze ()
  {
  }

  template <typename P>
  void
  Freeze<P>::operator () (const frozenArguments_t frozenArguments)
  {
    using namespace boost;
    typedef frozenArguments_t::const_iterator citer_t;

    for (citer_t it = frozenArguments.begin ();
	 it != frozenArguments.end (); ++it)
      {
	assert (it->first < problem_.function ().inputSize ());

	Function::interval_t& interval =
	  this->problem_.argumentBounds()[it->first];

	Function::value_type min =
	  std::max (Function::getLowerBound (interval), it->second);
	Function::value_type max =
	  std::min (Function::getUpperBound (interval), it->second);
	this->problem_.argumentBounds()[it->first] = Function::makeInterval (min, max);
      }
  }

  template <typename P>
  void
  Freeze<P>::operator () (const std::vector<Function::size_type>& indices,
			  const Function::vector_t& values)
  {
    frozenArguments_t fa;
    for (unsigned i = 0; i < indices.size (); ++i)
      {
	assert (indices[i] < values.size ());
	fa.push_back (std::make_pair (indices[i], values[indices[i]]));
      }
    (*this) (fa);
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_FREEZE_HXX
