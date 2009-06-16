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
# include <boost/make_shared.hpp>

# include <roboptim/core/numeric-linear-function.hh>

namespace roboptim
{
  template <typename F, typename CLIST, typename C>
  Freeze<F, CLIST, C>::Freeze (problem_t& problem,
			       const frozenArguments_t fa) throw ()
    : problem_ (problem),
      frozenArguments_ (fa)
  {
  }

  template <typename F, typename CLIST, typename C>
  Freeze<F, CLIST, C>::~Freeze () throw ()
  {
  }

  template <typename F, typename CLIST, typename C>
  void
  Freeze<F, CLIST, C>::operator () () throw ()
  {
    using namespace boost;
    typedef frozenArguments_t::const_iterator citer_t;

    for (citer_t it = frozenArguments_.begin ();
	 it != frozenArguments_.end (); ++it)
      {
	assert (it->first < problem_.function ().inputSize ());

	Function::matrix_t a (1, problem_.function ().inputSize ());
	Function::vector_t b (1);

	a.clear (), b.clear ();

	a(0, it->first) = 1.;

	b[0] = -it->second;

	shared_ptr<C> ptr = static_pointer_cast<C>
	  (make_shared<NumericLinearFunction> (a, b));
	this->problem_.addConstraint (ptr,
				      Function::makeInterval (0., 0.));
      }
  }

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_FREEZE_HXX
