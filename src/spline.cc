
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

/**
 * \brief Class Spline implementation.
 */

#include <roboptim-trajectory/spline.hh>

#include <spline/bspline.h>

namespace roboptim
{
  //FIXME: defined_lc_in has to be true (false untested).
  Spline::Spline (size_type m, const vector_t& p, int nbP) throw ()
    : Trajectory<4> (m, p),
      spline_ (),
      nbp_ (nbP)
  {
    //FIXME: check params here.
    spline_ = new bspline (m, nbP + 4, 1, true, true, true);

    vector_t pos_init (m);
    vector_t final_pos (m);
    double duree_mvt;
    matrix_t mp (m, nbP - 2);

    spline_->convert_parameters_x2P (&parameters_[0],
				     &mp,
				     pos_init,
				     final_pos,
				     duree_mvt);

    spline_->def_parameters(&mp, pos_init, final_pos, duree_mvt);
  }

  Spline::~Spline () throw ()
  {
    delete spline_;
  }

  Spline::vector_t
  Spline::operator () (double t) const throw ()
  {
    vector_t res (m);
    spline_->calc_fun (t, &res);
    return res;
  }

  Spline::vector_t
  Spline::derivative (double t, size_type order) const throw ()
  {
    vector_t res (m);
    switch (order)
      {
      case 0:
	return operator () (t);
      case 1:
	spline_->calc_dfun (t, &res);
	break;

      case 2:
	spline_->calc_ddfun(t, &res);
	break;
      default:
	assert (0);
      }
    return res;
  }

  Spline::jacobian_t
  Spline::variationConfigWrtParam (double t) const throw ()
  {
    jacobian_t jac (m, nbp_ * m);
    matrix_t fun (m, 1);

    //FIXME: change by two points if required.
    //vector_t all_t (2);
    vector_t all_t (1);
    all_t[0] = t;
    //res[1] = t + 1.;

    ublas::matrix<ublas::vector<double> > fun_grad (m, 1);
    for (size_type i = 0; i < m; ++i)
      fun_grad (i, 0).resize (5);

    spline_->calc_fun_grad(&all_t, &fun, &fun_grad, 1);

//     matrix_t tmp (m, 5);
//     uncompress_grad (t, &fun_grad, &grad);

    return jac;
  }

  Spline::jacobian_t
  Spline::variationDerivWrtParam (double t, size_type order)
    const throw ()
  {
    jacobian_t jac (m, m);
    return jac;
  }

  Spline::value_type
  Spline::singularPointAtRank (size_type rank) const
  {
    return 42.;
  }

  Spline::vector_t
  Spline::derivBeforeSingularPoint (size_type rank, size_type order) const
  {
    vector_t res (42);
    return res;
  }

  Spline::vector_t
  Spline::derivAfterSingularPoint (size_type rank, size_type order) const
  {
    vector_t res (42);
    return res;
  }

} // end of namespace roboptim.
