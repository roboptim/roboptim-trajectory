
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

#include <roboptim/core/indent.hh>
#include <roboptim/trajectory/spline.hh>

#include <spline/bspline.h>

namespace roboptim
{

  //FIXME: defined_lc_in has to be true (false untested).
  Spline::Spline (interval_t tr, size_type outputSize, const vector_t& p,
		  std::string name)
    throw ()
    : Trajectory<4> (tr, outputSize, p, name),
      spline_ (),
      nbp_ (p.size () / outputSize)
  {
    // Can not work with a smalle parameters vector.
    assert (parameters_.size () >= 4);

    //
    assert (parameters_.size () >= 2 * outputSize
	    && parameters_.size () % outputSize == 0);

    //FIXME: check params here.
    spline_ = new bspline (outputSize, nbp_ + 4, 1, true, true, true);

    setParameters (p);
  }

  Spline::Spline (const Spline& spline) throw ()
    : Trajectory<4> (spline.timeRange (), spline.outputSize (),
		     spline.parameters ()),
      spline_ (),
      nbp_ (spline.parameters ().size () / spline.outputSize ())
  {
    // Can not work with a smaller parameters vector.
    assert (parameters_.size () >= 4);

    assert (parameters_.size () >= 2 * spline.outputSize ()
	    && parameters_.size () % spline.outputSize () == 0);

    //FIXME: check params here.
    spline_ = new bspline (spline.outputSize (), nbp_ + 4, 1,
			   true, true, true);

    setParameters (spline.parameters ());
  }


  Spline::~Spline () throw ()
  {
    delete spline_;
  }

  void
  Spline::setParameters (const vector_t& p) throw ()
  {
    assert (p.size () == parameters_.size ());
    parameters_ = p;

    // Initialized by convert_parameters.
    vector_t pos_init (outputSize ());
    vector_t final_pos (outputSize ());
    double l = 0.;
    matrix_t mp (outputSize (), nbp_ - 2);

    vector_t splineParams = makeBSplineVector ();

    spline_->convert_parameters_x2P (&splineParams[0],
				     &mp,
				     pos_init,
				     final_pos,
				     l);
    spline_->def_parameters (&mp, pos_init, final_pos, l);
  }

  void
  Spline::impl_compute (result_t& derivative, double t) const throw ()
  {
    t = detail::fixTime (t, *this);
    assert (timeRange ().first <= t && t <= timeRange ().second);
    this->derivative (derivative, t, 0);
  }

  void
  Spline::impl_derivative (gradient_t& derivative, double t, size_type order)
    const throw ()
  {
    t = detail::fixTime (t, *this);
    assert (timeRange ().first <= t && t <= timeRange ().second);
    switch (order)
      {
      case 0:
	spline_->calc_fun (t, &derivative);
	break;
      case 1:
	spline_->calc_dfun (t, &derivative);
	break;
      case 2:
	spline_->calc_ddfun(t, &derivative);
	break;
      default:
	assert (0);
      }
  }

  Spline::jacobian_t
  Spline::variationConfigWrtParam (double t) const throw ()
  {
    return variationDerivWrtParam (t, 0.);
  }

  Spline::jacobian_t
  Spline::variationDerivWrtParam (double t, size_type order)
    const throw ()
  {
    t = detail::fixTime (t, *this);
    assert (timeRange ().first <= t && t <= timeRange ().second);
    matrix_t fun (outputSize (), 1);

    vector_t all_t (1);
    all_t[0] = t;

    ublas::matrix<ublas::vector<double> > fun_grad (outputSize (), 1);
    for (size_type i = 0; i < outputSize (); ++i)
      fun_grad (i, 0).resize (5);

    switch (order)
      {
      case 0:
	spline_->calc_fun_grad (&all_t, &fun, &fun_grad, 1);
	break;
      case 1:
	spline_->calc_dfun_grad (&all_t, &fun, &fun_grad, 1);
	break;

      case 2:
	spline_->calc_ddfun_grad (&all_t, &fun, &fun_grad, 1);
	break;
      default:
	assert (0);
      }

    ublas::matrix<ublas::vector<double> > tmp (outputSize (), 1);
    for (size_type i = 0; i < outputSize (); ++i)
      tmp (i, 0).resize (parameters_.size ()+1);

    spline_->uncompress_grad (&all_t,
	                      &fun_grad,
			      &tmp,
			      1);

    jacobian_t jac (outputSize (), nbp_ * outputSize ());
    for (size_type i = 0; i < outputSize (); ++i)
      for (size_type j = 0; j < nbp_ * outputSize (); ++j)
	jac (i, j) = tmp (i, 0)[j];

    return jac;
  }

  Spline::value_type
  Spline::singularPointAtRank (size_type rank) const
  {
    return rank * length () / spline_->get_Nint ();
  }

  Spline::vector_t
  Spline::derivBeforeSingularPoint (size_type rank, size_type order) const
  {
    return derivative (singularPointAtRank (rank), order);
  }

  Spline::vector_t
  Spline::derivAfterSingularPoint (size_type rank, size_type order) const
  {
    return derivative (singularPointAtRank (rank), order);
  }

  Spline::vector_t
  Spline::makeBSplineVector ()
  {
    vector_t res (parameters_.size () + 1);

    for (size_type i = 0; i < parameters_.size (); ++i)
      res[i] = parameters_[i];
    res[parameters_.size ()] = length ();
    return res;
  }

  std::ostream&
  Spline::print (std::ostream& o) const throw ()
  {
    o << "Spline" << incindent << std::endl
      << "Number of parameters per spline function: " << nbp_ << std::endl
      << "Length: " << length ()
      << decindent;
    return o;
  }
} // end of namespace roboptim.
