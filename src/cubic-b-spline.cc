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

#include <boost/numeric/conversion/converter.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <roboptim/trajectory/sys.hh>

#include <roboptim/core/indent.hh>
#include <roboptim/trajectory/cubic-b-spline.hh>

namespace roboptim
{

  //FIXME: defined_lc_in has to be true (false untested).
  CubicBSpline::CubicBSpline (interval_t tr, size_type outputSize,
			      const vector_t& p,
			      std::string name)
    throw ()
    : Trajectory<3> (tr, outputSize, p, name),
      nbp_ (p.size () / outputSize)
  {
    //Parameter size should be a multiple of spline dimension
    assert (parameters_.size () % outputSize == 0);
    // number of control points should be at least 4.
    assert (nbp_ >= 4);

    setParameters (p);
  }

  CubicBSpline::CubicBSpline (const CubicBSpline& spline) throw ()
    : Trajectory<3> (spline.timeRange (), spline.outputSize (),
		     spline.parameters ()),
      nbp_ (spline.parameters ().size () / spline.outputSize ())
  {
    //Parameter size should be a multiple of spline dimension
    assert (parameters_.size () % outputSize () == 0);
    // number of control points should be at least 4.
    assert (nbp_ >= 5);

    setParameters (spline.parameters ());
  }


  CubicBSpline::~CubicBSpline () throw ()
  {
  }

  void
  CubicBSpline::setParameters (const vector_t& p) throw ()
  {
    assert (p.size () == parameters_.size ());
    parameters_ = p;
  }

  void
  CubicBSpline::impl_compute (result_t& derivative, double t) const throw ()
  {
    t = detail::fixTime (t, *this);
    assert (timeRange ().first <= t && t <= timeRange ().second);
    this->derivative (derivative, t, 0);
  }

  CubicBSpline::value_type
  CubicBSpline::Dt () const
  {
    return length () / (nbp_ - 3.);
  }

  CubicBSpline::size_type
  CubicBSpline::interval (value_type t) const
  {
    typedef boost::numeric::converter<size_type, double> Double2SizeType;
    size_type i = Double2SizeType::convert
      (std::floor
       (3 + (t - getLowerBound (timeRange ())) / Dt ()));
    if (i == nbp_)
      i--;
    return i;
  }

  CubicBSpline::vector_t
  CubicBSpline::basisFunctions (value_type t, size_type order) const
  {
    t = detail::fixTime (t, *this);

    const double Dt = this->Dt ();
    const double t3 = getLowerBound (timeRange ());
    //const double tm = getUpperBound (timeRange ());

    const size_type i = interval (t);
    //const size_type n = outputSize ();

    // Non zero basis functions are b_{i-3,3}(t), b_{i-2,3}(t), b_{i-1,3}(t),
    // b_{i,3}(t).
    const value_type t_i = t3 + (i - 3) * Dt;

    const value_type tau_i = t - t_i;
    const value_type tau_i_2 = tau_i * tau_i;
    const value_type tau_i_3 = tau_i * tau_i_2;
    assert (tau_i <= Dt + tolerance ());

    const value_type Dt_2 = Dt * Dt;
    const value_type Dt_3 = Dt * Dt_2;

    // Evaluate basis functions, up to division by 6 Dt_3 that will
    // be performed at the end.
    vector_t b_i (4);
    switch (order)
      {
      case 0:
	{
	  b_i[3] = (  -tau_i_3 + 3*Dt*tau_i_2 - 3*Dt_2*tau_i +   Dt_3);
	  b_i[2] = ( 3*tau_i_3 - 6*Dt*tau_i_2                + 4*Dt_3);
	  b_i[1] = (-3*tau_i_3 + 3*Dt*tau_i_2 + 3*Dt_2*tau_i +   Dt_3);
	  b_i[0] = (   tau_i_3                                       );
	  b_i /= 6.;
	}
	break;
      case 1:
	{
	  b_i[3] = (-3*tau_i_2 +  6*Dt*tau_i - 3*Dt_2);
	  b_i[2] = ( 9*tau_i_2 - 12*Dt*tau_i         );
	  b_i[1] = (-9*tau_i_2 +  6*Dt*tau_i + 3*Dt_2);
	  b_i[0] = ( 3*tau_i_2                       );
	  b_i /= 6.;
	}
	break;
      case 2:
	{
	  b_i[3] = (  -tau_i +  Dt);
	  b_i[2] = ( 3*tau_i - 2*Dt);
	  b_i[1] = (-3*tau_i +   Dt);
	  b_i[0] = (   tau_i       );
	  // division by six not needed here.
	}
	break;
      default:
	assert (0);
      }
    return b_i;
  }

  void
  CubicBSpline::impl_derivative (gradient_t& derivative, double t,
				 size_type order)
    const throw ()
  {
    using boost::numeric::ublas::subrange;

    t = detail::fixTime (t, *this);
    const size_type i = interval (t);
    const size_type n = outputSize ();
    const value_type Dt_3 = Dt () * Dt () * Dt ();

    const vector_t& P_i_3 = subrange (parameters (), (i - 3) * n, (i - 2) * n);
    const vector_t& P_i_2 = subrange (parameters (), (i - 2) * n, (i - 1) * n);
    const vector_t& P_i_1 = subrange (parameters (), (i - 1) * n, (i + 0) * n);
    const vector_t& P_i   = subrange (parameters (), (i - 0) * n, (i + 1) * n);

    const vector_t b_i = basisFunctions (t, order);

    derivative = b_i[3] * P_i_3 + b_i[2] * P_i_2 + b_i[1] * P_i_1 + b_i[0] * P_i;
    derivative /= Dt_3;
  }

  void
  CubicBSpline::impl_derivative (gradient_t& derivative,
				 StableTimePoint stp,
				 size_type order) const throw ()
  {
    this->impl_derivative (derivative,
			   stp.getTime (this->timeRange ()),
			   order);
  }

  CubicBSpline::jacobian_t
  CubicBSpline::variationConfigWrtParam (double t) const throw ()
  {
    return variationDerivWrtParam (t, 0);
  }

  CubicBSpline::jacobian_t
  CubicBSpline::variationDerivWrtParam (double t, size_type order)
    const throw ()
  {
    using boost::numeric::ublas::noalias;
    using boost::numeric::ublas::subrange;

    const size_type i = interval (t);
    const size_type n = outputSize ();
    const vector_t b_i = basisFunctions (t, order);
    const value_type Dt_3 = Dt () * Dt () * Dt ();

    jacobian_t jac = ublas::zero_matrix<double> (n, nbp_ * n);
    const ublas::identity_matrix<double> In (n);

    noalias (subrange (jac, 0, n, (i - 3) * n, (i - 2) * n)) = b_i[3] * In;
    noalias (subrange (jac, 0, n, (i - 2) * n, (i - 1) * n)) = b_i[2] * In;
    noalias (subrange (jac, 0, n, (i - 1) * n,       i * n)) = b_i[1] * In;
    noalias (subrange (jac, 0, n,       i * n, (i + 1) * n)) = b_i[0] * In;
    jac  /= Dt_3;
    return jac;
  }

  CubicBSpline::jacobian_t
  CubicBSpline::variationConfigWrtParam (StableTimePoint stp)
    const throw ()
  {
    return this->variationConfigWrtParam (stp.getTime (this->timeRange ()));
  }


  CubicBSpline::jacobian_t
  CubicBSpline::variationDerivWrtParam (StableTimePoint stp, size_type order)
    const throw ()
  {
    return this->variationDerivWrtParam
      (stp.getTime (this->timeRange ()), order);
  }


  CubicBSpline::value_type
  CubicBSpline::singularPointAtRank (size_type rank) const
  {
    return rank * length () / (nbp_- 3);
  }

  CubicBSpline::vector_t
  CubicBSpline::derivBeforeSingularPoint (size_type rank, size_type order) const
  {
    return derivative (singularPointAtRank (rank), order);
  }

  CubicBSpline::vector_t
  CubicBSpline::derivAfterSingularPoint (size_type rank, size_type order) const
  {
    return derivative (singularPointAtRank (rank), order);
  }

  std::ostream&
  CubicBSpline::print (std::ostream& o) const throw ()
  {
    o << "CubicBSpline" << incindent << std::endl
      << "Number of parameters per spline function: " << nbp_ << std::endl
      << "Length: " << length () << std::endl
      << "Parameters: " << parameters ()
      << decindent;
    return o;
  }
} // end of namespace roboptim.
