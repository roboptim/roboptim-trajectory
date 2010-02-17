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

  void
  CubicBSpline::impl_derivative (gradient_t& derivative, double t,
				 size_type order)
    const throw ()
  {
    t = detail::fixTime (t, *this);
    double t3 = timeRange ().first;
    double tm = timeRange ().second;
    assert (t3 <= t && t <= tm);
    // determine on which interval t lies
    double Dt = (tm-t3)/(nbp_-3);
    unsigned int i=(unsigned int) floor(3+(t-t3)/Dt);
    assert(i>=3);
    assert(i<=nbp_);
    if (i == nbp_) i--;
    // Non zero basis functions are b_{i-3,3}(t), b_{i-2,3}(t), b_{i-1,3}(t),
    // b_{i,3}(t).
    double t_i = t3 + (i-3)*Dt;
    double tau_i = t - t_i;
    assert (tau_i <= Dt);
    double tau_i_2 = tau_i*tau_i;
    double tau_i_3 = tau_i*tau_i_2;
    double Dt_2 = Dt*Dt;
    double Dt_3 = Dt*Dt_2;
    unsigned int n = outputSize();
    try {
    ublas::vector<double> P_i_3 =
      boost::numeric::ublas::subrange(parameters(),
				      (i-3)*n,
				      (i-2)*n);
    ublas::vector<double> P_i_2 =
      boost::numeric::ublas::subrange(parameters(),
				      (i-2)*n,
				      (i-1)*n);
    ublas::vector<double> P_i_1 =
      boost::numeric::ublas::subrange(parameters(),
				      (i-1)*n,
				      (i)*n);
    ublas::vector<double> P_i =
      boost::numeric::ublas::subrange(parameters(),
				      (i)*n,
				      (i+1)*n);
    switch (order)
      {
      case 0:
	{
	  // Evaluate basis functions, up to division by 6 Dt_3 that will
	  // be performed at the end.
	  double b_i_3 = (  -tau_i_3 + 3*Dt*tau_i_2 - 3*Dt_2*tau_i +   Dt_3);
	  double b_i_2 = ( 3*tau_i_3 - 6*Dt*tau_i_2                + 4*Dt_3);
	  double b_i_1 = (-3*tau_i_3 + 3*Dt*tau_i_2 + 3*Dt_2*tau_i +   Dt_3);
	  double b_i   = (   tau_i_3                                       );
	  derivative =
	    (b_i_3*P_i_3 + b_i_2*P_i_2 + b_i_1*P_i_1 + b_i*P_i)/(6*Dt_3);
	}
	break;
      case 1:
	{
	  // Evaluate basis function derivatives, up to division by
	  // 6 Dt_3 that will be performed at the end.
	  double b_i_3 = (-3*tau_i_2 +  6*Dt*tau_i - 3*Dt_2);
	  double b_i_2 = ( 9*tau_i_2 - 12*Dt*tau_i         );
	  double b_i_1 = (-9*tau_i_2 +  6*Dt*tau_i + 3*Dt_2);
	  double b_i   = ( 3*tau_i_2                       );
	  derivative =
	    (b_i_3*P_i_3 + b_i_2*P_i_2 + b_i_1*P_i_1 + b_i*P_i)/(6*Dt_3);
	}
	break;
      case 2:
	{
	  // Evaluate basis function second derivatives, up to division by
	  // 6 Dt_3 that will be performed at the end.
	  double b_i_3 = (  -tau_i +  Dt);
	  double b_i_2 = ( 3*tau_i - 2*Dt);
	  double b_i_1 = (-3*tau_i +   Dt);
	  double b_i   = (   tau_i       );
	  derivative =
	    (b_i_3*P_i_3 + b_i_2*P_i_2 + b_i_1*P_i_1 + b_i*P_i)/(Dt_3);
	}
	break;
      default:
	assert (0);
      }
    } catch (...) {
      assert(0);
    }
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
    return variationDerivWrtParam (t, 0.);
  }

  CubicBSpline::jacobian_t
  CubicBSpline::variationDerivWrtParam (double t, size_type order)
    const throw ()
  {
    t = detail::fixTime (t, *this);
    double t3 = timeRange ().first;
    double tm = timeRange ().second;
    assert (t3 <= t && t <= tm);
    // determine on which interval t lies
    double Dt = (tm-t3)/(nbp_-3);
    unsigned int i=(unsigned int) floor(3+(t-t3)/Dt);
    assert(i>=3);
    // Non zero basis functions are b_{i-3,3}(t), b_{i-2,3}(t), b_{i-1,3}(t),
    // b_{i,3}(t).
    double t_i = t3 + (i-3)*Dt;
    double tau_i = t - t_i;
    assert (tau_i < Dt);
    double tau_i_2 = tau_i*tau_i;
    double tau_i_3 = tau_i*tau_i_2;
    double Dt_2 = Dt*Dt;
    double Dt_3 = Dt*Dt_2;
    unsigned int n = outputSize();
    ublas::vector<double> P_i_3 =
      boost::numeric::ublas::subrange(parameters(),
				      (i-3)*n,
				      (i-2)*n);
    ublas::vector<double> P_i_2 =
      boost::numeric::ublas::subrange(parameters(),
				      (i-2)*n,
				      (i-1)*n);
    ublas::vector<double> P_i_1 =
      boost::numeric::ublas::subrange(parameters(),
				      (i-1)*n,
				      (i)*n);
    ublas::vector<double> P_i =
      boost::numeric::ublas::subrange(parameters(),
				      (i)*n,
				      (i+1)*n);
    jacobian_t jac = ublas::zero_matrix<double>(n, nbp_ * n);
    ublas::identity_matrix<double> In(n);
    double b_i=0, b_i_1=0, b_i_2=0, b_i_3=0;

    switch (order)
      {
      case 0:
	{
	  // Evaluate basis functions, up to division by 6 Dt_3 that will
	  // be performed at the end.
	  b_i_3 = (  -tau_i_3 + 3*Dt*tau_i_2 - 3*Dt_2*tau_i +   Dt_3);
	  b_i_2 = ( 3*tau_i_3 - 6*Dt*tau_i_2                + 4*Dt_3);
	  b_i_1 = (-3*tau_i_3 + 3*Dt*tau_i_2 + 3*Dt_2*tau_i +   Dt_3);
	  b_i   = (   tau_i_3                                       );
	}
	break;
      case 1:
	{
	  // Evaluate basis function derivatives, up to division by
	  // 6 Dt_3 that will be performed at the end.
	  b_i_3 = (-3*tau_i_2 +  6*Dt*tau_i - 3*Dt_2);
	  b_i_2 = ( 9*tau_i_2 - 12*Dt*tau_i         );
	  b_i_1 = (-9*tau_i_2 +  6*Dt*tau_i + 3*Dt_2);
	  b_i   = ( 3*tau_i_2                       );
	}
	break;
      case 2:
	{
	  // Evaluate basis function second derivatives, up to division by
	  // 6 Dt_3 that will be performed at the end.
	  b_i_3 = (  -tau_i +  Dt);
	  b_i_2 = ( 3*tau_i - 2*Dt);
	  b_i_1 = (-3*tau_i +   Dt);
	  b_i   = (   tau_i       );
	}
	break;
      default:
	assert (0);
      }
    noalias(subrange(jac, 0, n, (i-3)*n, (i-2)*n)) = (b_i_3)/(6*Dt_3)*In;
    noalias(subrange(jac, 0, n, (i-2)*n, (i-1)*n)) = (b_i_2)/(6*Dt_3)*In;
    noalias(subrange(jac, 0, n, (i-1)*n, i*n))   = (b_i_1)/(6*Dt_3)*In;
    noalias(subrange(jac, 0, n, i*n, (i+1)*n)) = (b_i)/(6*Dt_3)*In;
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
