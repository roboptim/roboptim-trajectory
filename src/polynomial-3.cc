// Copyright (C) 2012 by Florent Lamiraux, CNRS.
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

#include "roboptim/trajectory/polynomial-3.hh"

#include <boost/lexical_cast.hpp>

namespace roboptim {

  Polynomial3::Polynomial3 ()
  {
    t0_ = 0;
    memset (coefs_, 0, 4*sizeof(double));
  }

  Polynomial3::Polynomial3 (double t0, double a0, double a1, double a2,
                            double a3)
  {
    t0_ = t0;
    coefs_ [0] = a0;
    coefs_ [1] = a1;
    coefs_ [2] = a2;
    coefs_ [3] = a3;
  }

  Polynomial3 Polynomial3::operator* (const Polynomial3& poly) const
  {
    Polynomial3 other = poly.translate (t0_);
    assert (coefs_ [3] * other.coefs_ [1] +
	    coefs_ [2] * other.coefs_ [2] +
	    coefs_ [1] * other.coefs_ [3] == 0);

    assert (coefs_ [3] * other.coefs_ [2] +
	    coefs_ [2] * other.coefs_ [3] == 0);

    assert (coefs_ [3] * other.coefs_ [3] == 0);

    return Polynomial3 (t0_, coefs_ [0] * other.coefs_ [0],
			coefs_ [1] * other.coefs_ [0] +
			coefs_ [0] * other.coefs_ [1],
			coefs_ [2] * other.coefs_ [0] +
			coefs_ [1] * other.coefs_ [1] +
			coefs_ [0] * other.coefs_ [2],
			coefs_ [3] * other.coefs_ [0] +
			coefs_ [2] * other.coefs_ [1] +
			coefs_ [1] * other.coefs_ [2] +
			coefs_ [0] * other.coefs_ [3]);
  }

  Polynomial3 Polynomial3::operator+ (const Polynomial3& poly) const
  {
    Polynomial3 other = poly.translate (t0_);
    return Polynomial3 (t0_,
			coefs_ [0] + other.coefs_ [0],
			coefs_ [1] + other.coefs_ [1],
			coefs_ [2] + other.coefs_ [2],
			coefs_ [3] + other.coefs_ [3]);
  }

  Polynomial3 Polynomial3::operator- (const Polynomial3& poly) const
  {
    Polynomial3 other = poly.translate (t0_);
    return Polynomial3 (t0_,
			coefs_ [0] - other.coefs_ [0],
			coefs_ [1] - other.coefs_ [1],
			coefs_ [2] - other.coefs_ [2],
			coefs_ [3] - other.coefs_ [3]);
  }

  Polynomial3 operator* (const double& lambda, const Polynomial3& poly)
  {
    return Polynomial3 (poly.t0_,
			lambda * poly.coefs_ [0],
			lambda * poly.coefs_ [1],
			lambda * poly.coefs_ [2],
			lambda * poly.coefs_ [3]);
  }

  Polynomial3 Polynomial3::translate (const double &t1) const
  {
    double dt = t1 - t0_;
    double dt2 = dt*dt;
    double dt3 = dt*dt2;
    double a0 = coefs_ [0];
    double a1 = coefs_ [1];
    double a2 = coefs_ [2];
    double a3 = coefs_ [3];

    return Polynomial3 (t1, a0 + a1*dt + a2*dt2 + a3*dt3,
			a1 + 2*dt*a2 + 3*dt2*a3,
			a2 + 3*dt*a3,
			a3);
  }

  double Polynomial3::operator () (const double& t) const
  {
    double dt = t - t0_;
    double dt2 = dt*dt;
    double dt3 = dt*dt2;
    double a0 = coefs_ [0];
    double a1 = coefs_ [1];
    double a2 = coefs_ [2];
    double a3 = coefs_ [3];

    return a0 + a1*dt + a2*dt2 + a3*dt3;
  }

  double Polynomial3::derivative (const double& t, size_type order) const
  {
    double dt = t - t0_;
    double dt2 = dt*dt;
    double a1 = coefs_ [1];
    double a2 = coefs_ [2];
    double a3 = coefs_ [3];

    switch (order) {
    case 0:
      return (*this) (t);
    case 1:
      return a1 + 2*a2*dt + 3*a3*dt2;
    case 2:
      return 2*a2 + 6*a3*dt;
    case 3:
      return 6*a3;
    default:
      return 0;
    }
  }

  std::ostream&
  Polynomial3::print (std::ostream& o) const throw ()
  {
    for (unsigned int i = 0; i < 4; ++i)
      {
        o << ((coefs_[i] >= 0)? " + " : " - ")
          << fabs(coefs_[i]) << " * "
          << ((fabs(t0_) < 1e-12)? "t":
              std::string("(t - " )
              + boost::lexical_cast<std::string>(t0_)
              + ")")
          << "^" << i;
      }
    return o;
  }

  std::ostream& operator << (std::ostream& o, const Polynomial3& p)
  {
    return p.print(o);
  }

} // namespace roboptim
