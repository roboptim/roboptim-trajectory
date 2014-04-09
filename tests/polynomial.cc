/*
 * polynomial_test.cc
 *
 *  Created on: Mar 13, 2013
 *      Author: wern_al
 */

# include <roboptim/trajectory/polynomial.hh>
# include <roboptim/trajectory/polynomial-3.hh>

# include <cstdlib>
# include <limits>
# include <cassert>
# include <cmath>
# include <iostream>

using namespace roboptim;
using namespace std;

#undef NDEBUG



template <int N>
void test_multiply(){
  Eigen::Matrix<double,N+1,1> params;
  params.setRandom();
  double t0 = (double)rand()/RAND_MAX;
  double t1 = (double)rand()/RAND_MAX;
  Polynomial<N> p_1(t0,params);
  p_1(0.1);

  params.setZero();
  params[0]=0.3;
  Polynomial<N> p_2(t1,params);
  p_2(0.1);

  Polynomial<N> p_3 = p_1 * p_2;
  p_3(0.1);
}




template <int N>
struct poly_props{
  static void check_derivative(Polynomial<N> const &, double , int , double){
    //do nothing, not implemented
  };
  static void check_translate (Polynomial<N> const &, const double &, Polynomial<N> const &){
    //do nothing, not implemented
  }
  static void check_evaluate(double t0, Eigen::Matrix<double,3+1,1> &params, double t){
    //do nothing, not implemented
  }
};

/**
 * code from Polynomial3 class
 */

template <>
void poly_props<3>::check_derivative(Polynomial<3> const & poly, double t, int order, double derivative){
  double dt = t - poly.t0_;
  double dt2 = dt*dt;
  double a1 = poly.coefs_ [1];
  double a2 = poly.coefs_ [2];
  double a3 = poly.coefs_ [3];

  double result;
  switch (order) {
  case 0:
    result = poly (t);break;
  case 1:
    result = a1 + 2*a2*dt + 3*a3*dt2;break;
  case 2:
    result = 2*a2 + 6*a3*dt;break;
  case 3:
    result = 6*a3;break;
  default:
    result = 0;break;
  }
  assert( std::abs(result-derivative) < std::numeric_limits<double>::epsilon()*1e2 );
}

template <>
void poly_props<5>::check_derivative(Polynomial<5> const & poly, double t, int order, double derivative){
  double dt = t - poly.t0_;
  double dt2 = dt*dt;
  double dt3 = dt*dt2;
  double dt4 = dt*dt3;
  double a1 = poly.coefs_ [1];
  double a2 = poly.coefs_ [2];
  double a3 = poly.coefs_ [3];
  double a4 = poly.coefs_ [4];
  double a5 = poly.coefs_ [5];

  double result;
  switch (order) {
  case 0:
    result =  poly (t);break;
  case 1:
    result =  a1 + 2*a2*dt + 3*a3*dt2 + 4*a4*dt3 + 5*a5*dt4;break;
  case 2:
    result =  2*a2 + 6*a3*dt + 12*a4*dt2 + 20*a5*dt3;break;
  case 3:
    result =  6*a3 + 24*a4*dt + 60*a5*dt2;break;
  case 4:
    result =  24*a4 + 120*a5*dt;break;
  case 5:
    result =  120*a5;break;
  default:
    result =  0;break;
  }
  assert( std::abs(result-derivative) < std::numeric_limits<double>::epsilon()*1e2 );
}

/**
 * degree specific operation on Polynomial<3>
 */

template <>
void poly_props<3>::check_translate (Polynomial<3> const &poly, const double &t1, Polynomial<3> const &ref_poly)
{
  double dt = t1 - poly.t0_;
  double dt2 = dt*dt;
  double dt3 = dt*dt2;
  double a0 = poly.coefs_ [0];
  double a1 = poly.coefs_ [1];
  double a2 = poly.coefs_ [2];
  double a3 = poly.coefs_ [3];

  Eigen::Matrix<double,3+1,1> temp;
  temp [0] = a0 + a1*dt + a2*dt2 + a3*dt3;
  temp [1] = a1 + 2*dt*a2 + 3*dt2*a3;
  temp [2] = a2 + 3*dt*a3;
  temp [3] = a3;
  Polynomial<3> p_new(t1, temp);

  assert(p_new.t0_== ref_poly.t0_);
  for(int idx=0;idx<3+1;idx++){
    const double delta = std::abs(p_new.coefs_[idx]-ref_poly.coefs_[idx]);
    assert(delta < std::numeric_limits<double>::epsilon()*1e2);
  }
}

//FIXME: import code from Polynomial3 once it is replaced by Polynomial
template <>
void poly_props<3>::check_evaluate(double t0, Eigen::Matrix<double,3+1,1> &params, double t){
  Polynomial<3> new_poly(t0,params);
  Polynomial3 old_poly(t0,params(0),params(1),params(2),params(3));
  double delta = old_poly(t) - new_poly(t);
  assert(delta < std::numeric_limits<double>::epsilon()*1e3 );
}

/**
 * degree specific operations for Polynomial<5>
 */
/*template <>
  void poly_props<5>::check_translate (Polynomial<5> const &poly, const double &t1, Polynomial<5> const &ref_poly)
  {
  double dt = t1 - poly.t0_;
  double dt2 = dt*dt;
  double dt3 = dt*dt2;
  double dt4 = dt*dt3;
  double dt5 = dt*dt4;
  double a0 = poly.coefs_ [0];
  double a1 = poly.coefs_ [1];
  double a2 = poly.coefs_ [2];
  double a3 = poly.coefs_ [3];
  double a4 = poly.coefs_ [4];
  double a5 = poly.coefs_ [5];


  Eigen::Matrix<double,5+1,1> temp;
  //have doubts about this, minus sign is probably not correct
  temp [0] = a0 -   a1*dt +   a2*dt2 -   a3*dt3 +   a4*dt4 - a5*dt5;
  temp [1] = a1 - 2*a2*dt + 3*a3*dt2 - 4*a4*dt3 + 5*a5*dt4;
  temp [2] = a2 - 3*a3*dt + 6*a4*dt2 - 10*a5*dt3;
  temp [3] = a3 - 4*a4*dt + 10*a5*dt2;
  temp [4] = a4 - 5*a5*dt;
  temp [5] = a5;
  Polynomial<5> p_new(t1,temp);

  assert(p_new.t0_== ref_poly.t0_);
  for(int idx=0;idx<5+1;idx++){
  const double delta = std::abs(p_new.coefs_[idx]-ref_poly.coefs_[idx]);
  assert(delta < std::numeric_limits<double>::epsilon()*1e2);
  }
  }*/



template <int N>
void test_derivative(){
  Eigen::Matrix<double,N+1,1> params;
  params.setRandom();
  double t0 = (double)rand()/RAND_MAX;
  Polynomial<N> p_1(t0,params);

  double t = (double)rand()/RAND_MAX;

  for(int order=1;order<N+1;order++){
    double derivative = p_1.derivative(t,order);
    if(order>N)assert(derivative==0);

    poly_props<N>::check_derivative(p_1,t,order,derivative);	//only for N=3
  }
}

template <int N>
void test_translate(){
  Eigen::Matrix<double,N+1,1> params;
  params.setRandom();
  double t0 = (double)rand()/RAND_MAX;
  double t1 = (double)rand()/RAND_MAX;
  Polynomial<N> p_1(t0,params);
  p_1(0.);
  Polynomial<N> p_2 = p_1.translate(t1);
  p_2(0.);

  poly_props<N>::check_translate(p_1,t1,p_2);
}


template <int N>
void test_evaluate(){
  Eigen::Matrix<double,N+1,1> params;
  params.setRandom();
  double t0 = (double)rand()/RAND_MAX;
  double t = (double)rand()/RAND_MAX;
  poly_props<N>::check_evaluate(t0,params,t);
}

int main(){
  srand (time(NULL));
  test_evaluate<3>();
  //test_evaluate<5>();

  test_derivative<3>();
  test_derivative<5>();

  test_translate<3>();
  test_translate<5>();

  test_multiply<3>();
  test_multiply<5>();


}
