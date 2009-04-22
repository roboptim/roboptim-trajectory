#include <iostream>

#include "bspline.h"


bspline::bspline() : basis_function(4)  //number of spline dependence at each instant
{
}

bspline::bspline(int nb_fun_in, int nb_P_in, int nb_t_m_in, bool grad_of_time_in) : basis_function(4)
{
	define_sizes(nb_fun_in, nb_P_in, nb_t_m_in, grad_of_time_in);
}

bspline::bspline(int nb_fun_in, int nb_P_in, int nb_t_m_in, bool grad_of_time_in, bool defined_lc_in, bool free_fun_lc_in) : basis_function(4)
{
	define_sizes(nb_fun_in, nb_P_in, nb_t_m_in, grad_of_time_in, defined_lc_in, free_fun_lc_in);
}

bspline::~bspline()
{
}

void bspline::def_Nint() //define Nint from nb_P
{
	Nint = nb_P - 3;
}

inline double bspline::f11(double tpow3, double tpow2, double t)
{
	return - tpow3 + 3*tpow2 - 3*t + 1;
}

inline double bspline::f12(double tpow3, double tpow2, double t)
{
	return 7./4.*tpow3 - 9./2.*tpow2 + 3*t;
}

inline double bspline::f13(double tpow3, double tpow2, double t)
{
	return - 11./12.*tpow3 + 3./2.*tpow2;
}

inline double bspline::f21(double tpow3, double tpow2, double t)
{
	return - 1./4.*tpow3 + 3./4.*tpow2 - 3./4.*t + 1./4.;
}

inline double bspline::f22(double tpow3, double tpow2, double t)
{
	return 7./12.*tpow3 - 5./4.*tpow2 + t/4. + 7./12.;
}

inline double bspline::f31(double tpow3, double tpow2, double t)
{
	return - 1./6.*tpow3 + tpow2/2. - t/2. + 1./6.;
}

inline double bspline::fi1(double tpow3, double tpow2, double t)
{
	return tpow3/6.;
}

inline double bspline::fi2(double tpow3, double tpow2, double t)
{
	return - tpow3/2. + tpow2/2. + t/2. + 1./6.;
}

inline double bspline::fi3(double tpow3, double tpow2, double t)
{
	return tpow3/2. - tpow2 + 2./3.;
}

inline double bspline::fi4(double tpow3, double tpow2, double t)
{
	return - tpow3/6. + tpow2/2. - t/2. + 1./6.;
}


void bspline::calc_fun(boost::numeric::ublas::vector<double> *all_t, boost::numeric::ublas::matrix<double> *fun, int nb_t)
{
	using namespace boost::numeric::ublas;
	
	if (nb_fun!=(*fun).size1()) {
		std::cout << "error in bspline, fun does not have nb_fun rows" << std::endl;
	}
	if (nb_t> (int) (*fun).size2()) {
		std::cout << "error in bspline, fun does not have nb_t columns" << std::endl;
	}	
	if (defined_lc && !free_fun_lc && !given_lc) {
		std::cout << "error in bspline::calc_fun, computation of fun values with fixed limit conditions not defined" << std::endl;
	}

	int j;
    double tpw3, tpw2, t, tpwend3, tpwend2, tend;

    //computation of vector N such that N[j] is the number of time sample in [0,j*Tint]
    calc_N(all_t, nb_t);
	
	//Normalisation of all_t
	(*all_tn) = (*all_t) / Tint;
	
	//Effect of first base function
	//i indice for actionnors
	//j indice of time
	//k indice of phases limits
	
	
	for (j=N[0] ; j<N[1] ; j++)
    { 
        t = (*all_tn)[j];
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        column(*fun,j) = column(*matP,0)*f11(tpw3, tpw2, t) + column(*matP,1)*f12(tpw3, tpw2, t)
                                + column(*matP,2)*f13(tpw3, tpw2, t) + column(*matP,3)*fi1(tpw3, tpw2, t);        
    }
    
    for (j=N[1] ; j<N[2] ; j++)
    { 
        t = (*all_tn)[j] - 1;
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        column(*fun,j) = column(*matP,1)*f21(tpw3, tpw2, t) + column(*matP,2)*f22(tpw3, tpw2, t)
                                + column(*matP,3)*fi2(tpw3, tpw2, t) + column(*matP,4)*fi1(tpw3, tpw2, t);
    }
    
    for (j=N[2] ; j<N[3] ; j++)
    { 
        t = (*all_tn)[j] - 2;
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        column(*fun,j) = column(*matP,2)*f31(tpw3, tpw2, t) + column(*matP,3)*fi3(tpw3, tpw2, t)
                                + column(*matP,4)*fi2(tpw3, tpw2, t) + column(*matP,5)*fi1(tpw3, tpw2, t);
    }
    
	
	for (int k=3 ; k<Nint-3 ; k++)
	{
		for (j=N[k] ; j<N[k+1] ; j++)
        { 
            t = (*all_tn)[j] - k;
            tpw2 = pow(t,2);
            tpw3 = pow(t,3);
            column(*fun,j) = column(*matP,k)*fi4(tpw3, tpw2, t) + column(*matP,k+1)*fi3(tpw3, tpw2, t)
                                    + column(*matP,k+2)*fi2(tpw3, tpw2, t) + column(*matP,k+3)*fi1(tpw3, tpw2, t);
        }
	}
	
	for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-3);
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        tpwend3 = pow(tend,3);
        column(*fun,j) = column(*matP,Nint-3)*fi4(tpw3, tpw2, t) + column(*matP,Nint-2)*fi3(tpw3, tpw2, t)
                                + column(*matP,Nint-1)*fi2(tpw3, tpw2, t) 
                                + column(*matP,Nint)*f31(tpwend3, tpwend2, tend);
        
    }
    for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-2);
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        tpwend3 = pow(tend,3);
        column(*fun,j) = column(*matP,Nint-2)*fi4(tpw3, tpw2, t) + column(*matP,Nint-1)*fi3(tpw3, tpw2, t)
                                + column(*matP,Nint)*f22(tpwend3, tpwend2, tend) 
                                + column(*matP,Nint+1)*f21(tpwend3, tpwend2, tend);
        
    }
    for (j=N[Nint-1] ; j<N[Nint] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-1);
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        tpwend3 = pow(tend,3);
        column(*fun,j) = column(*matP,Nint-1)*fi4(tpw3, tpw2, t) 
                                + column(*matP,Nint)*f13(tpwend3, tpwend2, tend)
                                + column(*matP,Nint+1)*f12(tpwend3, tpwend2, tend) 
                                + column(*matP,Nint+2)*f11(tpwend3, tpwend2, tend);
        
    }
	
	
}


void bspline::calc_fun(double t, boost::numeric::ublas::vector<double> *fun)
{

	if (nb_fun!=(*fun).size()) {
		std::cout << "error in bspline, fun does not have nb_fun elements" << std::endl;
	}
	if (defined_lc && !free_fun_lc && !given_lc) {
		std::cout << "error in bspline::calc_fun, computation of fun values with fixed limit conditions not defined" << std::endl;
	}

	int j;
    double tpw3, tpw2, tpwend3, tpwend2, tend;

    //computation of vector N such that N[j] is the number of time sample in [0,j*Tint]
    int k = intervalle_of_time(t);
	
	//Normalisation of all_t
	t /= Tint;
	
	//Effect of first base function
	//i indice for actionnors
	//j indice of time
	
	
	if (k==0)
    { 
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        for (j=0; j<nb_fun; j++) {
            (*fun)(j) = (*matP)(j,0)*f11(tpw3, tpw2, t) + (*matP)(j,1)*f12(tpw3, tpw2, t)
                                    + (*matP)(j,2)*f13(tpw3, tpw2, t) + (*matP)(j,3)*fi1(tpw3, tpw2, t);        
        }
    } else if (k==1)
    { 
        t -= 1.;
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        for (j=0; j<nb_fun; j++) {
            (*fun)(j) = (*matP)(j,1)*f21(tpw3, tpw2, t) + (*matP)(j,2)*f22(tpw3, tpw2, t)
                                    + (*matP)(j,3)*fi2(tpw3, tpw2, t) + (*matP)(j,4)*fi1(tpw3, tpw2, t);
        }
    } else if (k==2)
    { 
        t -= 2.;
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        for (j=0; j<nb_fun; j++) {
            (*fun)(j) = (*matP)(j,2)*f31(tpw3, tpw2, t) + (*matP)(j,3)*fi3(tpw3, tpw2, t)
                                    + (*matP)(j,4)*fi2(tpw3, tpw2, t) + (*matP)(j,5)*fi1(tpw3, tpw2, t);
        }
    } else if (k<Nint-3)
	{
        t -= (double)(k);
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        for (j=0; j<nb_fun; j++) {
            (*fun)(j) = (*matP)(j,k)*fi4(tpw3, tpw2, t) + (*matP)(j,k+1)*fi3(tpw3, tpw2, t)
                                    + (*matP)(j,k+2)*fi2(tpw3, tpw2, t) + (*matP)(j,k+3)*fi1(tpw3, tpw2, t);
        }
	} else if (k==Nint-3)
    { 
        t -= (Nint-3);
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        tpwend3 = pow(tend,3);
        for (j=0; j<nb_fun; j++) {
            (*fun)(j) = (*matP)(j,Nint-3)*fi4(tpw3, tpw2, t) + (*matP)(j,Nint-2)*fi3(tpw3, tpw2, t)
                                    + (*matP)(j,Nint-1)*fi2(tpw3, tpw2, t) 
                                    + (*matP)(j,Nint)*f31(tpwend3, tpwend2, tend);
        }
    } else if (k==Nint-2)
    { 
        t -= (Nint-2);
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        tpwend3 = pow(tend,3);
        for (j=0; j<nb_fun; j++) {
            (*fun)(j) = (*matP)(j,Nint-2)*fi4(tpw3, tpw2, t) + (*matP)(j,Nint-1)*fi3(tpw3, tpw2, t)
                                    + (*matP)(j,Nint)*f22(tpwend3, tpwend2, tend) 
                                    + (*matP)(j,Nint+1)*f21(tpwend3, tpwend2, tend);
        }
    } else if (k==Nint-1)
    { 
        t -= (Nint-1);
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        tpwend3 = pow(tend,3);
        for (j=0; j<nb_fun; j++) {
            (*fun)(j) = (*matP)(j,Nint-1)*fi4(tpw3, tpw2, t) 
                                    + (*matP)(j,Nint)*f13(tpwend3, tpwend2, tend)
                                    + (*matP)(j,Nint+1)*f12(tpwend3, tpwend2, tend) 
                                    + (*matP)(j,Nint+2)*f11(tpwend3, tpwend2, tend);
        }
    }
	
	
}


void bspline::calc_fun_grad(boost::numeric::ublas::vector<double> *all_t, 
                       boost::numeric::ublas::matrix<double> *fun, 
                       boost::numeric::ublas::matrix< boost::numeric::ublas::vector<double> > *fun_grad, int nb_t)
{
	using namespace boost::numeric::ublas;
	
	if (nb_fun!=(*fun).size1()) {
		std::cout << "error in bspline, fun does not have nb_fun rows" << std::endl;
	}
	if (nb_t> (int) (*fun).size2()) {
		std::cout << "error in bspline, fun does not have nb_t columns" << std::endl;
	}	
	if (nb_fun!= (int)(*fun_grad).size1()) {
		std::cout << "error in bspline, fun_grad does not have nb_fun rows" << std::endl;
	}
	if (nb_t>(int) ((*fun_grad).size2())) {
		std::cout << "error in bspline, fun_grad does not have nb_t columns" << std::endl;
	}	
	if ((4!=(*fun_grad)(0,0).size()) && (!grad_of_time)) {
		std::cout << "error in bspline, fun_grad does not have 4 components of gradient" << std::endl;
	}	
	if ((5!=(*fun_grad)(0,0).size()) && (grad_of_time)) {
		std::cout << "error in bspline, fun_grad does not have 5 components of gradient" << std::endl;
	}	
	if (defined_lc && !free_fun_lc && !given_lc) {
		std::cout << "error in bspline::calc_fun_grad, computation of fun values with fixed limit conditions not defined" << std::endl;
	}

	int i, j;
    double tpw3, tpw2, t, tpwend3, tpwend2, tend;

    //computation of vector N such that N[j] is the number of time sample in [0,j*Tint]
    calc_N(all_t, nb_t);
	
	//Normalisation of all_t
	(*all_tn) = (*all_t) / Tint;
	
	//Effect of first base function
	//i indice for actionnors
	//j indice of time
	//k indice of phases limits
	
	
	for (j=N[0] ; j<N[1] ; j++)
    { 
        t = (*all_tn)[j];
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        f11v = f11(tpw3, tpw2, t);
        f12v = f12(tpw3, tpw2, t);
        f13v = f13(tpw3, tpw2, t);
        fi1v = fi1(tpw3, tpw2, t);
        column(*fun,j) = column(*matP,0)*f11v + column(*matP,1)*f12v
                                + column(*matP,2)*f13v + column(*matP,3)*fi1v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*fun_grad)(i,j)(0) = f11v;
                (*fun_grad)(i,j)(1) = f12v;
                (*fun_grad)(i,j)(2) = f13v;
                (*fun_grad)(i,j)(3) = fi1v;
                (*fun_grad)(i,j)(4) = 0;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*fun_grad)(i,j)(0) = f11v;
                (*fun_grad)(i,j)(1) = f12v;
                (*fun_grad)(i,j)(2) = f13v;
                (*fun_grad)(i,j)(3) = fi1v;
            }
        }
    }
    
    for (j=N[1] ; j<N[2] ; j++)
    { 
        t = (*all_tn)[j] - 1;
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        f21v = f21(tpw3, tpw2, t);
        f22v = f22(tpw3, tpw2, t);
        fi1v = fi1(tpw3, tpw2, t);
        fi2v = fi2(tpw3, tpw2, t);
        column(*fun,j) = column(*matP,1)*f21v + column(*matP,2)*f22v
                                + column(*matP,3)*fi2v + column(*matP,4)*fi1v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*fun_grad)(i,j)(0) = f21v;
                (*fun_grad)(i,j)(1) = f22v;
                (*fun_grad)(i,j)(2) = fi2v;
                (*fun_grad)(i,j)(3) = fi1v;
                (*fun_grad)(i,j)(4) = 0;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*fun_grad)(i,j)(0) = f21v;
                (*fun_grad)(i,j)(1) = f22v;
                (*fun_grad)(i,j)(2) = fi2v;
                (*fun_grad)(i,j)(3) = fi1v;
            }
        }
    }
    
    for (j=N[2] ; j<N[3] ; j++)
    { 
        t = (*all_tn)[j] - 2;
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        f31v = f31(tpw3, tpw2, t);
        fi1v = fi1(tpw3, tpw2, t);
        fi2v = fi2(tpw3, tpw2, t);
        fi3v = fi3(tpw3, tpw2, t);
        column(*fun,j) = column(*matP,2)*f31v + column(*matP,3)*fi3v
                                + column(*matP,4)*fi2v + column(*matP,5)*fi1v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*fun_grad)(i,j)(0) = f31v;
                (*fun_grad)(i,j)(1) = fi3v;
                (*fun_grad)(i,j)(2) = fi2v;
                (*fun_grad)(i,j)(3) = fi1v;
                (*fun_grad)(i,j)(4) = 0;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*fun_grad)(i,j)(0) = f31v;
                (*fun_grad)(i,j)(1) = fi3v;
                (*fun_grad)(i,j)(2) = fi2v;
                (*fun_grad)(i,j)(3) = fi1v;
            }
        }
    }
    
	
	for (int k=3 ; k<Nint-3 ; k++)
	{
		for (j=N[k] ; j<N[k+1] ; j++)
        { 
            t = (*all_tn)[j] - k;
            tpw2 = pow(t,2);
            tpw3 = pow(t,3);
            fi1v = fi1(tpw3, tpw2, t);
            fi2v = fi2(tpw3, tpw2, t);
            fi3v = fi3(tpw3, tpw2, t);
            fi4v = fi4(tpw3, tpw2, t);
            column(*fun,j) = column(*matP,k)*fi4v + column(*matP,k+1)*fi3v
                                    + column(*matP,k+2)*fi2v + column(*matP,k+3)*fi1v;
            if (grad_of_time) {
                for (i=0; i<nb_fun; i++) {
                    (*fun_grad)(i,j)(0) = fi4v;
                    (*fun_grad)(i,j)(1) = fi3v;
                    (*fun_grad)(i,j)(2) = fi2v;
                    (*fun_grad)(i,j)(3) = fi1v;
                    (*fun_grad)(i,j)(4) = 0;
                }
            } else {
                for (i=0; i<nb_fun; i++) {
                    (*fun_grad)(i,j)(0) = fi4v;
                    (*fun_grad)(i,j)(1) = fi3v;
                    (*fun_grad)(i,j)(2) = fi2v;
                    (*fun_grad)(i,j)(3) = fi1v;
                }
            }
        }
	}
	
	for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-3);
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        tpwend3 = pow(tend,3);
        fi4v = fi4(tpw3, tpw2, t);
        fi3v = fi3(tpw3, tpw2, t);
        fi2v = fi2(tpw3, tpw2, t);
        f31v = f31(tpwend3, tpwend2, tend);
        column(*fun,j) = column(*matP,Nint-3)*fi4v + column(*matP,Nint-2)*fi3v
                                + column(*matP,Nint-1)*fi2v 
                                + column(*matP,Nint)*f31v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*fun_grad)(i,j)(0) = fi4v;
                (*fun_grad)(i,j)(1) = fi3v;
                (*fun_grad)(i,j)(2) = fi2v;
                (*fun_grad)(i,j)(3) = f31v;
                (*fun_grad)(i,j)(4) = 0;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*fun_grad)(i,j)(0) = fi4v;
                (*fun_grad)(i,j)(1) = fi3v;
                (*fun_grad)(i,j)(2) = fi2v;
                (*fun_grad)(i,j)(3) = f31v;
            }
        }
    }

    for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-2);
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        tpwend3 = pow(tend,3);
        fi4v = fi4(tpw3, tpw2, t);
        fi3v = fi3(tpw3, tpw2, t);
        f21v = f21(tpwend3, tpwend2, tend);
        f22v = f22(tpwend3, tpwend2, tend);
        column(*fun,j) = column(*matP,Nint-2)*fi4v + column(*matP,Nint-1)*fi3v
                                + column(*matP,Nint)*f22v 
                                + column(*matP,Nint+1)*f21v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*fun_grad)(i,j)(0) = fi4v;
                (*fun_grad)(i,j)(1) = fi3v;
                (*fun_grad)(i,j)(2) = f22v;
                (*fun_grad)(i,j)(3) = f21v;
                (*fun_grad)(i,j)(4) = 0;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*fun_grad)(i,j)(0) = fi4v;
                (*fun_grad)(i,j)(1) = fi3v;
                (*fun_grad)(i,j)(2) = f22v;
                (*fun_grad)(i,j)(3) = f21v;
            }
        }
    }

    for (j=N[Nint-1] ; j<N[Nint] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-1);
        tpw2 = pow(t,2);
        tpw3 = pow(t,3);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        tpwend3 = pow(tend,3);
        fi4v = fi4(tpw3, tpw2, t);
        f11v = f11(tpwend3, tpwend2, tend);
        f12v = f12(tpwend3, tpwend2, tend);
        f13v = f13(tpwend3, tpwend2, tend);
        column(*fun,j) = column(*matP,Nint-1)*fi4v 
                                + column(*matP,Nint)*f13v
                                + column(*matP,Nint+1)*f12v 
                                + column(*matP,Nint+2)*f11v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*fun_grad)(i,j)(0) = fi4v;
                (*fun_grad)(i,j)(1) = f13v;
                (*fun_grad)(i,j)(2) = f12v;
                (*fun_grad)(i,j)(3) = f11v;
                (*fun_grad)(i,j)(4) = 0;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*fun_grad)(i,j)(0) = fi4v;
                (*fun_grad)(i,j)(1) = f13v;
                (*fun_grad)(i,j)(2) = f12v;
                (*fun_grad)(i,j)(3) = f11v;
            }
        }
    }
	
	
}

inline double bspline::df11(double tpow2, double t)
{
	return - 3*tpow2 + 6*t - 3;
}

inline double bspline::df12(double tpow2, double t)
{
	return 21./4.*tpow2 - 9.*t + 3;
}

inline double bspline::df13(double tpow2, double t)
{
	return - 11./4.*tpow2 + 3.*t;
}

inline double bspline::df21(double tpow2, double t)
{
	return - 3./4.*tpow2 + 3./2.*t - 3./4.;
}

inline double bspline::df22(double tpow2, double t)
{
	return 7./4.*tpow2 - 5./2.*t + 1./4.;
}

inline double bspline::df31(double tpow2, double t)
{
	return - 1./2.*tpow2 + t - 1./2.;
}

inline double bspline::dfi1(double tpow2, double t)
{
	return tpow2/2.;
}

inline double bspline::dfi2(double tpow2, double t)
{
	return - 3*tpow2/2. + t + 1./2.;
}

inline double bspline::dfi3(double tpow2, double t)
{
	return 3*tpow2/2. - 2*t;
}

inline double bspline::dfi4(double tpow2, double t)
{
	return - tpow2/2. + t - 1./2.;
}


void bspline::calc_dfun(boost::numeric::ublas::vector<double> *all_t, boost::numeric::ublas::matrix<double> *dfun, int nb_t)
{
	using namespace boost::numeric::ublas;
	
	if (nb_fun!=(int) ((*dfun).size1())) {
		std::cout << "error in bspline, dfun does not have nb_fun rows" << std::endl;
	}
	if (nb_t>(int) ((*dfun).size2())) {
		std::cout << "error in bspline, dfun does not have nb_t columns" << std::endl;
	}	
	if (defined_lc && !free_fun_lc && !given_lc) {
		std::cout << "error in bspline::calc_dfun, computation of fun values with fixed limit conditions not defined" << std::endl;
	}

	int j;
    double tpw2, t, tpwend2, tend;

    //computation of vector N such that N[j] is the number of time sample in [0,j*Tint]
    calc_N(all_t, nb_t);
	
	//Normalisation of all_t
	(*all_tn) = (*all_t) / Tint;
	
	//Effect of first base function
	//i indice for actionnors
	//j indice of time
	//k indice of phases limits
	
	
	for (j=N[0] ; j<N[1] ; j++)
    { 
        t = (*all_tn)[j];
        tpw2 = pow(t,2);
        column(*dfun,j) = (column(*matP,0)*(df11(tpw2, t)/Tint) + column(*matP,1)*(df12(tpw2, t)/Tint)
                                + column(*matP,2)*(df13(tpw2, t)/Tint) + column(*matP,3)*(dfi1(tpw2, t)/Tint));        
    }
    
    for (j=N[1] ; j<N[2] ; j++)
    { 
        t = (*all_tn)[j] - 1;
        tpw2 = pow(t,2);
        column(*dfun,j) = (column(*matP,1)*(df21(tpw2, t)/Tint) + column(*matP,2)*(df22(tpw2, t)/Tint)
                                + column(*matP,3)*(dfi2(tpw2, t)/Tint) + column(*matP,4)*(dfi1(tpw2, t)/Tint));
    }
    
    for (j=N[2] ; j<N[3] ; j++)
    { 
        t = (*all_tn)[j] - 2;
        tpw2 = pow(t,2);
        column(*dfun,j) = (column(*matP,2)*(df31(tpw2, t)/Tint) + column(*matP,3)*(dfi3(tpw2, t)/Tint)
                                + column(*matP,4)*(dfi2(tpw2, t)/Tint) + column(*matP,5)*(dfi1(tpw2, t)/Tint));
    }
    
	
	for (int k=3 ; k<Nint-3 ; k++)
	{
		for (j=N[k] ; j<N[k+1] ; j++)
        { 
            t = (*all_tn)[j] - k;
            tpw2 = pow(t,2);
            column(*dfun,j) = (column(*matP,k)*(dfi4(tpw2, t)/Tint) + column(*matP,k+1)*(dfi3(tpw2, t)/Tint)
                                    + column(*matP,k+2)*(dfi2(tpw2, t)/Tint) + column(*matP,k+3)*(dfi1(tpw2, t)/Tint));
        }
	}
	
	for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-3);
        tpw2 = pow(t,2);
        tend = -t + 1;
        tpwend2 = pow(tend,2);
        column(*dfun,j) = (column(*matP,Nint-3)*(dfi4(tpw2, t)/Tint) + column(*matP,Nint-2)*(dfi3(tpw2, t)/Tint)
                                + column(*matP,Nint-1)*(dfi2(tpw2, t)/Tint) 
                                - column(*matP,Nint)*(df31(tpwend2, tend)/Tint));
        
    }
    for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-2);
        tpw2 = pow(t,2);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        column(*dfun,j) = (column(*matP,Nint-2)*(dfi4(tpw2, t)/Tint) + column(*matP,Nint-1)*(dfi3(tpw2, t)/Tint)
                                - column(*matP,Nint)*(df22(tpwend2, tend)/Tint)
                                - column(*matP,Nint+1)*(df21(tpwend2, tend)/Tint));
        
    }
    for (j=N[Nint-1] ; j<N[Nint] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-1);
        tpw2 = pow(t,2);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        column(*dfun,j) = (column(*matP,Nint-1)*(dfi4(tpw2, t)/Tint)
                                - column(*matP,Nint)*(df13(tpwend2, tend)/Tint)
                                - column(*matP,Nint+1)*(df12(tpwend2, tend)/Tint)
                                - column(*matP,Nint+2)*(df11(tpwend2, tend)/Tint));
        
    }
}

void bspline::calc_dfun(double t, boost::numeric::ublas::vector<double> *dfun)
{

	if (nb_fun!=(*dfun).size()) {
		std::cout << "error in bspline, dfun does not have nb_fun elements" << std::endl;
	}
	if (defined_lc && !free_fun_lc && !given_lc) {
		std::cout << "error in bspline::calc_dfun, computation of fun values with fixed limit conditions not defined" << std::endl;
	}

	int j;
    double tpw2, tpwend2, tend;

    //computation of vector N such that N[j] is the number of time sample in [0,j*Tint]
    int k = intervalle_of_time(t);
	
	//Normalisation of all_t
	t /= Tint;
	
	//Effect of first base function
	//i indice for actionnors
	//j indice of time
	
	
	if (k==0)
    { 
        tpw2 = pow(t,2);
        for (j=0; j<nb_fun; j++) { //possibility of optimizing the computation
            (*dfun)(j) = (*matP)(j,0)*df11(tpw2, t)/Tint + (*matP)(j,1)*df12(tpw2, t)/Tint
                                    + (*matP)(j,2)*df13(tpw2, t)/Tint + (*matP)(j,3)*dfi1(tpw2, t)/Tint;        
        }
    } else if (k==1)
    { 
        t -= 1.;
        tpw2 = pow(t,2);
        for (j=0; j<nb_fun; j++) {
            (*dfun)(j) = (*matP)(j,1)*df21(tpw2, t)/Tint + (*matP)(j,2)*df22(tpw2, t)/Tint
                                    + (*matP)(j,3)*dfi2(tpw2, t)/Tint + (*matP)(j,4)*dfi1(tpw2, t)/Tint;
        }
    } else if (k==2)
    { 
        t -= 2.;
        tpw2 = pow(t,2);
        for (j=0; j<nb_fun; j++) {
            (*dfun)(j) = (*matP)(j,2)*df31(tpw2, t)/Tint + (*matP)(j,3)*dfi3(tpw2, t)/Tint
                                    + (*matP)(j,4)*dfi2(tpw2, t)/Tint + (*matP)(j,5)*dfi1(tpw2, t)/Tint;
        }
    } else if (k<Nint-3)
	{
        t -= (double)(k);
        tpw2 = pow(t,2);
        for (j=0; j<nb_fun; j++) {
            (*dfun)(j) = (*matP)(j,k)*dfi4(tpw2, t)/Tint + (*matP)(j,k+1)*dfi3(tpw2, t)/Tint
                                    + (*matP)(j,k+2)*dfi2(tpw2, t)/Tint + (*matP)(j,k+3)*dfi1(tpw2, t)/Tint;
        }
	} else if (k==Nint-3)
    { 
        t -= (Nint-3);
        tpw2 = pow(t,2);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        for (j=0; j<nb_fun; j++) {
            (*dfun)(j) = (*matP)(j,Nint-3)*dfi4(tpw2, t)/Tint + (*matP)(j,Nint-2)*dfi3(tpw2, t)/Tint
                                    + (*matP)(j,Nint-1)*dfi2(tpw2, t)/Tint 
                                    - (*matP)(j,Nint)*df31(tpwend2, tend)/Tint;
        }
    } else if (k==Nint-2)
    { 
        t -= (Nint-2);
        tpw2 = pow(t,2);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        for (j=0; j<nb_fun; j++) {
            (*dfun)(j) = (*matP)(j,Nint-2)*dfi4(tpw2, t)/Tint + (*matP)(j,Nint-1)*dfi3(tpw2, t)/Tint
                                    - (*matP)(j,Nint)*df22(tpwend2, tend)/Tint 
                                    - (*matP)(j,Nint+1)*df21(tpwend2, tend)/Tint;
        }
    } else if (k==Nint-1)
    { 
        t -= (Nint-1);
        tpw2 = pow(t,2);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        for (j=0; j<nb_fun; j++) {
            (*dfun)(j) = (*matP)(j,Nint-1)*dfi4(tpw2, t)/Tint 
                                    - (*matP)(j,Nint)*df13(tpwend2, tend)/Tint
                                    - (*matP)(j,Nint+1)*df12(tpwend2, tend)/Tint 
                                    - (*matP)(j,Nint+2)*df11(tpwend2, tend)/Tint;
        }
    }
	
	
}

void bspline::calc_dfun_grad(boost::numeric::ublas::vector<double> *all_t, 
                       boost::numeric::ublas::matrix<double> *dfun, 
                       boost::numeric::ublas::matrix< boost::numeric::ublas::vector<double> > *dfun_grad, int nb_t)
{
	using namespace boost::numeric::ublas;
	
	if (nb_fun!=(int) ((*dfun).size1())) {
		std::cout << "error in bspline, dfun does not have nb_fun rows" << std::endl;
	}
	if (nb_t> (int) ((*dfun).size2())) {
		std::cout << "error in bspline, dfun does not have nb_t columns" << std::endl;
	}	
	if (nb_fun!=(*dfun_grad).size1()) {
		std::cout << "error in bspline, dfun_grad does not have nb_fun rows" << std::endl;
	}
	if (nb_t>(int) (*dfun_grad).size2()) {
		std::cout << "error in bspline, dfun_grad does not have nb_t columns" << std::endl;
	}	
	if ((4!=(*dfun_grad)(0,0).size()) && (!grad_of_time)) {
		std::cout << "error in bspline, dfun_grad does not have 4 components of gradient" << std::endl;
	}	
	if ((5!=(*dfun_grad)(0,0).size()) && (grad_of_time)) {
		std::cout << "error in bspline, dfun_grad does not have 5 components of gradient" << std::endl;
	}	
	if (defined_lc && !free_fun_lc && !given_lc) {
		std::cout << "error in bspline::calc_dfun_grad, computation of fun values with fixed limit conditions not defined" << std::endl;
	}

	int i, j;
    double tpw2, t, tpwend2, tend;

    //computation of vector N such that N[j] is the number of time sample in [0,j*Tint]
    calc_N(all_t, nb_t);
	
	//Normalisation of all_t
	(*all_tn) = (*all_t) / Tint;
	
	//Effect of first base function
	//i indice for actionnors
	//j indice of time
	//k indice of phases limits
	
	
	for (j=N[0] ; j<N[1] ; j++)
    { 
        t = (*all_tn)[j];
        tpw2 = pow(t,2);
        df11v = df11(tpw2, t)/Tint;
        df12v = df12(tpw2, t)/Tint;
        df13v = df13(tpw2, t)/Tint;
        dfi1v = dfi1(tpw2, t)/Tint;
        column(*dfun,j) = column(*matP,0)*df11v + column(*matP,1)*df12v
                                + column(*matP,2)*df13v + column(*matP,3)*dfi1v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*dfun_grad)(i,j)(0) = df11v;
                (*dfun_grad)(i,j)(1) = df12v;
                (*dfun_grad)(i,j)(2) = df13v;
                (*dfun_grad)(i,j)(3) = dfi1v;
                (*dfun_grad)(i,j)(4) = -(*dfun)(i,j)/T;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*dfun_grad)(i,j)(0) = df11v;
                (*dfun_grad)(i,j)(1) = df12v;
                (*dfun_grad)(i,j)(2) = df13v;
                (*dfun_grad)(i,j)(3) = dfi1v;
            }
        }
    }
    
    for (j=N[1] ; j<N[2] ; j++)
    { 
        t = (*all_tn)[j] - 1;
        tpw2 = pow(t,2);
        df21v = df21(tpw2, t)/Tint;
        df22v = df22(tpw2, t)/Tint;
        dfi1v = dfi1(tpw2, t)/Tint;
        dfi2v = dfi2(tpw2, t)/Tint;
        column(*dfun,j) = column(*matP,1)*df21v + column(*matP,2)*df22v
                                + column(*matP,3)*dfi2v + column(*matP,4)*dfi1v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*dfun_grad)(i,j)(0) = df21v;
                (*dfun_grad)(i,j)(1) = df22v;
                (*dfun_grad)(i,j)(2) = dfi2v;
                (*dfun_grad)(i,j)(3) = dfi1v;
                (*dfun_grad)(i,j)(4) = -(*dfun)(i,j)/T;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*dfun_grad)(i,j)(0) = df21v;
                (*dfun_grad)(i,j)(1) = df22v;
                (*dfun_grad)(i,j)(2) = dfi2v;
                (*dfun_grad)(i,j)(3) = dfi1v;
            }
        }
    }
    
    for (j=N[2] ; j<N[3] ; j++)
    { 
        t = (*all_tn)[j] - 2;
        tpw2 = pow(t,2);
        df31v = df31(tpw2, t)/Tint;
        dfi1v = dfi1(tpw2, t)/Tint;
        dfi2v = dfi2(tpw2, t)/Tint;
        dfi3v = dfi3(tpw2, t)/Tint;
        column(*dfun,j) = column(*matP,2)*df31v + column(*matP,3)*dfi3v
                                + column(*matP,4)*dfi2v + column(*matP,5)*dfi1v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*dfun_grad)(i,j)(0) = df31v;
                (*dfun_grad)(i,j)(1) = dfi3v;
                (*dfun_grad)(i,j)(2) = dfi2v;
                (*dfun_grad)(i,j)(3) = dfi1v;
                (*dfun_grad)(i,j)(4) = -(*dfun)(i,j)/T;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*dfun_grad)(i,j)(0) = df31v;
                (*dfun_grad)(i,j)(1) = dfi3v;
                (*dfun_grad)(i,j)(2) = dfi2v;
                (*dfun_grad)(i,j)(3) = dfi1v;
            }
        }
    }
    
	
	for (int k=3 ; k<Nint-3 ; k++)
	{
        for (j=N[k] ; j<N[k+1] ; j++)
        { 
            t = (*all_tn)[j] - k;
            tpw2 = pow(t,2);
            dfi1v = dfi1(tpw2, t)/Tint;
            dfi2v = dfi2(tpw2, t)/Tint;
            dfi3v = dfi3(tpw2, t)/Tint;
            dfi4v = dfi4(tpw2, t)/Tint;
            column(*dfun,j) = column(*matP,k)*dfi4v + column(*matP,k+1)*dfi3v
                                 + column(*matP,k+2)*dfi2v + column(*matP,k+3)*dfi1v;
            if (grad_of_time) {
                for (i=0; i<nb_fun; i++) {
                    (*dfun_grad)(i,j)(0) = dfi4v;
                    (*dfun_grad)(i,j)(1) = dfi3v;
                    (*dfun_grad)(i,j)(2) = dfi2v;
                    (*dfun_grad)(i,j)(3) = dfi1v;
                    (*dfun_grad)(i,j)(4) = -(*dfun)(i,j)/T;
                }
            } else {
                for (i=0; i<nb_fun; i++) {
                    (*dfun_grad)(i,j)(0) = dfi4v;
                    (*dfun_grad)(i,j)(1) = dfi3v;
                    (*dfun_grad)(i,j)(2) = dfi2v;
                    (*dfun_grad)(i,j)(3) = dfi1v;
                }
            }
        }
	}
	
	for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-3);
        tpw2 = pow(t,2);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        dfi4v = dfi4(tpw2, t)/Tint;
        dfi3v = dfi3(tpw2, t)/Tint;
        dfi2v = dfi2(tpw2, t)/Tint;
        df31v = -df31(tpwend2, tend)/Tint;
        column(*dfun,j) = column(*matP,Nint-3)*dfi4v + column(*matP,Nint-2)*dfi3v
                                + column(*matP,Nint-1)*dfi2v 
                                + column(*matP,Nint)*df31v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*dfun_grad)(i,j)(0) = dfi4v;
                (*dfun_grad)(i,j)(1) = dfi3v;
                (*dfun_grad)(i,j)(2) = dfi2v;
                (*dfun_grad)(i,j)(3) = df31v;
                (*dfun_grad)(i,j)(4) = -(*dfun)(i,j)/T;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*dfun_grad)(i,j)(0) = dfi4v;
                (*dfun_grad)(i,j)(1) = dfi3v;
                (*dfun_grad)(i,j)(2) = dfi2v;
                (*dfun_grad)(i,j)(3) = df31v;
            }
        }
    }

    for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-2);
        tpw2 = pow(t,2);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        dfi4v = dfi4(tpw2, t)/Tint;
        dfi3v = dfi3(tpw2, t)/Tint;
        df21v = -df21(tpwend2, tend)/Tint;
        df22v = -df22(tpwend2, tend)/Tint;
        column(*dfun,j) = column(*matP,Nint-2)*dfi4v + column(*matP,Nint-1)*dfi3v
                                + column(*matP,Nint)*df22v 
                                + column(*matP,Nint+1)*df21v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*dfun_grad)(i,j)(0) = dfi4v;
                (*dfun_grad)(i,j)(1) = dfi3v;
                (*dfun_grad)(i,j)(2) = df22v;
                (*dfun_grad)(i,j)(3) = df21v;
                (*dfun_grad)(i,j)(4) = -(*dfun)(i,j)/T;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*dfun_grad)(i,j)(0) = dfi4v;
                (*dfun_grad)(i,j)(1) = dfi3v;
                (*dfun_grad)(i,j)(2) = df22v;
                (*dfun_grad)(i,j)(3) = df21v;
            }
        }
    }

    for (j=N[Nint-1] ; j<N[Nint] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-1);
        tpw2 = pow(t,2);
        tend = - t + 1;
        tpwend2 = pow(tend,2);
        dfi4v = dfi4(tpw2, t)/Tint;
        df11v = -df11(tpwend2, tend)/Tint;
        df12v = -df12(tpwend2, tend)/Tint;
        df13v = -df13(tpwend2, tend)/Tint;
        column(*dfun,j) = column(*matP,Nint-1)*dfi4v 
                                + column(*matP,Nint)*df13v
                                + column(*matP,Nint+1)*df12v 
                                + column(*matP,Nint+2)*df11v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*dfun_grad)(i,j)(0) = dfi4v;
                (*dfun_grad)(i,j)(1) = df13v;
                (*dfun_grad)(i,j)(2) = df12v;
                (*dfun_grad)(i,j)(3) = df11v;
                (*dfun_grad)(i,j)(4) = -(*dfun)(i,j)/T;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*dfun_grad)(i,j)(0) = dfi4v;
                (*dfun_grad)(i,j)(1) = df13v;
                (*dfun_grad)(i,j)(2) = df12v;
                (*dfun_grad)(i,j)(3) = df11v;
            }
        }
    }
	
	
}

inline double bspline::d2f11(double t)
{
	return - 6*t + 6;
}

inline double bspline::d2f12(double t)
{
	return 21./2.*t - 9.;
}

inline double bspline::d2f13(double t)
{
	return - 11./2.*t + 3.;
}

inline double bspline::d2f21(double t)
{
	return - 3./2.*t + 3./2.;
}

inline double bspline::d2f22(double t)
{
	return 7./2.*t - 5./2.;
}

inline double bspline::d2f31(double t)
{
	return - t + 1.;
}

inline double bspline::d2fi1(double t)
{
	return t;
}

inline double bspline::d2fi2(double t)
{
	return - 3*t + 1.;
}

inline double bspline::d2fi3(double t)
{
	return 3*t - 2.;
}

inline double bspline::d2fi4(double t)
{
	return - t + 1.;
}

void bspline::calc_ddfun(boost::numeric::ublas::vector<double> *all_t, boost::numeric::ublas::matrix<double> *ddfun, int nb_t)
{
	using namespace boost::numeric::ublas;
	
	if (nb_fun!=(*ddfun).size1()) {
		std::cout << "error in bspline, ddfun does not have nb_fun rows" << std::endl;
	}
	if (nb_t>(int) ((*ddfun).size2())) {
		std::cout << "error in bspline, ddfun does not have nb_t columns" << std::endl;
	}	
	if (defined_lc && !free_fun_lc && !given_lc) {
		std::cout << "error in bspline::calc_ddfun, computation of fun values with fixed limit conditions not defined" << std::endl;
	}

	int j;
    double t, tend;
    double Tint2 = pow(Tint,2);

    //computation of vector N such that N[j] is the number of time sample in [0,j*Tint]
    calc_N(all_t, nb_t);
	
	//Normalisation of all_t
	(*all_tn) = (*all_t) / Tint;
	
	//Effect of first base function
	//i indice for actionnors
	//j indice of time
	//k indice of phases limits
	
	
	for (j=N[0] ; j<N[1] ; j++)
    { 
        t = (*all_tn)[j];
        column(*ddfun,j) = (column(*matP,0)*(d2f11(t)/Tint2) + column(*matP,1)*(d2f12(t)/Tint2)
                                + column(*matP,2)*(d2f13(t)/Tint2) + column(*matP,3)*(d2fi1(t)/Tint2));        
    }
    
    for (j=N[1] ; j<N[2] ; j++)
    { 
        t = (*all_tn)[j] - 1;
        column(*ddfun,j) = (column(*matP,1)*(d2f21(t)/Tint2) + column(*matP,2)*(d2f22(t)/Tint2)
                                + column(*matP,3)*(d2fi2(t)/Tint2) + column(*matP,4)*(d2fi1(t)/Tint2));
    }
    
    for (j=N[2] ; j<N[3] ; j++)
    { 
        t = (*all_tn)[j] - 2;
        column(*ddfun,j) = (column(*matP,2)*(d2f31(t)/Tint2) + column(*matP,3)*(d2fi3(t)/Tint2)
                                + column(*matP,4)*(d2fi2(t)/Tint2) + column(*matP,5)*(d2fi1(t)/Tint2));
    }
    
	
	for (int k=3 ; k<Nint-3 ; k++)
	{
		for (j=N[k] ; j<N[k+1] ; j++)
        { 
            t = (*all_tn)[j] - k;
            column(*ddfun,j) = (column(*matP,k)*(d2fi4(t)/Tint2) + column(*matP,k+1)*(d2fi3(t)/Tint2)
                                    + column(*matP,k+2)*(d2fi2(t)/Tint2) + column(*matP,k+3)*(d2fi1(t)/Tint2));
        }
	}
	
	for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-3);
        tend = -t + 1;
        column(*ddfun,j) = (column(*matP,Nint-3)*(d2fi4(t)/Tint2) + column(*matP,Nint-2)*(d2fi3(t)/Tint2)
                                + column(*matP,Nint-1)*(d2fi2(t)/Tint2)
                                + column(*matP,Nint)*(d2f31(tend)/Tint2));
        
    }
    for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-2);
        tend = - t + 1;
        column(*ddfun,j) = (column(*matP,Nint-2)*(d2fi4(t)/Tint2) + column(*matP,Nint-1)*(d2fi3(t)/Tint2)
                                + column(*matP,Nint)*(d2f22(tend)/Tint2)
                                + column(*matP,Nint+1)*(d2f21(tend)/Tint2));
        
    }
    for (j=N[Nint-1] ; j<N[Nint] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-1);
        tend = - t + 1;
        column(*ddfun,j) = (column(*matP,Nint-1)*(d2fi4(t)/Tint2)
                                + column(*matP,Nint)*(d2f13(tend)/Tint2)
                                + column(*matP,Nint+1)*(d2f12(tend)/Tint2)
                                + column(*matP,Nint+2)*(d2f11(tend)/Tint2));
        
    }
}

void bspline::calc_ddfun(double t, boost::numeric::ublas::vector<double> *ddfun)
{

	if (nb_fun!=(*ddfun).size()) {
		std::cout << "error in bspline, ddfun does not have nb_fun elements" << std::endl;
	}
	if (defined_lc && !free_fun_lc && !given_lc) {
		std::cout << "error in bspline::calc_ddfun, computation of fun values with fixed limit conditions not defined" << std::endl;
	}

	int j;
    double tend;
    double Tint2 = pow(Tint,2);

    //computation of vector N such that N[j] is the number of time sample in [0,j*Tint]
    int k = intervalle_of_time(t);
	
	//Normalisation of all_t
	t /= Tint;
	
	//Effect of first base function
	//i indice for actionnors
	//j indice of time
	
	
	if (k==0)
    { 
        for (j=0; j<nb_fun; j++) { //possibility of optimizing the computation
            (*ddfun)(j) = (*matP)(j,0)*d2f11(t)/Tint2 + (*matP)(j,1)*d2f12(t)/Tint2
                                    + (*matP)(j,2)*d2f13(t)/Tint2 + (*matP)(j,3)*d2fi1(t)/Tint2;        
        }
    } else if (k==1)
    { 
        t -= 1.;
        for (j=0; j<nb_fun; j++) {
            (*ddfun)(j) = (*matP)(j,1)*d2f21(t)/Tint2 + (*matP)(j,2)*d2f22(t)/Tint2
                                    + (*matP)(j,3)*d2fi2(t)/Tint2 + (*matP)(j,4)*d2fi1(t)/Tint2;
        }
    } else if (k==2)
    { 
        t -= 2.;
        for (j=0; j<nb_fun; j++) {
            (*ddfun)(j) = (*matP)(j,2)*d2f31(t)/Tint2 + (*matP)(j,3)*d2fi3(t)/Tint2
                                    + (*matP)(j,4)*d2fi2(t)/Tint2 + (*matP)(j,5)*d2fi1(t)/Tint2;
        }
    } else if (k<Nint-3)
	{
        t -= (double)(k);
        for (j=0; j<nb_fun; j++) {
            (*ddfun)(j) = (*matP)(j,k)*d2fi4(t)/Tint2 + (*matP)(j,k+1)*d2fi3(t)/Tint2
                                    + (*matP)(j,k+2)*d2fi2(t)/Tint2 + (*matP)(j,k+3)*d2fi1(t)/Tint2;
        }
	} else if (k==Nint-3)
    { 
        t -= (Nint-3);
        tend = - t + 1;
        for (j=0; j<nb_fun; j++) {
            (*ddfun)(j) = (*matP)(j,Nint-3)*d2fi4(t)/Tint2 + (*matP)(j,Nint-2)*d2fi3(t)/Tint2
                                    + (*matP)(j,Nint-1)*d2fi2(t)/Tint2 
                                    + (*matP)(j,Nint)*d2f31(tend)/Tint2;
        }
    } else if (k==Nint-2)
    { 
        t -= (Nint-2);
        tend = - t + 1;
        for (j=0; j<nb_fun; j++) {
            (*ddfun)(j) = (*matP)(j,Nint-2)*d2fi4(t)/Tint2 + (*matP)(j,Nint-1)*d2fi3(t)/Tint2
                                    + (*matP)(j,Nint)*d2f22(tend)/Tint2 
                                    + (*matP)(j,Nint+1)*d2f21(tend)/Tint2;
        }
    } else if (k==Nint-1)
    { 
        t -= (Nint-1);
        tend = - t + 1;
        for (j=0; j<nb_fun; j++) {
            (*ddfun)(j) = (*matP)(j,Nint-1)*d2fi4(t)/Tint2 
                                    + (*matP)(j,Nint)*d2f13(tend)/Tint2
                                    + (*matP)(j,Nint+1)*d2f12(tend)/Tint2 
                                    + (*matP)(j,Nint+2)*d2f11(tend)/Tint2;
        }
    }
	
	
}

void bspline::calc_ddfun_grad(boost::numeric::ublas::vector<double> *all_t, 
                       boost::numeric::ublas::matrix<double> *ddfun, 
                       boost::numeric::ublas::matrix< boost::numeric::ublas::vector<double> > *ddfun_grad, int nb_t)
{
	using namespace boost::numeric::ublas;
	
	if (nb_fun!=(*ddfun).size1()) {
		std::cout << "error in bspline, ddfun does not have nb_fun rows" << std::endl;
	}
	if (nb_t>(int) ((*ddfun).size2())) {
		std::cout << "error in bspline, ddfun does not have nb_t columns" << std::endl;
	}	
	if (nb_fun!=(*ddfun_grad).size1()) {
		std::cout << "error in bspline, ddfun_grad does not have nb_fun rows" << std::endl;
	}
	if (nb_t> (int) ((*ddfun_grad).size2())) {
		std::cout << "error in bspline, ddfun_grad does not have nb_t columns" << std::endl;
	}	
	if ((4!=(*ddfun_grad)(0,0).size()) && (!grad_of_time)) {
		std::cout << "error in bspline, ddfun_grad does not have 4 components of gradient" << std::endl;
	}	
	if ((5!=(*ddfun_grad)(0,0).size()) && (grad_of_time)) {
		std::cout << "error in bspline, ddfun_grad does not have 5 components of gradient" << std::endl;
	}	
	if (defined_lc && !free_fun_lc && !given_lc) {
		std::cout << "error in bspline::calc_ddfun_grad, computation of fun values with fixed limit conditions not defined" << std::endl;
	}

	int i, j;
    double t, tend;
    double Tint2 = pow(Tint,2);

    //computation of vector N such that N[j] is the number of time sample in [0,j*Tint]
    calc_N(all_t, nb_t);
	
	//Normalisation of all_t
	(*all_tn) = (*all_t) / Tint;
	
	//Effect of first base function
	//i indice for actionnors
	//j indice of time
	//k indice of phases limits
	
	
	for (j=N[0] ; j<N[1] ; j++)
    { 
        t = (*all_tn)[j];
        d2f11v = d2f11(t)/Tint2;
        d2f12v = d2f12(t)/Tint2;
        d2f13v = d2f13(t)/Tint2;
        d2fi1v = d2fi1(t)/Tint2;
        column(*ddfun,j) = column(*matP,0)*d2f11v + column(*matP,1)*d2f12v
                                + column(*matP,2)*d2f13v + column(*matP,3)*d2fi1v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*ddfun_grad)(i,j)(0) = d2f11v;
                (*ddfun_grad)(i,j)(1) = d2f12v;
                (*ddfun_grad)(i,j)(2) = d2f13v;
                (*ddfun_grad)(i,j)(3) = d2fi1v;
                (*ddfun_grad)(i,j)(4) = -2.*(*ddfun)(i,j)/T;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*ddfun_grad)(i,j)(0) = d2f11v;
                (*ddfun_grad)(i,j)(1) = d2f12v;
                (*ddfun_grad)(i,j)(2) = d2f13v;
                (*ddfun_grad)(i,j)(3) = d2fi1v;
            }
        }
    }
    
    for (j=N[1] ; j<N[2] ; j++)
    { 
        t = (*all_tn)[j] - 1;
        d2f21v = d2f21(t)/Tint2;
        d2f22v = d2f22(t)/Tint2;
        d2fi1v = d2fi1(t)/Tint2;
        d2fi2v = d2fi2(t)/Tint2;
        column(*ddfun,j) = column(*matP,1)*d2f21v + column(*matP,2)*d2f22v
                                + column(*matP,3)*d2fi2v + column(*matP,4)*d2fi1v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*ddfun_grad)(i,j)(0) = d2f21v;
                (*ddfun_grad)(i,j)(1) = d2f22v;
                (*ddfun_grad)(i,j)(2) = d2fi2v;
                (*ddfun_grad)(i,j)(3) = d2fi1v;
                (*ddfun_grad)(i,j)(4) = -2.*(*ddfun)(i,j)/T;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*ddfun_grad)(i,j)(0) = d2f21v;
                (*ddfun_grad)(i,j)(1) = d2f22v;
                (*ddfun_grad)(i,j)(2) = d2fi2v;
                (*ddfun_grad)(i,j)(3) = d2fi1v;
            }
        }
    }
    
    for (j=N[2] ; j<N[3] ; j++)
    { 
        t = (*all_tn)[j] - 2;
        d2f31v = d2f31(t)/Tint2;
        d2fi1v = d2fi1(t)/Tint2;
        d2fi2v = d2fi2(t)/Tint2;
        d2fi3v = d2fi3(t)/Tint2;
        column(*ddfun,j) = column(*matP,2)*d2f31v + column(*matP,3)*d2fi3v
                                + column(*matP,4)*d2fi2v + column(*matP,5)*d2fi1v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*ddfun_grad)(i,j)(0) = d2f31v;
                (*ddfun_grad)(i,j)(1) = d2fi3v;
                (*ddfun_grad)(i,j)(2) = d2fi2v;
                (*ddfun_grad)(i,j)(3) = d2fi1v;
                (*ddfun_grad)(i,j)(4) = -2.*(*ddfun)(i,j)/T;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*ddfun_grad)(i,j)(0) = d2f31v;
                (*ddfun_grad)(i,j)(1) = d2fi3v;
                (*ddfun_grad)(i,j)(2) = d2fi2v;
                (*ddfun_grad)(i,j)(3) = d2fi1v;
            }
        }
    }
    
	
	for (int k=3 ; k<Nint-3 ; k++)
	{
		for (j=N[k] ; j<N[k+1] ; j++)
        { 
            t = (*all_tn)[j] - k;
            d2fi1v = d2fi1(t)/Tint2;
            d2fi2v = d2fi2(t)/Tint2;
            d2fi3v = d2fi3(t)/Tint2;
            d2fi4v = d2fi4(t)/Tint2;
            column(*ddfun,j) = column(*matP,k)*d2fi4v + column(*matP,k+1)*d2fi3v
                                    + column(*matP,k+2)*d2fi2v + column(*matP,k+3)*d2fi1v;
            if (grad_of_time) {
                for (i=0; i<nb_fun; i++) {
                    (*ddfun_grad)(i,j)(0) = d2fi4v;
                    (*ddfun_grad)(i,j)(1) = d2fi3v;
                    (*ddfun_grad)(i,j)(2) = d2fi2v;
                    (*ddfun_grad)(i,j)(3) = d2fi1v;
                    (*ddfun_grad)(i,j)(4) = -2.*(*ddfun)(i,j)/T;
                }
            } else {
                for (i=0; i<nb_fun; i++) {
                    (*ddfun_grad)(i,j)(0) = d2fi4v;
                    (*ddfun_grad)(i,j)(1) = d2fi3v;
                    (*ddfun_grad)(i,j)(2) = d2fi2v;
                    (*ddfun_grad)(i,j)(3) = d2fi1v;
                }
            }
        }
	}
	
	for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-3);
        tend = - t + 1;
        d2fi4v = d2fi4(t)/Tint2;
        d2fi3v = d2fi3(t)/Tint2;
        d2fi2v = d2fi2(t)/Tint2;
        d2f31v = d2f31(tend)/Tint2;
        column(*ddfun,j) = column(*matP,Nint-3)*d2fi4v + column(*matP,Nint-2)*d2fi3v
                                + column(*matP,Nint-1)*d2fi2v 
                                + column(*matP,Nint)*d2f31v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*ddfun_grad)(i,j)(0) = d2fi4v;
                (*ddfun_grad)(i,j)(1) = d2fi3v;
                (*ddfun_grad)(i,j)(2) = d2fi2v;
                (*ddfun_grad)(i,j)(3) = d2f31v;
                (*ddfun_grad)(i,j)(4) = -2.*(*ddfun)(i,j)/T;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*ddfun_grad)(i,j)(0) = d2fi4v;
                (*ddfun_grad)(i,j)(1) = d2fi3v;
                (*ddfun_grad)(i,j)(2) = d2fi2v;
                (*ddfun_grad)(i,j)(3) = d2f31v;
            }
        }
    }

    for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-2);
        tend = - t + 1;
        d2fi4v = d2fi4(t)/Tint2;
        d2fi3v = d2fi3(t)/Tint2;
        d2f21v = d2f21(tend)/Tint2;
        d2f22v = d2f22(tend)/Tint2;
        column(*ddfun,j) = column(*matP,Nint-2)*d2fi4v + column(*matP,Nint-1)*d2fi3v
                                + column(*matP,Nint)*d2f22v 
                                + column(*matP,Nint+1)*d2f21v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*ddfun_grad)(i,j)(0) = d2fi4v;
                (*ddfun_grad)(i,j)(1) = d2fi3v;
                (*ddfun_grad)(i,j)(2) = d2f22v;
                (*ddfun_grad)(i,j)(3) = d2f21v;
                (*ddfun_grad)(i,j)(4) = -2.*(*ddfun)(i,j)/T;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*ddfun_grad)(i,j)(0) = d2fi4v;
                (*ddfun_grad)(i,j)(1) = d2fi3v;
                (*ddfun_grad)(i,j)(2) = d2f22v;
                (*ddfun_grad)(i,j)(3) = d2f21v;
            }
        }
    }

    for (j=N[Nint-1] ; j<N[Nint] ; j++)
    { 
        t = (*all_tn)[j] - (Nint-1);
        tend = - t + 1;
        d2fi4v = d2fi4(t)/Tint2;
        d2f11v = d2f11(tend)/Tint2;
        d2f12v = d2f12(tend)/Tint2;
        d2f13v = d2f13(tend)/Tint2;
        column(*ddfun,j) = column(*matP,Nint-1)*d2fi4v 
                                + column(*matP,Nint)*d2f13v
                                + column(*matP,Nint+1)*d2f12v 
                                + column(*matP,Nint+2)*d2f11v;
        if (grad_of_time) {
            for (i=0; i<nb_fun; i++) {
                (*ddfun_grad)(i,j)(0) = d2fi4v;
                (*ddfun_grad)(i,j)(1) = d2f13v;
                (*ddfun_grad)(i,j)(2) = d2f12v;
                (*ddfun_grad)(i,j)(3) = d2f11v;
                (*ddfun_grad)(i,j)(4) = -2.*(*ddfun)(i,j)/T;
            }
        } else {
            for (i=0; i<nb_fun; i++) {
                (*ddfun_grad)(i,j)(0) = d2fi4v;
                (*ddfun_grad)(i,j)(1) = d2f13v;
                (*ddfun_grad)(i,j)(2) = d2f12v;
                (*ddfun_grad)(i,j)(3) = d2f11v;
            }
        }
    }
	
}



// give the position of the nonzero elements
void bspline::get_rank_nz_uncomp(boost::numeric::ublas::vector<int> &time_points_in,
                                 boost::numeric::ublas::vector<int> &funorcomp_type_in,
                                 boost::numeric::ublas::vector<int> &nb_nz_per_row_out,
                                 boost::numeric::ublas::matrix<int> &rank_grad_nz_out, int nb_t) 
{
	if (time_points_in.size()!=rank_grad_nz_out.size1()) {
		std::cout << "error in basis_function::get_rank_grad_u2uncomp, rank_grad_u2uncomp_in does not have as much rows as time_points_in" << std::endl;
	}
	if (time_points_in.size()!=nb_nz_per_row_out.size()) {
		std::cout << "error in basis_function::get_rank_grad_u2uncomp, nb_nz_per_row_out does not have as much elements as time_points_in" << std::endl;
	}
	if (nb_grad_comp!=rank_grad_nz_out.size2()) {
		std::cout << "error in basis_function::get_rank_grad_u2uncomp, rank_grad_u2uncomp_in does not have nb_grad_comp columns" << std::endl;
	}

    int i,kk;
    
    boost::numeric::ublas::vector<int> intervalles(time_points_in.size());
    for (i=0; i<(int) (time_points_in.size()); i++) {
    	if (time_points_in(i)==0) {
    		intervalles(i) = 0;
    	} else {
    		intervalles(i) = (time_points_in(i)*Nint-1)/(nb_t-1);
    	}
    	//std::cout << Nint << " " << intervalles(i) << std::endl;
    }
    
    if (!defined_lc) { // the limit conditions are not defined

        if (grad_of_time) {
            for (i=0; i<(int) (time_points_in.size()); i++) {
            	if (funorcomp_type_in(i)<0) {
                    for (kk=0; kk<nb_grad_comp-1; kk++) {
                        rank_grad_nz_out(i,kk) = kk+intervalles(i)*nb_fun;
                    }
                    rank_grad_nz_out(i,nb_grad_comp-1) = nb_grad_uncomp-1;
                    nb_nz_per_row_out(i) = nb_grad_comp;
            	} else {
                    for (kk=0; kk<nb_grad_fun-1; kk++) {
                        rank_grad_nz_out(i,kk) = funorcomp_type_in(i) + kk*nb_fun + intervalles(i)*nb_fun;
                    }
                    rank_grad_nz_out(i,nb_grad_comp-1) = nb_grad_uncomp-1;
                    nb_nz_per_row_out(i) = nb_grad_fun;
            	}
            }
        } else {
            for (i=0; i<(int) time_points_in.size(); i++) {
            	if (funorcomp_type_in(i)<0) {
                    for (kk=0; kk<nb_grad_comp; kk++) {
                        rank_grad_nz_out(i,kk) = kk+intervalles(i)*nb_fun;
                    }
                    nb_nz_per_row_out(i) = nb_grad_comp;
            	} else {
                    for (kk=0; kk<nb_grad_fun; kk++) {
                        rank_grad_nz_out(i,kk) = funorcomp_type_in(i) + kk*nb_fun + intervalles(i)*nb_fun;
                    }
                    nb_nz_per_row_out(i) = nb_grad_fun;
            	}
            }
        }
        
    //---------------------------------------------------------------------------
    } else if (defined_lc & !free_fun_lc) { // the limit conditions are completely defined
    
        if (grad_of_time) { // gradient of time is computed
            
            for (i=0; i<(int) time_points_in.size(); i++)
            {
                if (funorcomp_type_in(i)<0)
                {
					if (intervalles(i)<3) 
					{ 
						for (kk=0; kk<(intervalles(i)+1)*nb_fun; kk++) {
							rank_grad_nz_out(i,kk) = kk;
						}
						rank_grad_nz_out(i,(intervalles(i)+1)*nb_fun) = nb_grad_uncomp-1;
						nb_nz_per_row_out(i) = (intervalles(i)+1)*nb_fun + 1;
					}
					else if (intervalles(i)<Nint-3)
					{
						for (kk=0; kk<nb_grad_comp-1; kk++) {
							rank_grad_nz_out(i,kk) = kk+(intervalles(i)-3)*nb_fun;
						}
						rank_grad_nz_out(i,nb_grad_comp-1) = nb_grad_uncomp-1;
						nb_nz_per_row_out(i) = nb_grad_comp;
					}
					else 
					{
						for (kk=0; kk<(Nint-intervalles(i))*nb_fun; kk++) {
							rank_grad_nz_out(i,kk) = kk+(intervalles(i)-3)*nb_fun;
						}
						rank_grad_nz_out(i,(Nint-intervalles(i))*nb_fun) = nb_grad_uncomp-1;
						nb_nz_per_row_out(i) = (Nint-intervalles(i))*nb_fun + 1;
					}
				} else { //it is a gradient of a joint angle
                    if (intervalles(i)==0)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i);
                        rank_grad_nz_out(i,1) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 2;
                    }
                    else if (intervalles(i)==1)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i);
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+nb_fun;
                        rank_grad_nz_out(i,2) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 3;
                    }
                    else if (intervalles(i)==2)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i);
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+nb_fun;
                        rank_grad_nz_out(i,2) = funorcomp_type_in(i)+2*nb_fun;
                        rank_grad_nz_out(i,3) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 4;
                    }
                    else if (intervalles(i)<Nint-3)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i)+(intervalles(i)-3)*nb_fun;
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+(intervalles(i)-2)*nb_fun;
                        rank_grad_nz_out(i,2) = funorcomp_type_in(i)+(intervalles(i)-1)*nb_fun;
                        rank_grad_nz_out(i,3) = funorcomp_type_in(i)+ intervalles(i)   *nb_fun;
                        rank_grad_nz_out(i,4) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 5;
                    }
                    else if (intervalles(i)==Nint-3)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i)+(Nint-6)*nb_fun;
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+(Nint-5)*nb_fun;
                        rank_grad_nz_out(i,2) = funorcomp_type_in(i)+(Nint-4)*nb_fun;
                        rank_grad_nz_out(i,3) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 4;
                    }
                    else if (intervalles(i)==Nint-2) 
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i)+(Nint-5)*nb_fun;
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+(Nint-4)*nb_fun;
                        rank_grad_nz_out(i,2) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 3;
                    }
                    else if (intervalles(i)==Nint-1)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i)+(Nint-4)*nb_fun;
                        rank_grad_nz_out(i,1) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 2;
                    }
				}
            }
    //---------------------------------------------------------------------------
        } else { // gradient of time is not computed
            for (i=0; i<(int) time_points_in.size(); i++)
            {
                if (funorcomp_type_in(i)<0)
				{
					if (intervalles(i)<3) 
					{ 
						for (kk=0; kk<(intervalles(i)+1)*nb_fun; kk++) {
							rank_grad_nz_out(i,kk) = kk;
						}
						nb_nz_per_row_out(i) = (intervalles(i)+1)*nb_fun;
					}
					else if (intervalles(i)<Nint-3)
					{
						for (kk=0; kk<nb_grad_comp; kk++) {
							rank_grad_nz_out(i,kk) = kk+(intervalles(i)-3)*nb_fun;
						}
						nb_nz_per_row_out(i) = nb_grad_comp;
					}
					else 
					{
						for (kk=0; kk<(Nint-intervalles(i))*nb_fun; kk++) {
							rank_grad_nz_out(i,kk) = kk+(intervalles(i)-3)*nb_fun;
						}
						nb_nz_per_row_out(i) = (Nint-intervalles(i))*nb_fun;
					}
				} else { //it is a gradient of a joint angle
                    if (intervalles(i)==0)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i);
                        nb_nz_per_row_out(i) = 1;
                    }
                    else if (intervalles(i)==1)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i);
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+nb_fun;
                        nb_nz_per_row_out(i) = 2;
                    }
                    else if (intervalles(i)==2)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i);
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+nb_fun;
                        rank_grad_nz_out(i,2) = funorcomp_type_in(i)+2*nb_fun;
                        nb_nz_per_row_out(i) = 3;
                    }
                    else if (intervalles(i)<Nint-3)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i)+(intervalles(i)-3)*nb_fun;
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+(intervalles(i)-2)*nb_fun;
                        rank_grad_nz_out(i,2) = funorcomp_type_in(i)+(intervalles(i)-1)*nb_fun;
                        rank_grad_nz_out(i,3) = funorcomp_type_in(i)+ intervalles(i)   *nb_fun;
                        nb_nz_per_row_out(i) = 4;
                    }
                    else if (intervalles(i)==Nint-3)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i)+(Nint-6)*nb_fun;
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+(Nint-5)*nb_fun;
                        rank_grad_nz_out(i,2) = funorcomp_type_in(i)+(Nint-4)*nb_fun;
                        nb_nz_per_row_out(i) = 3;
                    }
                    else if (intervalles(i)==Nint-2) 
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i)+(Nint-5)*nb_fun;
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+(Nint-4)*nb_fun;
                        nb_nz_per_row_out(i) = 2;
                    }
                    else if (intervalles(i)==Nint-1)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i)+(Nint-4)*nb_fun;
                        nb_nz_per_row_out(i) = 1;
                    }
				}
            }
        } // end gradient of time is not computed

    //---------------------------------------------------------------------------
    } else { // the limit conditions are defined, limit conditions of fun are free
        
        if (grad_of_time) { // gradient of time is computed
            for (i=0; i<(int) time_points_in.size(); i++)
            {
                if (funorcomp_type_in(i)<0)
                {
                    if (intervalles(i)==0)
                    {
                        for (kk=0; kk<nb_fun; kk++) {
                            rank_grad_nz_out(i,kk) = kk;
                            rank_grad_nz_out(i,kk+nb_fun) = kk+nb_fun;
                        }
                        rank_grad_nz_out(i,2*nb_fun) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 2*nb_fun + 1;
                    }
                    else if (intervalles(i)==1)
                    {
                        for (kk=0; kk<nb_fun; kk++) {
                            rank_grad_nz_out(i,kk) = kk;
                        }
                        for (kk=0; kk<2*nb_fun; kk++) {
                            rank_grad_nz_out(i,kk+nb_fun) = kk+nb_fun;
                        }
                        rank_grad_nz_out(i,3*nb_fun) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 3*nb_fun + 1;
                    }
                    else if (intervalles(i)==2)
                    {
                        for (kk=0; kk<nb_fun; kk++) {
                            rank_grad_nz_out(i,kk) = kk;
                        }
                        for (kk=0; kk<3*nb_fun; kk++) {
                            rank_grad_nz_out(i,kk+nb_fun) = kk+nb_fun;
                        }
                        rank_grad_nz_out(i,nb_grad_comp-1) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = nb_grad_comp;
                    }
                    else if (intervalles(i)<Nint-3)
                    {
                        for (kk=0; kk<nb_grad_comp-1; kk++) {
                            rank_grad_nz_out(i,kk) = kk+(intervalles(i)-2)*nb_fun;
                        }
                        rank_grad_nz_out(i,nb_grad_comp-1) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = nb_grad_comp;
                    }
                    else if (intervalles(i)==Nint-3)
                    {
                        for (kk=0; kk<3*nb_fun; kk++) {
                            rank_grad_nz_out(i,kk) = kk+(Nint-5)*nb_fun;
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            rank_grad_nz_out(i,kk+3*nb_fun) = kk+(Nint-2)*nb_fun;
                        }
                        rank_grad_nz_out(i,nb_grad_comp-1) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = nb_grad_comp;
                    }
                    else if (intervalles(i)==Nint-2) 
                    {
                        for (kk=0; kk<2*nb_fun; kk++) {
                            rank_grad_nz_out(i,kk) = kk+(Nint-4)*nb_fun;
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            rank_grad_nz_out(i,kk+2*nb_fun) = kk+(Nint-2)*nb_fun;
                        }
                        rank_grad_nz_out(i,3*nb_fun) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 3*nb_fun + 1;
                    }
                    else if (intervalles(i)==Nint-1)
                    {
                        for (kk=0; kk<nb_fun; kk++) {
                            rank_grad_nz_out(i,kk) = kk+(Nint-3)*nb_fun;
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            rank_grad_nz_out(i,kk+nb_fun) = kk+(Nint-2)*nb_fun;
                        }
                        rank_grad_nz_out(i,2*nb_fun) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 2*nb_fun + 1;
                    }
                } else {
                    if (intervalles(i)==0)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i);
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+nb_fun;
                        rank_grad_nz_out(i,2) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 3;
                    }
                    else if (intervalles(i)==1)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i);
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+nb_fun;
                        rank_grad_nz_out(i,2) = funorcomp_type_in(i)+2*nb_fun;
                        rank_grad_nz_out(i,3) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 4;
                    }
                    else if (intervalles(i)==2)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i);
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+nb_fun;
                        rank_grad_nz_out(i,2) = funorcomp_type_in(i)+2*nb_fun;
                        rank_grad_nz_out(i,3) = funorcomp_type_in(i)+3*nb_fun;
                        rank_grad_nz_out(i,4) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 5;
                    }
                    else if (intervalles(i)<Nint-3)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i)+(intervalles(i)-2)*nb_fun;
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+(intervalles(i)-1)*nb_fun;
                        rank_grad_nz_out(i,2) = funorcomp_type_in(i)+ intervalles(i)   *nb_fun;
                        rank_grad_nz_out(i,3) = funorcomp_type_in(i)+(intervalles(i)+1)*nb_fun;
                        rank_grad_nz_out(i,4) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 5;
                    }
                    else if (intervalles(i)==Nint-3)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i)+(Nint-5)*nb_fun;
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+(Nint-4)*nb_fun;
                        rank_grad_nz_out(i,2) = funorcomp_type_in(i)+(Nint-3)*nb_fun;
                        rank_grad_nz_out(i,3) = funorcomp_type_in(i)+(Nint-2)*nb_fun;
                        rank_grad_nz_out(i,4) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 5;
                    }
                    else if (intervalles(i)==Nint-2) 
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i)+(Nint-4)*nb_fun;
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+(Nint-3)*nb_fun;
                        rank_grad_nz_out(i,2) = funorcomp_type_in(i)+(Nint-2)*nb_fun;
                        rank_grad_nz_out(i,3) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 4;
                    }
                    else if (intervalles(i)==Nint-1)
                    {
                        rank_grad_nz_out(i,0) = funorcomp_type_in(i)+(Nint-3)*nb_fun;
                        rank_grad_nz_out(i,1) = funorcomp_type_in(i)+(Nint-2)*nb_fun;
                        rank_grad_nz_out(i,2) = nb_grad_uncomp-1;
                        nb_nz_per_row_out(i) = 3;
                    }
                }
            }
    //---------------------------------------------------------------------------
        } else { // gradient of time is not computed
            for (i=0; i<(int) time_points_in.size(); i++)
            {
                if (intervalles(i)==0)
                {
                    for (kk=0; kk<nb_fun; kk++) {
                        rank_grad_nz_out(i,kk) = kk;
                        rank_grad_nz_out(i,kk+nb_fun) = kk+nb_fun;
                    }
                    nb_nz_per_row_out(i) = 2*nb_fun;
                }
                else if (intervalles(i)==1)
                {
                    for (kk=0; kk<nb_fun; kk++) {
                        rank_grad_nz_out(i,kk) = kk;
                    }
                    for (kk=0; kk<2*nb_fun; kk++) {
                        rank_grad_nz_out(i,kk+nb_fun) = kk+nb_fun;
                    }
                    nb_nz_per_row_out(i) = 3*nb_fun;
                }
                else if (intervalles(i)==2)
                {
                    for (kk=0; kk<nb_fun; kk++) {
                        rank_grad_nz_out(i,kk) = kk;
                    }
                    for (kk=0; kk<3*nb_fun; kk++) {
                        rank_grad_nz_out(i,kk+nb_fun) = kk+nb_fun;
                    }
                    nb_nz_per_row_out(i) = nb_grad_comp;
                }
                else if (intervalles(i)<Nint-3)
                {
                    for (kk=0; kk<nb_grad_comp-1; kk++) {
                        rank_grad_nz_out(i,kk) = kk+(intervalles(i)-2)*nb_fun;
                    }
                    nb_nz_per_row_out(i) = nb_grad_comp;
                }
                else if (intervalles(i)==Nint-3)
                {
                    for (kk=0; kk<3*nb_fun; kk++) {
                        rank_grad_nz_out(i,kk) = kk+(Nint-5)*nb_fun;
                    }
                    for (kk=0; kk<nb_fun; kk++) {
                        rank_grad_nz_out(i,kk+3*nb_fun) = kk+(Nint-2)*nb_fun;
                    }
                    nb_nz_per_row_out(i) = nb_grad_comp;
                }
                else if (intervalles(i)==Nint-2) 
                {
                    for (kk=0; kk<2*nb_fun; kk++) {
                        rank_grad_nz_out(i,kk) = kk+(Nint-4)*nb_fun;
                    }
                    for (kk=0; kk<nb_fun; kk++) {
                        rank_grad_nz_out(i,kk+2*nb_fun) = kk+(Nint-2)*nb_fun;
                    }
                    nb_nz_per_row_out(i) = 3*nb_fun;
                }
                else if (intervalles(i)==Nint-1)
                {
                    for (kk=0; kk<nb_fun; kk++) {
                        rank_grad_nz_out(i,kk) = kk+(Nint-3)*nb_fun;
                    }
                    for (kk=0; kk<nb_fun; kk++) {
                        rank_grad_nz_out(i,kk+nb_fun) = kk+(Nint-2)*nb_fun;
                    }
                    nb_nz_per_row_out(i) = 2*nb_fun;
                }
            }
        } // end gradient of time is not computed
    } // end limit conditions are defined

}
    
//convert the gradient with respect to all splines coefficients to the considered parameters
// and put all non zero values to the left
// and perform the scaling
void bspline::grad_comp2sparse(boost::numeric::ublas::vector<int> &time_points_in,
                               boost::numeric::ublas::vector<int> &funorcomp_type_in,
                               boost::numeric::ublas::matrix<double> &data_grad, int nb_t)
{
	if (time_points_in.size()!=data_grad.size1()) {
		std::cout << "error in bspline::grad_comp2sparse, data_grad does not have as much rows as there is time_points_in" << std::endl;
	}
	if (nb_grad_comp!=data_grad.size2()) {
		std::cout << "error in bspline::grad_comp2sparse, data_grad does not have nb_grad_comp columns" << std::endl;
	}

    int i,kk;
    
    boost::numeric::ublas::vector<int> intervalles(time_points_in.size());
    for (i=0; i<(int) time_points_in.size(); i++) {
    	if (time_points_in(i)==0) {
    		intervalles(i) = 0;
    	} else {
    		intervalles(i) = (time_points_in(i)*Nint-1)/(nb_t-1);
    	}
    	//std::cout << Nint << " " << intervalles(i) << std::endl;
    }
    

    if (!defined_lc) { // the limit conditions are not defined

        //nothing to be done
        
//---------------------------------------------------------------------------
    } else if (defined_lc & !free_fun_lc) { // the limit conditions are completely defined
    
        if (grad_of_time) { // gradient of time is computed
            
            for (i=0; i<(int) time_points_in.size(); i++)
            {
            	if (funorcomp_type_in(i)<0) {
					if (intervalles(i)<3) 
					{ 
						//few translation for the first three intervalls so that gradiant dependance are on the left
						// since there is no dependance with first spline coefficient (due to fixed position, volicity, acc)
						for (kk=0; kk<(intervalles(i)+1)*nb_fun; kk++) {
							data_grad(i,kk) = data_grad(i,kk+(3-intervalles(i))*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
						}
						data_grad(i,(intervalles(i)+1)*nb_fun) = data_grad(i,nb_grad_comp-1)/T_scale;
					}
					else if (intervalles(i)<Nint-3)
					{
						//perform only scaling
						for (kk=0; kk<4*nb_fun; kk++) {
							data_grad(i,kk) /= P_scale*funs_scale(kk%nb_fun);
						}
						data_grad(i,nb_grad_comp-1) /= T_scale;
					}
					else 
					{
						//perform scaling and translation of gradient with respect to duration
						for (kk=0; kk<(Nint-intervalles(i))*nb_fun; kk++) {
							data_grad(i,kk) /= P_scale*funs_scale(kk%nb_fun);
						}
						data_grad(i,(Nint-intervalles(i))*nb_fun) = data_grad(i,nb_grad_comp-1)/T_scale;
					}
				} else { //it is a gradient of a joint angle
            		double fun_scale = funs_scale(funorcomp_type_in(i));
					if (intervalles(i)<3) 
					{ 
						//few translation for the first three intervalls so that gradiant dependance are on the left
						// since there is no dependance with first spline coefficient (due to fixed position, volicity, acc)
						for (kk=0; kk<intervalles(i)+1; kk++) {
							data_grad(i,kk) = data_grad(i,kk+3-intervalles(i))/(P_scale*fun_scale);
						}
						data_grad(i,intervalles(i)+1) = data_grad(i,4)/T_scale;
					}
					else if (intervalles(i)<Nint-3)
					{
						//perform only scaling
						data_grad(i,0) /= P_scale*fun_scale;
						data_grad(i,1) /= P_scale*fun_scale;
						data_grad(i,2) /= P_scale*fun_scale;
						data_grad(i,3) /= P_scale*fun_scale;
						data_grad(i,4) /= T_scale;
					}
					else 
					{
						//perform scaling and translation of gradient with respect to duration
						for (kk=0; kk<Nint-intervalles(i); kk++) {
							data_grad(i,kk) /= P_scale*fun_scale;
						}
						data_grad(i,Nint-intervalles(i)) = data_grad(i,4)/T_scale;
					}
				}
            }
//---------------------------------------------------------------------------
        } else { // gradient of time is not computed
            for (i=0; i<(int) time_points_in.size(); i++)
            {
            	if (funorcomp_type_in(i)<0) {
					if (intervalles(i)<3) 
					{ 
						//few translation for the first three intervalls so that gradiant dependance are on the left
						// since there is no dependance with first spline coefficient (due to fixed position, volicity, acc)
						for (kk=0; kk<(intervalles(i)+1)*nb_fun; kk++) {
							data_grad(i,kk) = data_grad(i,kk+(3-intervalles(i))*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
						}
					}
					else if (intervalles(i)<Nint-3)
					{
						//perform only scaling
						for (kk=0; kk<4*nb_fun; kk++) {
							data_grad(i,kk) /= P_scale*funs_scale(kk%nb_fun);
						}
					}
					else 
					{
						//perform only scaling of necessary components
						for (kk=0; kk<(Nint-intervalles(i))*nb_fun; kk++) {
							data_grad(i,kk) /= P_scale*funs_scale(kk%nb_fun);
						}
					}
				} else { //it is a gradient of a joint angle
            		double fun_scale = funs_scale(funorcomp_type_in(i));
					if (intervalles(i)<3) 
					{ 
						//few translation for the first three intervalls so that gradiant dependance are on the left
						// since there is no dependance with first spline coefficient (due to fixed position, volicity, acc)
						for (kk=0; kk<intervalles(i)+1; kk++) {
							data_grad(i,kk) = data_grad(i,kk+3-intervalles(i))/(P_scale*fun_scale);
						}
					}
					else if (intervalles(i)<Nint-3)
					{
						//perform only scaling
						data_grad(i,0) /= P_scale*fun_scale;
						data_grad(i,1) /= P_scale*fun_scale;
						data_grad(i,2) /= P_scale*fun_scale;
						data_grad(i,3) /= P_scale*fun_scale;
					}
					else 
					{
						//perform scaling and translation of gradient with respect to duration
						for (kk=0; kk<Nint-intervalles(i); kk++) {
							data_grad(i,kk) /= P_scale*fun_scale;
						}
					}
				}
            }
        } // end gradient of time is not computed

//---------------------------------------------------------------------------
    } else { // the limit conditions are defined, limit conditions of fun are free
        calc_P_fun_0();
        calc_P_fun_f();
        
        if (grad_of_time) { // gradient of time is computed
            for (i=0; i<(int) time_points_in.size(); i++)
            { 
            	if (funorcomp_type_in(i)<0) {
                    if (intervalles(i)==0)
                    {
                        for (kk=0; kk<nb_fun; kk++) {
                            data_grad(i,kk) *= P0_fun/(funif_scale*funs_scale(kk%nb_fun));
                            data_grad(i,kk) += (data_grad(i,kk+nb_fun)*P1_fun
                                                  + data_grad(i,kk+2*nb_fun)*P2_fun)/(funif_scale*funs_scale(kk%nb_fun));
                            data_grad(i,kk+nb_fun) = data_grad(i,kk+3*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        data_grad(i,2*nb_fun) = data_grad(i,nb_grad_comp-1)/T_scale;
                    }
                    else if (intervalles(i)==1)
                    {
                        for (kk=0; kk<nb_fun; kk++) {
                            data_grad(i,kk) *= P1_fun/(funif_scale*funs_scale(kk%nb_fun));
                            data_grad(i,kk) += data_grad(i,kk+nb_fun)*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=0; kk<2*nb_fun; kk++) {
                            data_grad(i,kk+nb_fun) = data_grad(i,kk+2*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        data_grad(i,3*nb_fun) = data_grad(i,nb_grad_comp-1)/T_scale;
                    }
                    else if (intervalles(i)==2)
                    {
                        for (kk=0; kk<nb_fun; kk++) {
                            data_grad(i,kk) *= P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=0; kk<3*nb_fun; kk++) {
                            data_grad(i,kk+nb_fun) /= P_scale*funs_scale(kk%nb_fun);
                        }
                        data_grad(i,4*nb_fun) /= T_scale;
                    }
                    else if (intervalles(i)<Nint-3)
                    {
                        for (kk=0; kk<4*nb_fun; kk++) {
                            data_grad(i,kk) /= P_scale*funs_scale(kk%nb_fun);
                        }
                        data_grad(i,4*nb_fun) /= T_scale;
                    }
                    else if (intervalles(i)==Nint-3)
                    {
                        for (kk=0; kk<3*nb_fun; kk++) {
                            data_grad(i,kk) /= P_scale*funs_scale(kk%nb_fun);
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            data_grad(i,kk+3*nb_fun) = data_grad(i,kk+3*nb_fun)*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        }
                        data_grad(i,4*nb_fun) /= T_scale;
                    }
                    else if (intervalles(i)==Nint-2) 
                    {
                        for (kk=0; kk<2*nb_fun; kk++) {
                            data_grad(i,kk) /= P_scale*funs_scale(kk%nb_fun);
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            data_grad(i,kk+2*nb_fun) *= Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                            data_grad(i,kk+2*nb_fun) += data_grad(i,kk+3*nb_fun)*Pendm1_fun/(funif_scale*funs_scale(kk%nb_fun));
                        }
                        data_grad(i,3*nb_fun) = data_grad(i,nb_grad_comp-1)/T_scale;
                    }
                    else if (intervalles(i)==Nint-1)
                    {
                        for (kk=0; kk<nb_fun; kk++) {
                            data_grad(i,kk) /= P_scale*funs_scale(kk%nb_fun);
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            data_grad(i,kk+nb_fun) *= Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                            data_grad(i,kk+nb_fun) += (data_grad(i,kk+2*nb_fun)*Pendm1_fun
                                                            + data_grad(i,kk+3*nb_fun)*Pend_fun)/(funif_scale*funs_scale(kk%nb_fun));
                        }
                        data_grad(i,2*nb_fun) = data_grad(i,nb_grad_comp-1)/T_scale;
                    }
            	} else {
            		double fun_scale = funs_scale(funorcomp_type_in(i));
                    if (intervalles(i)==0)
                    {
                        data_grad(i,0) = (data_grad(i,0)*P0_fun
                                            + data_grad(i,1)*P1_fun
                                            + data_grad(i,2)*P2_fun)/(funif_scale*fun_scale);
                        data_grad(i,1) = data_grad(i,3)/(P_scale*fun_scale);
                        data_grad(i,2) = data_grad(i,4)/T_scale;
                    }
                    else if (intervalles(i)==1)
                    {
                        data_grad(i,0) = (data_grad(i,0)*P1_fun
                                            + data_grad(i,1)*P2_fun)/(funif_scale*fun_scale);
                        data_grad(i,1) = data_grad(i,2)/(P_scale*fun_scale);
                        data_grad(i,2) = data_grad(i,3)/(P_scale*fun_scale);
                        data_grad(i,3) = data_grad(i,4)/T_scale;
                    }
                    else if (intervalles(i)==2)
                    {
                        data_grad(i,0) = data_grad(i,0)*P2_fun/(funif_scale*fun_scale);
                        data_grad(i,1) /= P_scale*fun_scale;
                        data_grad(i,2) /= P_scale*fun_scale;
                        data_grad(i,3) /= P_scale*fun_scale;
                        data_grad(i,4) /= T_scale;
                    }
                    else if (intervalles(i)<Nint-3)
                    {
                        data_grad(i,0) /= P_scale*fun_scale;
                        data_grad(i,1) /= P_scale*fun_scale;
                        data_grad(i,2) /= P_scale*fun_scale;
                        data_grad(i,3) /= P_scale*fun_scale;
                        data_grad(i,4) /= T_scale;
                    }
                    else if (intervalles(i)==Nint-3)
                    {
                        data_grad(i,0) /= P_scale*fun_scale;
                        data_grad(i,1) /= P_scale*fun_scale;
                        data_grad(i,2) /= P_scale*fun_scale;
                        data_grad(i,3) = data_grad(i,3)*Pendm2_fun/(funif_scale*fun_scale);
                        data_grad(i,4) /= T_scale;
                    }
                    else if (intervalles(i)==Nint-2) 
                    {
                        data_grad(i,0) /= P_scale*fun_scale;
                        data_grad(i,1) /= P_scale*fun_scale;
                        data_grad(i,2) = (data_grad(i,2)*Pendm2_fun
                                            + data_grad(i,3)*Pendm1_fun)/(funif_scale*fun_scale);
                        data_grad(i,3) = data_grad(i,4)/T_scale;
                    }
                    else if (intervalles(i)==Nint-1)
                    {
                        data_grad(i,0) /= P_scale*fun_scale;
                        data_grad(i,1) = (data_grad(i,1)*Pendm2_fun
                                            + data_grad(i,2)*Pendm1_fun
                                            + data_grad(i,3)*Pend_fun)/(funif_scale*fun_scale);
                        data_grad(i,2) = data_grad(i,4)/T_scale;
                    }
            	}
            }
//---------------------------------------------------------------------------
        } else { // gradient of time is not computed
            for (i=0; i<(int) time_points_in.size(); i++)
            { 
                if (intervalles(i)==0)
                {
                    for (kk=0; kk<nb_fun; kk++) {
                        data_grad(i,kk) = (data_grad(i,kk)*P0_fun
                                              + data_grad(i,kk+nb_fun)*P1_fun
                                              + data_grad(i,kk+2*nb_fun)*P2_fun)/(funif_scale*funs_scale(kk%nb_fun));
                        data_grad(i,kk+nb_fun) = data_grad(i,kk+3*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                    }
                }
                else if (intervalles(i)==1)
                {
                    for (kk=0; kk<nb_fun; kk++) {
                        data_grad(i,kk) = (data_grad(i,kk)*P1_fun
                                                + data_grad(i,kk+nb_fun)*P2_fun)/(funif_scale*funs_scale(kk%nb_fun));
                    }
                    for (kk=0; kk<2*nb_fun; kk++) {
                        data_grad(i,kk+nb_fun) = data_grad(i,kk+2*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                    }
                }
                else if (intervalles(i)==2)
                {
                    for (kk=0; kk<nb_fun; kk++) {
                        data_grad(i,kk) = data_grad(i,kk)*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                    }
                    for (kk=0; kk<3*nb_fun; kk++) {
                        data_grad(i,kk+nb_fun) /= P_scale*funs_scale(kk%nb_fun);
                    }
                }
                else if (intervalles(i)<Nint-3)
                {
                    for (kk=0; kk<4*nb_fun; kk++) {
                        data_grad(i,kk) /= P_scale*funs_scale(kk%nb_fun);
                    }
                }
                else if (intervalles(i)==Nint-3)
                {
                    for (kk=0; kk<3*nb_fun; kk++) {
                        data_grad(i,kk) /= P_scale*funs_scale(kk%nb_fun);
                    }
                    for (kk=0; kk<nb_fun; kk++) {
                        data_grad(i,kk+3*nb_fun) = data_grad(i,kk+3*nb_fun)*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                    }
                }
                else if (intervalles(i)==Nint-2) 
                {
                    for (kk=0; kk<2*nb_fun; kk++) {
                        data_grad(i,kk) /= P_scale*funs_scale(kk%nb_fun);
                    }
                    for (kk=0; kk<nb_fun; kk++) {
                        data_grad(i,kk+2*nb_fun) = (data_grad(i,kk+2*nb_fun)*Pendm2_fun
                                                        + data_grad(i,kk+3*nb_fun)*Pendm1_fun)/(funif_scale*funs_scale(kk%nb_fun));
                    }
                }
                else if (intervalles(i)==Nint-1)
                {
                    for (kk=0; kk<nb_fun; kk++) {
                        data_grad(i,kk) /= P_scale*funs_scale(kk%nb_fun);
                    }
                    for (kk=0; kk<nb_fun; kk++) {
                        data_grad(i,kk+nb_fun) = (data_grad(i,kk+1*nb_fun)*Pendm2_fun
                                                        + data_grad(i,kk+2*nb_fun)*Pendm1_fun
                                                        + data_grad(i,kk+3*nb_fun)*Pend_fun)/(funif_scale*funs_scale(kk%nb_fun));
                    }
                }
            }
        } // end gradient of time is not computed
    } // end limit conditions are defined

}
    

void bspline::uncompress_grad(boost::numeric::ublas::vector<double> *all_t, 
	                 boost::numeric::ublas::matrix< boost::numeric::ublas::vector<double> > *data_grad_comp, 
	                 boost::numeric::ublas::matrix< boost::numeric::ublas::vector<double> > *data_grad_uncomp, int nb_t)
{
	using namespace boost::numeric::ublas;
	
	int nb_rows = (int) (*data_grad_comp).size1();

	if (nb_fun!=nb_rows) {
		std::cout << "in bspline::uncompress_grad, decompression of data having " << nb_rows << " rows" << std::endl;
	}
	if (nb_t>(int) (*data_grad_comp).size2()) {
		std::cout << "error in bspline::uncompress_grad, data_grad_comp does not have nb_t columns" << std::endl;
	}	
	if (nb_rows!=(*data_grad_uncomp).size1()) {
		std::cout << "error in bspline::uncompress_grad, data_grad_uncomp does not have as much rows as data_grad_comp" << std::endl;
	}
	if (nb_t>(int) (*data_grad_uncomp).size2()) {
		std::cout << "error in bspline::uncompress_grad, data_grad_uncomp does not have nb_t columns" << std::endl;
	}	
	if (nb_grad_fun==(*data_grad_comp)(0,0).size()) {
		if (nb_fun!=(*data_grad_uncomp).size1()) {
			std::cout << "error in bspline::uncompress_grad, data_grad_uncomp does not have nb_fun rows" << std::endl;
		}	
	} else if (nb_grad_comp!=(*data_grad_comp)(0,0).size()) {
		std::cout << "error in bspline::uncompress_grad, data_grad_comp does not have nb_grad_fun or nb_grad_comp components of gradient" << std::endl;
	}	
	if (nb_grad_uncomp!=(*data_grad_uncomp)(0,0).size()) {
		std::cout << "error in bspline::uncompress_grad, data_grad_uncomp does not have nb_grad_uncomp components of gradient" << std::endl;
		std::cout << "nb_grad_uncomp " << nb_grad_uncomp << std::endl;
		std::cout << "number of components of gradient of data_grad_uncomp " << (*data_grad_uncomp)(0,0).size() << std::endl;
	}	

	int i, j, k, kk;

    //computation of vector N such that N[j] is the number of time sample in [0,j*Tint]
    calc_N(all_t, nb_t);
    
    //---------------------------------------------------------------------------
    //uncompression of joint angles q, dq, ddq
    if (nb_grad_fun==(*data_grad_comp)(0,0).size()) {
    	
    //---------------------------------------------------------------------------
        if (!defined_lc) { // the limit conditions are not defined
    
            if (grad_of_time) {
                for (k=0 ; k<Nint ; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_fun; i++) {
                            for (kk=0; kk<nb_grad_uncomp-1; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            (*data_grad_uncomp)(i,j)(k*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k+1)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k+2)*nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k+3)*nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                        }
                    }
                }
            } else {
                for (k=0 ; k<Nint ; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_fun; i++) {
                            for (kk=0; kk<nb_grad_uncomp; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            (*data_grad_uncomp)(i,j)(k*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k+1)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k+2)*nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k+3)*nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                        }
                    }
                }
            }
            
    //---------------------------------------------------------------------------
        } else if (defined_lc & !free_fun_lc) { // the limit conditions are completely defined

    //---------------------------------------------------------------------------
            if (grad_of_time) { // gradient of time is computed
                for (j=N[0] ; j<N[1] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                    }
                }
                for (j=N[1] ; j<N[2] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                    }
                }
                for (j=N[2] ; j<N[3] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(2*nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                    }
                }
                for (k=3 ; k<Nint-3 ; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_fun; i++) {
                            for (kk=0; kk<nb_grad_uncomp; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            (*data_grad_uncomp)(i,j)((k-3)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k-2)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k-1)*nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)(k*nb_fun+i)     = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                        }
                    }
                }
                for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)((Nint-6)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-5)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-4)*nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                   }
                }
                for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)((Nint-5)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-4)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                    }
                }
                for (j=N[Nint-1] ; j<N[Nint] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)((Nint-4)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                    }
                }
    //---------------------------------------------------------------------------
            } else { // gradient of time is not computed
                for (j=N[0] ; j<N[1] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                    }
                }
                for (j=N[1] ; j<N[2] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                    }
                }
                for (j=N[2] ; j<N[3] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(2*nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                    }
                }
                for (k=3 ; k<Nint-3 ; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_fun; i++) {
                            for (kk=0; kk<nb_grad_uncomp; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            (*data_grad_uncomp)(i,j)((k-3)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k-2)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k-1)*nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)(k*nb_fun+i)     = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                        }
                    }
                }
                for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)((Nint-6)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-5)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-4)*nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                   }
                }
                for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)((Nint-5)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-4)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                    }
                }
                for (j=N[Nint-1] ; j<N[Nint] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)((Nint-4)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                    }
                }
            } // end gradient of time is not computed
        
    //---------------------------------------------------------------------------
        } else { // the limit conditions are defined, limit conditions of fun are free
            calc_P_fun_0();
            calc_P_fun_f();
            
            if (grad_of_time) { // gradient of time is computed
                for (j=N[0] ; j<N[1] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(i) = ((*data_grad_comp)(i,j)(0)*P0_fun
                                                        + (*data_grad_comp)(i,j)(1)*P1_fun
                                                         + (*data_grad_comp)(i,j)(2)*P2_fun)/(funif_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                    }
                }
                for (j=N[1] ; j<N[2] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(i) = ((*data_grad_comp)(i,j)(0)*P1_fun
                                                        + (*data_grad_comp)(i,j)(1)*P2_fun)/(funif_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(2*nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                    }
                }
                for (j=N[2] ; j<N[3] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(i) = (*data_grad_comp)(i,j)(0)*P2_fun/(funif_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(2*nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(3*nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                    }
                }
                for (k=3 ; k<Nint-3 ; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_fun; i++) {
                            for (kk=0; kk<nb_grad_uncomp; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            (*data_grad_uncomp)(i,j)((k-2)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k-1)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)(k*nb_fun+i)     = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k+1)*nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                        }
                    }
                }
                for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)((Nint-5)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-4)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-3)*nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-2)*nb_fun+i) = (*data_grad_comp)(i,j)(3)*Pendm2_fun/(funif_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                   }
                }
                for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)((Nint-4)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-3)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-2)*nb_fun+i) = ((*data_grad_comp)(i,j)(2)*Pendm2_fun
                                                                       + (*data_grad_comp)(i,j)(3)*Pendm1_fun)/(funif_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                    }
                }
                for (j=N[Nint-1] ; j<N[Nint] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)((Nint-3)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-2)*nb_fun+i) = ((*data_grad_comp)(i,j)(1)*Pendm2_fun
                                                                       + (*data_grad_comp)(i,j)(2)*Pendm1_fun
                                                                        + (*data_grad_comp)(i,j)(3)*Pend_fun)/(funif_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(4)/T_scale;
                    }
                }
    //---------------------------------------------------------------------------
            } else { // gradient of time is not computed
                 for (j=N[0] ; j<N[1] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(i) = ((*data_grad_comp)(i,j)(0)*P0_fun
                                                        + (*data_grad_comp)(i,j)(1)*P1_fun
                                                         + (*data_grad_comp)(i,j)(2)*P2_fun)/(funif_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                    }
                }
                for (j=N[1] ; j<N[2] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(i) = ((*data_grad_comp)(i,j)(0)*P1_fun
                                                        + (*data_grad_comp)(i,j)(1)*P2_fun)/(funif_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(2*nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                    }
                }
                for (j=N[2] ; j<N[3] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(i) = (*data_grad_comp)(i,j)(0)*P2_fun/(funif_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(2*nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)(3*nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                    }
                }
                for (k=3 ; k<Nint-3 ; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_fun; i++) {
                            for (kk=0; kk<nb_grad_uncomp; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            (*data_grad_uncomp)(i,j)((k-2)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k-1)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)(k*nb_fun+i)     = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                            (*data_grad_uncomp)(i,j)((k+1)*nb_fun+i) = (*data_grad_comp)(i,j)(3)/(P_scale*funs_scale(i));
                        }
                    }
                }
                for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)((Nint-5)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-4)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-3)*nb_fun+i) = (*data_grad_comp)(i,j)(2)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-2)*nb_fun+i) = (*data_grad_comp)(i,j)(3)*Pendm2_fun/(funif_scale*funs_scale(i));
                   }
                }
                for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)((Nint-4)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-3)*nb_fun+i) = (*data_grad_comp)(i,j)(1)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-2)*nb_fun+i) = ((*data_grad_comp)(i,j)(2)*Pendm2_fun
                                                                       + (*data_grad_comp)(i,j)(3)*Pendm1_fun)/(funif_scale*funs_scale(i));
                    }
                }
                for (j=N[Nint-1] ; j<N[Nint] ; j++)
                { 
                    for (i=0; i<nb_fun; i++) {
                        for (kk=0; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)((Nint-3)*nb_fun+i) = (*data_grad_comp)(i,j)(0)/(P_scale*funs_scale(i));
                        (*data_grad_uncomp)(i,j)((Nint-2)*nb_fun+i) = ((*data_grad_comp)(i,j)(1)*Pendm2_fun
                                                                       + (*data_grad_comp)(i,j)(2)*Pendm1_fun
                                                                        + (*data_grad_comp)(i,j)(3)*Pend_fun)/(funif_scale*funs_scale(i));
                    }
                }
            } // end gradient of time is not computed
        } // end limit conditions are defined
        
    //---------------------------------------------------------------------------
    } else { //uncompression of data depending on all angles
        
        if (!defined_lc) { // the limit conditions are not defined
    
            if (grad_of_time) { // gradient of time is computed
                for (k=0 ; k<Nint ; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_rows; i++) {
                            for (kk=0; kk<k*nb_fun; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            for (kk=0; kk<nb_grad_comp-1; kk++) {
                                (*data_grad_uncomp)(i,j)(kk+k*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                            }
                            for (kk=k*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(nb_grad_comp-1)/T_scale;
                        }
                    }
                }
            } else {
                for (k=0 ; k<Nint ; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_rows; i++) {
                            for (kk=0; kk<k*nb_fun; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            for (kk=0; kk<nb_grad_comp; kk++) {
                                (*data_grad_uncomp)(i,j)(kk+k*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                            }
                            for (kk=k*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                        }
                    }
                }
            }
            
    //---------------------------------------------------------------------------
        } else if (defined_lc & !free_fun_lc) { // the limit conditions are completely defined
        
            if (grad_of_time) { // gradient of time is computed
                for (k=0; k<3; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_rows; i++) {
                            for (kk=0; kk<(k+1)*nb_fun; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = (*data_grad_comp)(i,j)(kk+(3-k)*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                            }
                            for (kk=(k+1)*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(nb_grad_comp-1)/T_scale;
                        }
                    }
                }
                for (k=3 ; k<Nint-3 ; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_rows; i++) {
                            for (kk=0; kk<(k-3)*nb_fun; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            for (kk=0; kk<nb_grad_comp-1; kk++) {
                                (*data_grad_uncomp)(i,j)(kk+(k-3)*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                            }
                            for (kk=(k-3)*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(nb_grad_comp-1)/T_scale;
                        }
                    }
                }
                for (k=Nint-3; k<Nint; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_rows; i++) {
                            for (kk=0; kk<(k-3)*nb_fun; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            for (kk=0; kk<(Nint-k)*nb_fun; kk++) {
                                (*data_grad_uncomp)(i,j)(kk+(k-3)*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                            }
                            (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(nb_grad_comp-1)/T_scale;
                        }
                    }
                }
    //---------------------------------------------------------------------------
            } else { // gradient of time is not computed
                for (k=0; k<3; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_rows; i++) {
                            for (kk=0; kk<(k+1)*nb_fun; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = (*data_grad_comp)(i,j)(kk+(3-k)*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                            }
                            for (kk=(k+1)*nb_fun; kk<nb_grad_uncomp; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                        }
                    }
                }
                for (k=3 ; k<Nint-3 ; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_rows; i++) {
                            for (kk=0; kk<(k-3)*nb_fun; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            for (kk=0; kk<nb_grad_comp; kk++) {
                                (*data_grad_uncomp)(i,j)(kk+(k-3)*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                            }
                            for (kk=(k-3)*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                        }
                    }
                }
                for (k=Nint-3; k<Nint; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_rows; i++) {
                            for (kk=0; kk<(k-3)*nb_fun; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            for (kk=0; kk<(Nint-k)*nb_fun; kk++) {
                                (*data_grad_uncomp)(i,j)(kk+(k-3)*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                            }
                        }
                    }
                }
            } // end gradient of time is not computed

    //---------------------------------------------------------------------------
        } else { // the limit conditions are defined, limit conditions of fun are free
            calc_P_fun_0();
            calc_P_fun_f();
            
            if (grad_of_time) { // gradient of time is computed
                for (j=N[0] ; j<N[1] ; j++)
                { 
                    for (i=0; i<nb_rows; i++) {
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = ((*data_grad_comp)(i,j)(kk)*P0_fun
                                                            + (*data_grad_comp)(i,j)(kk+nb_fun)*P1_fun
                                                             + (*data_grad_comp)(i,j)(kk+2*nb_fun)*P2_fun)/(funif_scale*funs_scale(kk%nb_fun));
                            (*data_grad_uncomp)(i,j)(kk+nb_fun) = (*data_grad_comp)(i,j)(kk+3*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=2*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(nb_grad_comp-1)/T_scale;
                    }
                }
                for (j=N[1] ; j<N[2] ; j++)
                { 
                    for (i=0; i<nb_rows; i++) {
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = ((*data_grad_comp)(i,j)(kk)*P1_fun
                                                             + (*data_grad_comp)(i,j)(kk+nb_fun)*P2_fun)/(funif_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=0; kk<2*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+nb_fun) = (*data_grad_comp)(i,j)(kk+2*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=3*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(nb_grad_comp-1)/T_scale;
                    }
                }
                for (j=N[2] ; j<N[3] ; j++)
                { 
                    for (i=0; i<nb_rows; i++) {
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = (*data_grad_comp)(i,j)(kk)*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=0; kk<3*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+nb_fun) = (*data_grad_comp)(i,j)(kk+nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=4*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(nb_grad_comp-1)/T_scale;
                    }
                }
                for (k=3 ; k<Nint-3 ; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_rows; i++) {
                            for (kk=0; kk<(k-2)*nb_fun; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            for (kk=0; kk<nb_grad_comp-1; kk++) {
                                (*data_grad_uncomp)(i,j)(kk+(k-2)*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                            }
                            for (kk=(k-2)*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(nb_grad_comp-1)/T_scale;
                        }
                    }
                }
                for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
                { 
                    for (i=0; i<nb_rows; i++) {
                        for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        for (kk=0; kk<3*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+(Nint-5)*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+(Nint-2)*nb_fun) = (*data_grad_comp)(i,j)(kk+3*nb_fun)*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        }
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(nb_grad_comp-1)/T_scale;
                    }
                }
                for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
                { 
                    for (i=0; i<nb_rows; i++) {
                        for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        for (kk=0; kk<2*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+(Nint-4)*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+(Nint-2)*nb_fun) = ((*data_grad_comp)(i,j)(kk+2*nb_fun)*Pendm2_fun
                                                                            + (*data_grad_comp)(i,j)(kk+3*nb_fun)*Pendm1_fun)/(funif_scale*funs_scale(kk%nb_fun));
                        }
                        (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(nb_grad_comp-1)/T_scale;
                    }
                }
                for (j=N[Nint-1] ; j<N[Nint] ; j++)
                { 
                    for (i=0; i<nb_rows; i++) {
                        for (kk=0; kk<(Nint-3)*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+(Nint-3)*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+(Nint-2)*nb_fun) = ((*data_grad_comp)(i,j)(kk+1*nb_fun)*Pendm2_fun
                                                                            + (*data_grad_comp)(i,j)(kk+2*nb_fun)*Pendm1_fun
                                                                             + (*data_grad_comp)(i,j)(kk+3*nb_fun)*Pend_fun)/(funif_scale*funs_scale(kk%nb_fun));
                        }
                         (*data_grad_uncomp)(i,j)(nb_grad_uncomp-1) = (*data_grad_comp)(i,j)(nb_grad_comp-1)/T_scale;
                    }
                }
    //---------------------------------------------------------------------------
            } else { // gradient of time is not computed
                for (j=N[0] ; j<N[1] ; j++)
                { 
                    for (i=0; i<nb_rows; i++) {
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = ((*data_grad_comp)(i,j)(kk)*P0_fun
                                                            + (*data_grad_comp)(i,j)(kk+nb_fun)*P1_fun
                                                             + (*data_grad_comp)(i,j)(kk+2*nb_fun)*P2_fun)/(funif_scale*funs_scale(kk%nb_fun));
                            (*data_grad_uncomp)(i,j)(kk+nb_fun) = (*data_grad_comp)(i,j)(kk+3*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=2*nb_fun; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                    }
                }
                for (j=N[1] ; j<N[2] ; j++)
                { 
                    for (i=0; i<nb_rows; i++) {
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = ((*data_grad_comp)(i,j)(kk)*P1_fun
                                                             + (*data_grad_comp)(i,j)(kk+nb_fun)*P2_fun)/(funif_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=0; kk<2*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+nb_fun) = (*data_grad_comp)(i,j)(kk+2*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=3*nb_fun; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                    }
                }
                for (j=N[2] ; j<N[3] ; j++)
                { 
                    for (i=0; i<nb_rows; i++) {
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = (*data_grad_comp)(i,j)(kk)*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=0; kk<3*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+nb_fun) = (*data_grad_comp)(i,j)(kk+nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=4*nb_fun; kk<nb_grad_uncomp; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                    }
                }
                for (k=3 ; k<Nint-3 ; k++)
                {
                    for (j=N[k] ; j<N[k+1] ; j++)
                    { 
                        for (i=0; i<nb_rows; i++) {
                            for (kk=0; kk<(k-2)*nb_fun; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                            for (kk=0; kk<nb_grad_comp; kk++) {
                                (*data_grad_uncomp)(i,j)(kk+(k-2)*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                            }
                            for (kk=(k-2)*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) {
                                (*data_grad_uncomp)(i,j)(kk) = 0;
                            }
                        }
                    }
                }
                for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
                { 
                    for (i=0; i<nb_rows; i++) {
                        for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        for (kk=0; kk<3*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+(Nint-5)*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+(Nint-2)*nb_fun) = (*data_grad_comp)(i,j)(kk+3*nb_fun)*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        }
                    }
                }
                for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
                { 
                    for (i=0; i<nb_rows; i++) {
                        for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        for (kk=0; kk<2*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+(Nint-4)*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+(Nint-2)*nb_fun) = ((*data_grad_comp)(i,j)(kk+2*nb_fun)*Pendm2_fun
                                                                            + (*data_grad_comp)(i,j)(kk+3*nb_fun)*Pendm1_fun)/(funif_scale*funs_scale(kk%nb_fun));
                        }
                    }
                }
                for (j=N[Nint-1] ; j<N[Nint] ; j++)
                { 
                    for (i=0; i<nb_rows; i++) {
                        for (kk=0; kk<(Nint-3)*nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk) = 0;
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+(Nint-3)*nb_fun) = (*data_grad_comp)(i,j)(kk)/(P_scale*funs_scale(kk%nb_fun));
                        }
                        for (kk=0; kk<nb_fun; kk++) {
                            (*data_grad_uncomp)(i,j)(kk+(Nint-2)*nb_fun) = ((*data_grad_comp)(i,j)(kk+1*nb_fun)*Pendm2_fun
                                                                            + (*data_grad_comp)(i,j)(kk+2*nb_fun)*Pendm1_fun
                                                                             + (*data_grad_comp)(i,j)(kk+3*nb_fun)*Pend_fun)/(funif_scale*funs_scale(kk%nb_fun));
                        }
                    }
                }            
            } // end gradient of time is not computed
        } // end limit conditions are defined
    } // end uncompression of which type of data
    
}


void bspline::uncompress_grad_int(boost::numeric::ublas::vector<double> *all_t, 
	                 boost::numeric::ublas::vector<double> *data, 
	                 boost::numeric::ublas::matrix<double> *data_grad_comp, 
	                 boost::numeric::ublas::vector<double> *data_grad_uncomp, int nb_t, double Tfinal)
{
	using namespace boost::numeric::ublas;
	
	if (nb_t>(int) (*data).size()) {
		std::cout << "error in bspline::uncompress_grad_int, data does not have nb_t elements" << std::endl;
	}	
	if (nb_t>(int) (*data_grad_comp).size1()) {
		std::cout << "error in bspline::uncompress_grad_int, data_grad_comp does not have nb_t lines" << std::endl;
	}	
	if (nb_grad_comp!=(*data_grad_comp).size2()) {
		std::cout << "error in bspline::uncompress_grad_int, data_grad_comp does not have nb_grad_comp columns" << std::endl;
	}	
	if (nb_grad_uncomp!=(*data_grad_uncomp).size()) {
		std::cout << "error in bspline::uncompress_grad_int, data_grad_uncomp does not have nb_grad_uncomp components of gradient" << std::endl;
		std::cout << "nb_grad_uncomp " << nb_grad_uncomp << std::endl;
		std::cout << "number of components of gradient of data_grad_uncomp " << (*data_grad_uncomp).size() << std::endl;
	}	
    
	int j, k, kk;

    //computation of vector N such that N[j] is the number of time sample in [0,j*Tint]
    calc_N(all_t, nb_t);
    
    if (N[0]==N[1]) {
    	std::cout << "error in bspline::uncompress_grad_int, the first intervalle must have at least one point" << std::endl;
		std::cout << N[0] << " " << N[1] << std::endl;
    }
    if (N[Nint-1]==N[Nint]) {
    	std::cout << "error in bspline::uncompress_grad_int, the last intervalle must have at least one point" << std::endl;
    }
    
    //---------------------------------------------------------------------------
    if (!defined_lc) { // the limit conditions are not defined

        if (grad_of_time) { // gradient of time is computed
            for (kk=0; kk<nb_grad_uncomp; kk++) {
                (*data_grad_uncomp)(kk) = 0;
            }
            for (kk=0; kk<nb_grad_comp-1; kk++) {
                (*data_grad_uncomp)(kk) += (*data_grad_comp)(0,kk) * ((*all_t)(1)-(*all_t)(0))/2.;
            }
            (*data_grad_uncomp)(nb_grad_uncomp-1) = (*data_grad_comp)(0,nb_grad_comp-1) * ((*all_t)(1)-(*all_t)(0))/2.
                                                     + (*data)(0) * ((*all_t)(1)-(*all_t)(0))/2./Tfinal;
            for (j=N[0]+1 ; j<N[1] ; j++)
            { 
                for (kk=0; kk<nb_grad_comp-1; kk++) {
                    (*data_grad_uncomp)(kk) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                         + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
            }
            for (k=1 ; k<Nint-1 ; k++)
            {
                for (j=N[k] ; j<N[k+1] ; j++)
                { 
                    for (kk=0; kk<nb_grad_comp-1; kk++) {
                        (*data_grad_uncomp)(kk+k*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                    }
                    (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                             + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
                }
            }
            for (j=N[Nint-1] ; j<N[Nint]-1 ; j++) // for the last intervalle we do not take the last point
            { 
                for (kk=0; kk<nb_grad_comp-1; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-1)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                         + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
            }
            for (kk=0; kk<nb_grad_comp-1; kk++) {
                (*data_grad_uncomp)(kk+(Nint-1)*nb_fun) += (*data_grad_comp)(N[Nint]-1,kk) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2.;
            }
            (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(N[Nint]-1,nb_grad_comp-1) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2.
                                                      + (*data)(N[Nint]-1) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2./Tfinal;
            for (kk=0; kk<nb_grad_uncomp-1; kk++) {
            	(*data_grad_uncomp)(kk) /= P_scale*funs_scale(kk%nb_fun);
            }
            (*data_grad_uncomp)(nb_grad_uncomp-1) /= T_scale;
//---------------------------------------------------------------------------
        } else { // gradient of time is not computed
            for (kk=0; kk<nb_grad_uncomp; kk++) {
                (*data_grad_uncomp)(kk) = 0;
            }
            for (kk=0; kk<nb_grad_comp; kk++) {
                (*data_grad_uncomp)(kk) += (*data_grad_comp)(0,kk) * ((*all_t)(1)-(*all_t)(0))/2.;
            }
            for (j=N[0]+1 ; j<N[1] ; j++)
            { 
                for (kk=0; kk<nb_grad_comp; kk++) {
                    (*data_grad_uncomp)(kk) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
            }
            for (k=1 ; k<Nint-1 ; k++)
            {
                for (j=N[k] ; j<N[k+1] ; j++)
                { 
                    for (kk=0; kk<nb_grad_comp; kk++) {
                        (*data_grad_uncomp)(kk+k*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                    }
                }
            }
            for (j=N[Nint-1] ; j<N[Nint]-1 ; j++) // for the last intervalle we do not take the last point
            { 
                for (kk=0; kk<nb_grad_comp; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-1)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
            }
            for (kk=0; kk<nb_grad_comp; kk++) {
                (*data_grad_uncomp)(kk+(Nint-1)*nb_fun) += (*data_grad_comp)(N[Nint]-1,kk) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2.;
            }
            for (kk=0; kk<nb_grad_uncomp-1; kk++) {
            	(*data_grad_uncomp)(kk) /= P_scale*funs_scale(kk%nb_fun);
            }
        }
        
//---------------------------------------------------------------------------
    } else if (defined_lc & !free_fun_lc) { // the limit conditions are completely defined
    
        if (grad_of_time) { // gradient of time is computed
            for (kk=0; kk<nb_grad_uncomp; kk++) {
                (*data_grad_uncomp)(kk) = 0;
            }
            for (kk=0; kk<nb_fun; kk++) { //for the first point
                (*data_grad_uncomp)(kk) += (*data_grad_comp)(0,kk+3*nb_fun) * ((*all_t)(1)-(*all_t)(0))/2.;
            }
            (*data_grad_uncomp)(nb_grad_uncomp-1) = (*data_grad_comp)(0,nb_grad_comp-1) * ((*all_t)(1)-(*all_t)(0))/2.
                                                     + (*data)(0) * ((*all_t)(1)-(*all_t)(0))/2./Tfinal;
            for (j=N[0]+1 ; j<N[1] ; j++) //for the other points of the first intervalle
            { 
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk) += (*data_grad_comp)(j,kk+3*nb_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                         + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
            }
            for (k=1; k<3; k++) //for the second and third intervalles
            {
                for (j=N[k] ; j<N[k+1] ; j++)
                { 
                    for (kk=0; kk<(k+1)*nb_fun; kk++) {
                        (*data_grad_uncomp)(kk) += (*data_grad_comp)(j,kk+(3-k)*nb_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                    }
                    (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                             + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
                }
            }
            for (k=3 ; k<Nint-3 ; k++) // for intermediate intervalles
            {
                for (j=N[k] ; j<N[k+1] ; j++)
                { 
                    for (kk=0; kk<nb_grad_comp-1; kk++) {
                        (*data_grad_uncomp)(kk+(k-3)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                    }
                    (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                             + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
                }
            }
            for (k=Nint-3; k<Nint-1; k++) // for Nint-3 and Nint-2 intervalles
            {
                for (j=N[k] ; j<N[k+1] ; j++)
                { 
                    for (kk=0; kk<(Nint-k)*nb_fun; kk++) {
                        (*data_grad_uncomp)(kk+(k-3)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                    }
                    (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                             + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
                }
            }
            for (j=N[Nint-1] ; j<N[Nint]-1 ; j++) // for Nint-1 intervalle without last point
            { 
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-4)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                         + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
            }
            for (kk=0; kk<nb_fun; kk++) { // for last point
                (*data_grad_uncomp)(kk+(Nint-4)*nb_fun) += (*data_grad_comp)(N[Nint]-1,kk) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2.;
            }
            (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(N[Nint]-1,nb_grad_comp-1) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2.
                                                     + (*data)(N[Nint]-1) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2./Tfinal;
            for (kk=0; kk<nb_grad_uncomp-1; kk++) {
            	(*data_grad_uncomp)(kk) /= P_scale*funs_scale(kk%nb_fun);
            }
            (*data_grad_uncomp)(nb_grad_uncomp-1) /= T_scale;
//---------------------------------------------------------------------------
        } else { // gradient of time is not computed
            for (kk=0; kk<nb_grad_uncomp; kk++) {
                (*data_grad_uncomp)(kk) = 0;
            }
            for (kk=0; kk<nb_fun; kk++) { //for the first point
                (*data_grad_uncomp)(kk) += (*data_grad_comp)(0,kk+3*nb_fun) * ((*all_t)(1)-(*all_t)(0))/2.;
            }
            for (j=N[0]+1 ; j<N[1] ; j++) //for the other points of the first intervalle
            { 
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk) += (*data_grad_comp)(j,kk+3*nb_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
            }
            for (k=1; k<3; k++) //for the second and third intervalles
            {
                for (j=N[k] ; j<N[k+1] ; j++)
                { 
                    for (kk=0; kk<(k+1)*nb_fun; kk++) {
                        (*data_grad_uncomp)(kk) += (*data_grad_comp)(j,kk+(3-k)*nb_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                    }
                }
            }
            for (k=3 ; k<Nint-3 ; k++) // for intermediate intervalles
            {
                for (j=N[k] ; j<N[k+1] ; j++)
                { 
                    for (kk=0; kk<nb_grad_comp; kk++) {
                        (*data_grad_uncomp)(kk+(k-3)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                    }
                }
            }
            for (k=Nint-3; k<Nint-1; k++) // for Nint-3 and Nint-2 intervalles
            {
                for (j=N[k] ; j<N[k+1] ; j++)
                { 
                    for (kk=0; kk<(Nint-k)*nb_fun; kk++) {
                        (*data_grad_uncomp)(kk+(k-3)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                    }
                }
            }
            for (j=N[Nint-1] ; j<N[Nint]-1 ; j++) // for Nint-1 intervalle without last point
            { 
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-4)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
            }
            for (kk=0; kk<nb_fun; kk++) { // for last point
                (*data_grad_uncomp)(kk+(Nint-4)*nb_fun) += (*data_grad_comp)(N[Nint]-1,kk) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2.;
            }
            for (kk=0; kk<nb_grad_uncomp-1; kk++) {
            	(*data_grad_uncomp)(kk) /= P_scale*funs_scale(kk%nb_fun);
            }
        } // end gradient of time is not computed

//---------------------------------------------------------------------------
    } else { // the limit conditions are defined, limit conditions of fun are free
        calc_P_fun_0();
        calc_P_fun_f();
        
        if (grad_of_time) { // gradient of time is computed
            for (kk=0; kk<nb_grad_uncomp; kk++) {
                (*data_grad_uncomp)(kk) = 0;
            }
            for (kk=0; kk<nb_fun; kk++) { // for first point
                (*data_grad_uncomp)(kk) += ((*data_grad_comp)(0,kk)*P0_fun
                                                + (*data_grad_comp)(0,kk+nb_fun)*P1_fun
                                                 + (*data_grad_comp)(0,kk+2*nb_fun)*P2_fun) * ((*all_t)(1)-(*all_t)(0))/2.;
                (*data_grad_uncomp)(kk+nb_fun) += (*data_grad_comp)(0,kk+3*nb_fun) * ((*all_t)(1)-(*all_t)(0))/2.;
            }
            (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(0,nb_grad_comp-1) * ((*all_t)(1)-(*all_t)(0))/2.
                                                     + (*data)(0) * ((*all_t)(1)-(*all_t)(0))/2./Tfinal;
            for (j=N[0]+1 ; j<N[1] ; j++) // for first intervalle without first point
            { 
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk) += ((*data_grad_comp)(j,kk)*P0_fun
                                                    + (*data_grad_comp)(j,kk+nb_fun)*P1_fun
                                                     + (*data_grad_comp)(j,kk+2*nb_fun)*P2_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                    (*data_grad_uncomp)(kk+nb_fun) += (*data_grad_comp)(j,kk+3*nb_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                         + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
            }
            for (j=N[1] ; j<N[2] ; j++)  //for second intervalle
            { 
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk) += ((*data_grad_comp)(j,kk)*P1_fun
                                                     + (*data_grad_comp)(j,kk+nb_fun)*P2_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                for (kk=0; kk<2*nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+nb_fun) += (*data_grad_comp)(j,kk+2*nb_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                         + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
            }
            for (j=N[2] ; j<N[3] ; j++)
            { 
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk) += (*data_grad_comp)(j,kk)*P2_fun * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                for (kk=0; kk<3*nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+nb_fun) += (*data_grad_comp)(j,kk+nb_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                         + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
            }
            for (k=3 ; k<Nint-3 ; k++)
            {
                for (j=N[k] ; j<N[k+1] ; j++)
                { 
                    for (kk=0; kk<nb_grad_comp-1; kk++) {
                        (*data_grad_uncomp)(kk+(k-2)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                    }
                    (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                             + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
                }
            }
            for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
            { 
                for (kk=0; kk<3*nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-5)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-2)*nb_fun) += (*data_grad_comp)(j,kk+3*nb_fun)*Pendm2_fun * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                         + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
            }
            for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
            { 
                for (kk=0; kk<2*nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-4)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-2)*nb_fun) += ((*data_grad_comp)(j,kk+2*nb_fun)*Pendm2_fun
                                                                    + (*data_grad_comp)(j,kk+3*nb_fun)*Pendm1_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                         + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
            }
            for (j=N[Nint-1] ; j<N[Nint]-1 ; j++)
            { 
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-3)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-2)*nb_fun) += ((*data_grad_comp)(j,kk+1*nb_fun)*Pendm2_fun
                                                                    + (*data_grad_comp)(j,kk+2*nb_fun)*Pendm1_fun
                                                                     + (*data_grad_comp)(j,kk+3*nb_fun)*Pend_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(j,nb_grad_comp-1) * ((*all_t)(j+1)-(*all_t)(j-1))/2.
                                                         + (*data)(j) * ((*all_t)(j+1)-(*all_t)(j-1))/2./Tfinal;
            }
            for (kk=0; kk<nb_fun; kk++) {
                (*data_grad_uncomp)(kk+(Nint-3)*nb_fun) += (*data_grad_comp)(N[Nint]-1,kk) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2.;
            }
            for (kk=0; kk<nb_fun; kk++) {
                (*data_grad_uncomp)(kk+(Nint-2)*nb_fun) += ((*data_grad_comp)(N[Nint]-1,kk+1*nb_fun)*Pendm2_fun
                                                                + (*data_grad_comp)(N[Nint]-1,kk+2*nb_fun)*Pendm1_fun
                                                                 + (*data_grad_comp)(N[Nint]-1,kk+3*nb_fun)*Pend_fun) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2.;
            }
            (*data_grad_uncomp)(nb_grad_uncomp-1) += (*data_grad_comp)(N[Nint]-1,nb_grad_comp-1) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2.
                                                     + (*data)(N[Nint]-1) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2./Tfinal;
            for (kk=0; kk<nb_fun; kk++) {
            	(*data_grad_uncomp)(kk) /= funif_scale*funs_scale(kk%nb_fun);
            }
            for (kk=nb_fun; kk<(Nint-2)*nb_fun; kk++) {
            	(*data_grad_uncomp)(kk) /= P_scale*funs_scale(kk%nb_fun);
            }
            for (kk=(Nint-2)*nb_fun; kk<(Nint-1)*nb_fun; kk++) {
            	(*data_grad_uncomp)(kk) /= funif_scale*funs_scale(kk%nb_fun);
            }
            (*data_grad_uncomp)(nb_grad_uncomp-1) /= T_scale;
//---------------------------------------------------------------------------
        } else { // gradient of time is not computed
            for (kk=0; kk<nb_grad_uncomp; kk++) {
                (*data_grad_uncomp)(kk) = 0;
            }
            for (kk=0; kk<nb_fun; kk++) { // for first point
                (*data_grad_uncomp)(kk) += ((*data_grad_comp)(0,kk)*P0_fun
                                                + (*data_grad_comp)(0,kk+nb_fun)*P1_fun
                                                 + (*data_grad_comp)(0,kk+2*nb_fun)*P2_fun) * ((*all_t)(1)-(*all_t)(0))/2.;
                (*data_grad_uncomp)(kk+nb_fun) += (*data_grad_comp)(0,kk+3*nb_fun) * ((*all_t)(1)-(*all_t)(0))/2.;
            }
            for (j=N[0]+1 ; j<N[1] ; j++) // for first intervalle without first point
            { 
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk) += ((*data_grad_comp)(j,kk)*P0_fun
                                                    + (*data_grad_comp)(j,kk+nb_fun)*P1_fun
                                                     + (*data_grad_comp)(j,kk+2*nb_fun)*P2_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                    (*data_grad_uncomp)(kk+nb_fun) += (*data_grad_comp)(j,kk+3*nb_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
            }
            for (j=N[1] ; j<N[2] ; j++)  //for second intervalle
            { 
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk) += ((*data_grad_comp)(j,kk)*P1_fun
                                                     + (*data_grad_comp)(j,kk+nb_fun)*P2_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                for (kk=0; kk<2*nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+nb_fun) += (*data_grad_comp)(j,kk+2*nb_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
            }
            for (j=N[2] ; j<N[3] ; j++)
            { 
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk) += (*data_grad_comp)(j,kk)*P2_fun * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                for (kk=0; kk<3*nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+nb_fun) += (*data_grad_comp)(j,kk+nb_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
            }
            for (k=3 ; k<Nint-3 ; k++)
            {
                for (j=N[k] ; j<N[k+1] ; j++)
                { 
                    for (kk=0; kk<nb_grad_comp; kk++) {
                        (*data_grad_uncomp)(kk+(k-2)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                    }
                }
            }
            for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
            { 
                for (kk=0; kk<3*nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-5)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-2)*nb_fun) += (*data_grad_comp)(j,kk+3*nb_fun)*Pendm2_fun * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
            }
            for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
            { 
                for (kk=0; kk<2*nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-4)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-2)*nb_fun) += ((*data_grad_comp)(j,kk+2*nb_fun)*Pendm2_fun
                                                                    + (*data_grad_comp)(j,kk+3*nb_fun)*Pendm1_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
            }
            for (j=N[Nint-1] ; j<N[Nint]-1 ; j++)
            { 
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-3)*nb_fun) += (*data_grad_comp)(j,kk) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
                for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(kk+(Nint-2)*nb_fun) += ((*data_grad_comp)(j,kk+1*nb_fun)*Pendm2_fun
                                                                    + (*data_grad_comp)(j,kk+2*nb_fun)*Pendm1_fun
                                                                     + (*data_grad_comp)(j,kk+3*nb_fun)*Pend_fun) * ((*all_t)(j+1)-(*all_t)(j-1))/2.;
                }
            }
            for (kk=0; kk<nb_fun; kk++) {
                (*data_grad_uncomp)(kk+(Nint-3)*nb_fun) += (*data_grad_comp)(N[Nint]-1,kk) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2.;
            }
            for (kk=0; kk<nb_fun; kk++) {
                (*data_grad_uncomp)(kk+(Nint-2)*nb_fun) += ((*data_grad_comp)(N[Nint]-1,kk+1*nb_fun)*Pendm2_fun
                                                                + (*data_grad_comp)(N[Nint]-1,kk+2*nb_fun)*Pendm1_fun
                                                                 + (*data_grad_comp)(N[Nint]-1,kk+3*nb_fun)*Pend_fun) * ((*all_t)(N[Nint]-1)-(*all_t)(N[Nint]-2))/2.;
            }
            for (kk=0; kk<nb_fun; kk++) {
            	(*data_grad_uncomp)(kk) /= funif_scale*funs_scale(kk%nb_fun);
            }
            for (kk=nb_fun; kk<(Nint-2)*nb_fun; kk++) {
            	(*data_grad_uncomp)(kk) /= P_scale*funs_scale(kk%nb_fun);
            }
            for (kk=(Nint-2)*nb_fun; kk<(Nint-1)*nb_fun; kk++) {
            	(*data_grad_uncomp)(kk) /= funif_scale*funs_scale(kk%nb_fun);
            }
        } // end gradient of time is not computed
    } // end limit conditions are defined
    
}


void bspline::solve_P_0(const boost::numeric::ublas::vector<double> &fun0)
{
    calc_P_fun_0();
    
	for (int i=0; i<nb_fun; i++) {
		(*matP) (i,0) = P0_fun * fun0(i);
		(*matP) (i,1) = P1_fun * fun0(i);
		(*matP) (i,2) = P2_fun * fun0(i);
	}
	
}

void bspline::solve_P_f(const boost::numeric::ublas::vector<double> &funf)
{
    calc_P_fun_f();

	for (int i=0; i<nb_fun; i++) {
		(*matP) (i,nb_P-3) = Pendm2_fun * funf(i);
		(*matP) (i,nb_P-2) = Pendm1_fun * funf(i);
		(*matP) (i,nb_P-1) = Pend_fun * funf(i);
	}
	
}
	

void bspline::calc_P_fun_0()
{
    //double Tint2 = pow(Tint,2);

    f11v = f11(0., 0., 0.);
    df11v = df11(0., 0.);//Tint; remove division by Tint since it is not known 
    df12v = df12(0., 0.);//Tint;  it can be since they simplify later
    d2f11v = d2f11(0.);//Tint2;
    d2f12v = d2f12(0.);//Tint2;
    d2f13v = d2f13(0.);//Tint2;
    P0_fun = 1. / f11v;
    P1_fun = - df11v / f11v / df12v;
    P2_fun = - (df12v * d2f11v - df11v * d2f12v) / df12v / f11v / d2f13v;
    
}

void bspline::calc_P_fun_f()
{
    //double Tint2 = pow(Tint,2);

    f11v = f11(0., 0., 0.);
    df11v = - df11(0., 0.);//Tint; remove division by Tint since it is not known 
    df12v = - df12(0., 0.);//Tint;  it can be since they simplify later
    d2f11v = d2f11(0.);//Tint2;
    d2f12v = d2f12(0.);//Tint2;
    d2f13v = d2f13(0.);//Tint2;
    Pend_fun = 1. / f11v;
    Pendm1_fun = - df11v / f11v / df12v;
    Pendm2_fun = - (df12v * d2f11v - df11v * d2f12v) / df12v / f11v / d2f13v;
	
}
	

