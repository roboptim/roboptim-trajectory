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

#include <iostream>
#include <fstream>

#include <roboptim-trajectory/basis-function.hh>

#ifdef _MSC_VER
# include <minmax.h>
# define fmin min
# define fmax max
#endif

#ifdef __GNUC__
# define min std::min
# define max std::max
#endif

namespace roboptim
{
  BasisFunction::BasisFunction(int nb_P_per_fun_in) : Nint(0), nb_call_def_sizes(0), nb_P_per_fun(nb_P_per_fun_in), all_tn(NULL), matP(NULL)
  {
  }

  BasisFunction::BasisFunction(int nb_fun_in, int nb_P_in, int nb_t_m_in, bool grad_of_time_in, int nb_P_per_fun_in) : Nint(0), nb_call_def_sizes(0), nb_P_per_fun(nb_P_per_fun_in), all_tn(NULL), matP(NULL)
  {
    define_sizes(nb_fun_in, nb_P_in, nb_t_m_in, grad_of_time_in);
  }

  BasisFunction::BasisFunction(int nb_fun_in, int nb_P_in, int nb_t_m_in, bool grad_of_time_in, bool defined_lc_in, bool free_fun_lc_in, int nb_P_per_fun_in) : Nint(0), nb_call_def_sizes(0), nb_P_per_fun(nb_P_per_fun_in), all_tn(NULL), matP(NULL)
  {
    define_sizes(nb_fun_in, nb_P_in, nb_t_m_in, grad_of_time_in, defined_lc_in, free_fun_lc_in);
  }

  BasisFunction::~BasisFunction()
  {
    //dtor
    if (all_tn != NULL)
      delete all_tn;
    if (matP != NULL)
      delete matP;
  }

  void BasisFunction::define_sizes(int nb_fun_in, int nb_P_in, int nb_t_m_in, bool grad_of_time_in)
  {
    grad_of_time = grad_of_time_in;
    defined_lc = false;
    free_fun_lc = false;
    nb_fun = nb_fun_in;
    nb_P = nb_P_in;
    nb_t_m = nb_t_m_in;

    define_sizes();

    _T_fixed = 1.;
  }


  void BasisFunction::define_sizes(int nb_fun_in, int nb_P_in, int nb_t_m_in, bool grad_of_time_in, bool defined_lc_in, bool free_fun_lc_in)
  {
    grad_of_time = grad_of_time_in;
    defined_lc = defined_lc_in;
    free_fun_lc = free_fun_lc_in;
    nb_fun = nb_fun_in;
    nb_P = nb_P_in;
    nb_t_m = nb_t_m_in;

    define_sizes();
  }


  void BasisFunction::define_sizes()
  {
    using namespace boost::numeric;

    nb_call_def_sizes++;

    given_lc = false;
    given_T = false;
    nb_grad_fun = nb_P_per_fun;
    nb_grad_comp = nb_fun*nb_grad_fun;
    nb_grad_uncomp = nb_fun*nb_P;
    if (grad_of_time) { //if computation of time gradient, adding another component
      nb_grad_fun++;
      nb_grad_comp++;
      nb_grad_uncomp++;
    }
    if (defined_lc & !free_fun_lc) { //if initial and final states are fixed
      // there is 6 components of gradient less per function
      nb_grad_uncomp -= 6*nb_fun;
    }
    if (defined_lc & free_fun_lc) { //if initial and final first and second derivatives are fixed
      // there is 4 components of gradient less per function
      nb_grad_uncomp -= 4*nb_fun;
    }
    def_Nint();
    N.resize(Nint+1,false);
    all_tn = new ublas::vector<double> (nb_t_m);
    matP = new ublas::matrix<double> (nb_fun, nb_P);
    fun0_fixed.resize(nb_fun,false);
    funf_fixed.resize(nb_fun,false);

    //for the scaling
    funif_scale = 1.;
    P_scale = 1.;
    T_scale = 1.;
    funs_scale.resize(nb_fun,false);
    for (int i=0; i<nb_fun; i++) {
      funs_scale(i) = 1.;
    }
  }


  void BasisFunction::copy(const BasisFunction &bf)
  {
    define_sizes(bf.nb_fun, bf.nb_P, bf.nb_t_m, bf.grad_of_time, bf.defined_lc, bf.free_fun_lc);

    *matP = *(bf.matP);
    fun0_fixed = bf.fun0_fixed;
    funf_fixed = bf.funf_fixed;
    _T_fixed = bf._T_fixed;
    T = bf.T;
    Tint = bf.Tint;

    funif_scale = bf.funif_scale;
    P_scale = bf.P_scale;
    T_scale = bf.T_scale;
    funs_scale = bf.funs_scale;

    given_lc = bf.given_lc;
    given_T = bf.given_T;

  }


  int BasisFunction::get_nb_call_def_sizes() //give the number of call to define_sizes, so that other class getting sizes from this class can track if sizes have changed
  {
    return nb_call_def_sizes;
  }


  void BasisFunction::def_parameters(boost::numeric::ublas::matrix<double> *P_in, double T_in)
  {
    if (defined_lc && free_fun_lc) {
      std::cout << "error in BasisFunction, def_parameters without definition of funi and funf has been called while it should" << std::endl;
    }

    check_T(T_in);

    T = T_in;

    Tint = T/Nint;

    if (nb_fun!=(*P_in).size1()) {
      std::cout << "error in BasisFunction, P_in does not have nb_fun rows" << std::endl;
    }

    if (!defined_lc) {
      if (nb_P!=(*P_in).size2()) {
        std::cout << "error in BasisFunction, P_in does not have nb_P columns" << std::endl;
      }
      for (int i=0 ; i<nb_fun ; i++)
        {
          for (int j=0 ; j<nb_P ; j++ ) {
            (*matP) (i,j) = (*P_in) (i,j);
          }
        }
    }

    if (defined_lc & !free_fun_lc) {
      if (nb_P-6!=(*P_in).size2()) {
        std::cout << "error in BasisFunction, P_in does not have nb_P-6 columns" << std::endl;
      }
      for (int i=0 ; i<nb_fun ; i++) //columns 0,1,2,nb_P-3,nb_P-2,nb_P-1 of matP fixed and than already defined
        {
          for (int j=3 ; j<nb_P-3 ; j++ ) {
            (*matP) (i,j) = (*P_in) (i,j-3);
          }
        }
    }
  }


  void BasisFunction::def_parameters(boost::numeric::ublas::matrix<double> *P_in,
                                     const boost::numeric::ublas::vector<double> &fun0,
                                     const boost::numeric::ublas::vector<double> &funf, double T_in)
  {
    if (!(defined_lc && free_fun_lc)) {
      std::cout << "error in BasisFunction, def_parameters with definition of fun0 and funf has been called while it should not" << std::endl;
    }

    check_T(T_in);

    T = T_in;

    Tint = T/Nint;

    if (nb_fun!=(*P_in).size1()) {
      std::cout << "error in BasisFunction, P_in does not have nb_fun rows" << std::endl;
    }
    if (nb_P-6!=(*P_in).size2()) {
      std::cout << "error in BasisFunction, P_in does not have nb_P-6 columns" << std::endl;
    }
    solve_P_0(fun0);
    solve_P_f(funf);
    for (int i=0 ; i<nb_fun ; i++)
      {
    	for (int j=3 ; j<nb_P-3 ; j++ ) {
          (*matP) (i,j) = (*P_in) (i,j-3);
    	}
      }

  }


  void BasisFunction::define_fixed_lc(
                                      const boost::numeric::ublas::vector<double> &fun0,
                                      const boost::numeric::ublas::vector<double> &funf)
  {
    if (!(defined_lc && !free_fun_lc)) {
      std::cout << "error in BasisFunction, define_fixed_lc has been called while it should not" << std::endl;
    }

    given_lc = true;

    fun0_fixed = fun0;
    funf_fixed = funf;

    solve_P_0(fun0);
    solve_P_f(funf);

  }


  void BasisFunction::get_fixed_lc(boost::numeric::ublas::vector<double> &fun0,
                                   boost::numeric::ublas::vector<double> &funf) const
  {
    fun0 = fun0_fixed;
    funf = funf_fixed;
  }

  void BasisFunction::define_fixed_time(double T_in)
  {
    given_T = true;

    _T_fixed = T_in;
  }

  void BasisFunction::get_fixed_time(double &T_out)
  {
    T_out = _T_fixed;
  }

  double BasisFunction::get_free_time()
  {
    return T;
  }

  bool BasisFunction::get_defined_lc()
  {
    return defined_lc;
  }

  bool BasisFunction::get_free_fun_lc()
  {
    return free_fun_lc;
  }

  int BasisFunction::get_nb_fun()
  {
    return nb_fun;
  }

  int BasisFunction::get_nb_P_free() //number of free parameters for each joint (if limit conditions are defined it is less than nb_P)
  {
    if (defined_lc) {
      return nb_P-6;
    } else {
      return nb_P;
    }
  }


  int BasisFunction::get_Nint()
  {
    if (Nint == 0) {
      def_Nint();
    }
    return Nint;
  }


  void BasisFunction::save_parameters(const char* file_name, boost::numeric::ublas::matrix<double> *P_in, double T_in)
  {
    if (defined_lc) {
      std::cout << "error in BasisFunction, bad function save_parameters has been called" << std::endl;
    }
    if (nb_fun!=(*P_in).size1()) {
      std::cout << "error in BasisFunction, P_in does not have nb_fun rows" << std::endl;
    }
    if (nb_P!=(*P_in).size2() && !defined_lc) {
      std::cout << "error in BasisFunction, P_in does not have nb_P columns" << std::endl;
    }
    if (nb_P-6!=(*P_in).size2() && defined_lc) {
      std::cout << "error in BasisFunction, P_in does not have nb_P-6 columns" << std::endl;
    }

    check_T(T_in);

  }


  void BasisFunction::save_parameters(
                                      const char* file_name,
                                      boost::numeric::ublas::matrix<double> *P_in,
                                      const boost::numeric::ublas::vector<double> &fun0_in,
                                      const boost::numeric::ublas::vector<double> &funf_in, double T_in)
  {
    if (!defined_lc) {
      std::cout << "error in BasisFunction, bad function save_parameters has been called" << std::endl;
    }
    if (nb_fun!=(*P_in).size1()) {
      std::cout << "error in BasisFunction, P_in does not have nb_fun rows" << std::endl;
    }
    if (nb_P-6!=(*P_in).size2()) {
      std::cout << "error in BasisFunction, P_in does not have nb_P-6 columns" << std::endl;
    }

    check_T(T_in);

    std::ofstream file(file_name);

    file.precision(16); // be carefull to write the data with full precision
    file.flags(std::ios_base::scientific);

    file << nb_fun << std::endl;

    file << nb_P << std::endl;

    file << nb_t_m << std::endl;

    file << grad_of_time << std::endl;

    file << defined_lc << std::endl;

    file << free_fun_lc << std::endl;

    for (int j=0; j<nb_fun; j++) {
      file << fun0_in(j) << " ";
    }
    file << std::endl;

    for (int i=0; i<nb_P-6; i++) {
      for (int j=0; j<nb_fun; j++) {
        file << (*P_in)(j,i) << " ";
      }
      file << std::endl;
    }

    for (int j=0; j<nb_fun; j++) {
      file << funf_in(j) << " ";
    }
    file << std::endl;

    //if (grad_of_time) {
    file << T_in << std::endl;
    //}

    file.close();

  }


  void BasisFunction::load_parameters(const char* file_name, boost::numeric::ublas::matrix<double> *P_out, double &T_out)
  {
    std::cout << "error in BasisFunction::load_parameters, this overloaded function is not implemented" << std::endl;

    if (defined_lc) {
      std::cout << "error in BasisFunction, bad function load_parameters has been called" << std::endl;
    }
    if (nb_fun!=(*P_out).size1()) {
      std::cout << "error in BasisFunction, P_out does not have nb_fun rows" << std::endl;
    }
    if (nb_P!=(*P_out).size2() && !defined_lc) {
      std::cout << "error in BasisFunction, P_out does not have nb_P columns" << std::endl;
    }
    if (nb_P-6!=(*P_out).size2() && defined_lc) {
      std::cout << "error in BasisFunction, P_out does not have nb_P-6 columns" << std::endl;
    }

    //check_T(T_out);

  }


  void BasisFunction::load_parameters(
                                      const char* file_name,
                                      boost::numeric::ublas::matrix<double> *P_out,
                                      boost::numeric::ublas::vector<double> &fun0_out,
                                      boost::numeric::ublas::vector<double> &funf_out, double &T_out)
  {
    if (!defined_lc) {
      std::cout << "error in BasisFunction, bad function load_parameters has been called" << std::endl;
    }
    if (nb_fun!=(*P_out).size1()) {
      std::cout << "error in BasisFunction, P_out does not have nb_fun rows" << std::endl;
    }
    if (nb_P-6!=(*P_out).size2()) {
      std::cout << "error in BasisFunction, P_out does not have nb_P-6 columns" << std::endl;
    }

    std::ifstream file(file_name);

    if (!file.is_open()) {
      std::cout << "error in BasisFunction::load_parameters, could not open the file " << file_name << std::endl;
    }

    int nb_fun_in, nb_P_in, nb_t_m_in;
    bool grad_of_time_in, defined_lc_in, free_fun_lc_in;

    file >> nb_fun_in;
    if (nb_fun_in!=nb_fun) {
      std::cout << "error in BasisFunction::load_parameters, the data read does not have the same nb_fun parameter" << std::endl;
    }

    file >> nb_P_in;
    if (nb_P_in!=nb_P) {
      std::cout << "error in BasisFunction::load_parameters, the data read does not have the same nb_P parameter" << std::endl;
    }

    file >> nb_t_m_in;
    if (nb_t_m_in!=nb_t_m) {
      std::cout << "error in BasisFunction::load_parameters, the data read does not have the same nb_t_m parameter" << std::endl;
    }

    file >> grad_of_time_in;
    if (grad_of_time_in!=grad_of_time) {
      std::cout << "error in BasisFunction::load_parameters, the data read does not have the same grad_of_time parameter" << std::endl;
    }

    file >> defined_lc_in;
    if (defined_lc_in!=defined_lc) {
      std::cout << "error in BasisFunction::load_parameters, the data read does not have the same defined_lc parameter" << std::endl;
    }

    file >> free_fun_lc_in;
    if (free_fun_lc_in!=free_fun_lc) {
      std::cout << "error in BasisFunction::load_parameters, the data read does not have the same free_fun_lc parameter" << std::endl;
    }

    for (int j=0; j<nb_fun; j++) {
      file >> fun0_out(j);
    }

    for (int i=0; i<nb_P-6; i++) {
      for (int j=0; j<nb_fun; j++) {
        file >> (*P_out)(j,i);
      }
    }

    for (int j=0; j<nb_fun; j++) {
      file >> funf_out(j);
    }

    //if (grad_of_time) {
    file >> T_out;
    //}

    check_T(T_out);

    file.close();

  }


  void BasisFunction::load_parameters_and_sizes(
                                                const char* file_name,
                                                boost::numeric::ublas::matrix<double> *P_out,
                                                boost::numeric::ublas::vector<double> &fun0_out,
                                                boost::numeric::ublas::vector<double> &funf_out, double &T_out)
  {

    std::ifstream file(file_name);

    if (!file.is_open()) {
      std::cout << "error in BasisFunction::load_parameters, could not open the file " << file_name << std::endl;
    }

    file >> nb_fun;

    file >> nb_P;

    file >> nb_t_m;

    file >> grad_of_time;

    file >> defined_lc;
    if (!defined_lc) {
      std::cout << "error in BasisFunction, bad function load_parameters has been called" << std::endl;
      file.close();
      return;
    }

    file >> free_fun_lc;

    (*P_out).resize(nb_fun,nb_P-6, false);
    fun0_out.resize(nb_fun, false);
    funf_out.resize(nb_fun, false);

    for (int j=0; j<nb_fun; j++) {
      file >> fun0_out(j);
    }

    for (int i=0; i<nb_P-6; i++) {
      for (int j=0; j<nb_fun; j++) {
        file >> (*P_out)(j,i);
      }
    }

    for (int j=0; j<nb_fun; j++) {
      file >> funf_out(j);
    }

    //if (grad_of_time) {
    file >> T_out;
    //}

    check_T(T_out);

    file.close();

    define_sizes();

    // We define the fixed parameters if there are some
    if (defined_lc && !free_fun_lc) {
      define_fixed_lc(fun0_out, funf_out);
    }

    if (!grad_of_time) {
      define_fixed_time(T_out);
    }

    // We define the free parameters
    if (defined_lc && free_fun_lc) {
      def_parameters(P_out, fun0_out, funf_out, T_out);
    } else {
      def_parameters(P_out, T_out);
    }

  }

  void BasisFunction::load_parameters_and_sizes(const char* file_name)
  {
    boost::numeric::ublas::matrix<double> P_local;
    boost::numeric::ublas::vector<double> fun0_local;
    boost::numeric::ublas::vector<double> funf_local;
    double T_local;

    load_parameters_and_sizes(file_name, &P_local, fun0_local, funf_local, T_local);
  }


  void BasisFunction::define_parameters_scale(double P_scale_in, double T_scale_in, double funif_scale_in,
                                              const boost::numeric::ublas::vector<double> & funs_scale_in)
  {
    if (funs_scale_in.size() != nb_fun) {
      std::cout << "error in BasisFunction::define_parameters_scale, funs_scale_in does not have nb_fun elements" << std::endl;
    }

    P_scale = P_scale_in;
    T_scale = T_scale_in;
    funif_scale = funif_scale_in;
    for (int i=0; i<nb_fun; i++) {
      funs_scale(i) = funs_scale_in(i);
    }
  }


  void BasisFunction::get_parameters_scale(double &P_scale_out, double &T_scale_out, double &funif_scale_out,
                                           boost::numeric::ublas::vector<double> &funs_scale_out)
  {
    if (funs_scale_out.size() != nb_fun) {
      std::cout << "error in BasisFunction::get_parameters_scale, funs_scale_in does not have nb_fun elements" << std::endl;
    }

    P_scale_out = P_scale;
    T_scale_out = T_scale;
    funif_scale_out = funif_scale;
    for (int i=0; i<nb_fun; i++) {
      funs_scale_out(i) = funs_scale(i);
    }
  }


  void BasisFunction::get_parameters_scale(boost::numeric::ublas::vector<double> &x_scale_out)
  {

    std::cout << "warning in BasisFunction, function get_parameters_scale never debugged" << std::endl;


    if (nb_grad_uncomp!=x_scale_out.size()) {
      std::cout << "error in BasisFunction, x_scale_out is not of size nb_grad_uncomp" << std::endl;
    }

    int rank = 0;

    if (!defined_lc) {
      for (int i=0; i<nb_P; i++) {
        for (int j=0; j<nb_fun; j++) {
          x_scale_out(rank) = P_scale*funs_scale(j);
          rank++;
        }
      }
    } else {
      if (free_fun_lc) {
        for (int j=0; j<nb_fun; j++) {
          x_scale_out(rank) = funif_scale*funs_scale(j);
          rank++;
        }
      }
      for (int i=0; i<nb_P-6; i++) {
        for (int j=0; j<nb_fun; j++) {
          x_scale_out(rank) = P_scale*funs_scale(j);
          rank++;
        }
      }
      if (free_fun_lc) {
        for (int j=0; j<nb_fun; j++) {
          x_scale_out(rank) = funif_scale*funs_scale(j);
          rank++;
        }
      }
    }

    if (grad_of_time) {
      x_scale_out(rank) = T_scale;
    }

  }


  void BasisFunction::get_scale_for_auto_scale(int &nb_x_scale_out,
                                               boost::numeric::ublas::vector<int> &nb_p_x_scale_p_out,
                                               boost::numeric::ublas::matrix<int> &p_x_scale_p_out)
  {
    //scaling parameters ordering : funif_scale, P_scale (if free_fun_lc) and T_scale

    if (grad_of_time) {
      nb_x_scale_out = 1;
    } else {
      nb_x_scale_out = 0;
    }

    if (!defined_lc) {
      nb_x_scale_out += nb_fun;
      nb_p_x_scale_p_out.resize(nb_grad_uncomp,false);
      p_x_scale_p_out.resize(nb_grad_uncomp,1,false);
      for (int i=0; i<nb_P; i++) {
        for (int j=0; j<nb_fun; j++) {
          nb_p_x_scale_p_out(i*nb_fun+j) = 1;
          p_x_scale_p_out(i*nb_fun+j,0) = j;
        }
      }
    } else if (!free_fun_lc) {
      nb_x_scale_out += nb_fun;
      nb_p_x_scale_p_out.resize(nb_grad_uncomp,false);
      p_x_scale_p_out.resize(nb_grad_uncomp,1,false);
      for (int i=0; i<nb_P-6; i++) {
        for (int j=0; j<nb_fun; j++) {
          nb_p_x_scale_p_out(i*nb_fun+j) = 1;
          p_x_scale_p_out(i*nb_fun+j,0) = j;
        }
      }
    } else {
      nb_x_scale_out += nb_fun + 1;
      nb_p_x_scale_p_out.resize(nb_grad_uncomp,false);
      p_x_scale_p_out.resize(nb_grad_uncomp,2,false);
      for (int j=0; j<nb_fun; j++) {
        nb_p_x_scale_p_out(j) = 1;
        p_x_scale_p_out(j,0) = j;
      }
      for (int i=1; i<nb_P-5; i++) {
        for (int j=0; j<nb_fun; j++) {
          nb_p_x_scale_p_out(i*nb_fun+j) = 2;
          p_x_scale_p_out(i*nb_fun+j,0) = j;
          p_x_scale_p_out(i*nb_fun+j,1) = nb_fun;
        }
      }
      for (int j=0; j<nb_fun; j++) {
        nb_p_x_scale_p_out((nb_P-5)*nb_fun+j) = 1;
        p_x_scale_p_out((nb_P-5)*nb_fun+j,0) = j;
      }
    }

    if (grad_of_time) {
      nb_p_x_scale_p_out(nb_grad_uncomp-1) = 1;
      p_x_scale_p_out(nb_grad_uncomp-1,0) = nb_x_scale_out-1;
    }

  }


  void BasisFunction::set_scale_from_auto_scale(boost::numeric::ublas::vector<double> &p_x_scale)
  {
    //scaling parameters ordering : funif_scale, P_scale (if free_fun_lc) and T_scale

    for (int i=0; i<nb_fun; i++) {
      funs_scale(i) *= p_x_scale(i);
    }

    if (free_fun_lc) {
      P_scale *= p_x_scale(nb_fun);
      if (grad_of_time) {
        T_scale *= p_x_scale(nb_fun+1);
      }
    } else {
      if (grad_of_time) {
        T_scale *= p_x_scale(nb_fun);
      }
    }

  }


  void BasisFunction::convert_parameters_P2x(const boost::numeric::ublas::matrix<double> *P_in, const double T_in, double* x_out)
  {
    if ((defined_lc && free_fun_lc)) {
      std::cout << "error in BasisFunction, bad function convert_parameters_P2x has been called" << std::endl;
    }
    if (nb_fun!=(*P_in).size1()) {
      std::cout << "error in BasisFunction, P_in does not have nb_fun rows" << std::endl;
    }
    if (nb_P!=(*P_in).size2() && !defined_lc) {
      std::cout << "error in BasisFunction, P_in does not have nb_P columns" << std::endl;
    }
    if (nb_P-6!=(*P_in).size2() && defined_lc) {
      std::cout << "error in BasisFunction, P_in does not have nb_P-6 columns" << std::endl;
    }

    //scale_parameters(P_in,T_in);

    if (!defined_lc) {
      for (int i=0; i<nb_P; i++) {
        for (int j=0; j<nb_fun; j++) {
          *x_out = (*P_in)(j,i)*(P_scale*funs_scale(j));
          ++x_out;
        }
      }
    } else {
      for (int i=0; i<nb_P-6; i++) {
        for (int j=0; j<nb_fun; j++) {
          *x_out = (*P_in)(j,i)*(P_scale*funs_scale(j));
          ++x_out;
        }
      }
    }

    if (grad_of_time) {
      *x_out = T_in*T_scale;
    }

  }


  void BasisFunction::convert_parameters_P2x(const boost::numeric::ublas::matrix<double> *P_in,
                                             const boost::numeric::ublas::vector<double> &fun0_in,
                                             const boost::numeric::ublas::vector<double> &funf_in,
                                             const double T_in, double* x_out)
  {
    if (!defined_lc || (defined_lc && !free_fun_lc)) {
      std::cout << "error in BasisFunction, bad function convert_parameters_P2x has been called" << std::endl;
    }
    if (nb_fun!=(*P_in).size1()) {
      std::cout << "error in BasisFunction, P_in does not have nb_fun rows" << std::endl;
    }
    if (nb_P-6!=(*P_in).size2()) {
      std::cout << "error in BasisFunction, P_in does not have nb_P-6 columns" << std::endl;
    }

    //scale_parameters(P_in,fun0_in,funf_in,T_in);

    for (int j=0; j<nb_fun; j++) {
      *x_out = fun0_in(j)*(funif_scale*funs_scale(j));
      ++x_out;
    }

    for (int i=0; i<nb_P-6; i++) {
      for (int j=0; j<nb_fun; j++) {
        *x_out = (*P_in)(j,i)*(P_scale*funs_scale(j));
        ++x_out;
      }
    }

    for (int j=0; j<nb_fun; j++) {
      *x_out = funf_in(j)*(funif_scale*funs_scale(j));
      ++x_out;
    }

    if (grad_of_time) {
      *x_out = T_in*T_scale;
    }

  }


  void BasisFunction::convert_parameters_P2x(
                                             const boost::numeric::ublas::matrix<bool> *P_in,
                                             const bool T_in,
                                             bool* x_out)
  {
    if ((defined_lc && free_fun_lc)) {
      std::cout << "error in BasisFunction, bad function convert_parameters_P2x has been called" << std::endl;
    }
    if (nb_fun!=(*P_in).size1()) {
      std::cout << "error in BasisFunction, P_in does not have nb_fun rows" << std::endl;
    }
    if (nb_P!=(*P_in).size2() && !defined_lc) {
      std::cout << "error in BasisFunction, P_in does not have nb_P columns" << std::endl;
    }
    if (nb_P-6!=(*P_in).size2() && defined_lc) {
      std::cout << "error in BasisFunction, P_in does not have nb_P-6 columns" << std::endl;
    }

    //scale_parameters(P_in,T_in);

    if (!defined_lc) {
      for (int i=0; i<nb_P; i++) {
        for (int j=0; j<nb_fun; j++) {
          *x_out = (*P_in)(j,i);
          ++x_out;
        }
      }
    } else {
      for (int i=0; i<nb_P-6; i++) {
        for (int j=0; j<nb_fun; j++) {
          *x_out = (*P_in)(j,i);
          ++x_out;
        }
      }
    }

    if (grad_of_time) {
      *x_out = T_in;
    }

  }


  void BasisFunction::convert_parameters_P2x(const boost::numeric::ublas::matrix<bool> *P_in,
                                             const boost::numeric::ublas::vector<bool> &fun0_in,
                                             const boost::numeric::ublas::vector<bool> &funf_in, const bool T_in, bool* x_out)
  {
    if (!defined_lc || (defined_lc && !free_fun_lc)) {
      std::cout << "error in BasisFunction, bad function convert_parameters_P2x has been called" << std::endl;
    }
    if (nb_fun!=(*P_in).size1()) {
      std::cout << "error in BasisFunction, P_in does not have nb_fun rows" << std::endl;
    }
    if (nb_P-6!=(*P_in).size2()) {
      std::cout << "error in BasisFunction, P_in does not have nb_P-6 columns" << std::endl;
    }

    //scale_parameters(P_in,fun0_in,funf_in,T_in);

    for (int j=0; j<nb_fun; j++) {
      *x_out = fun0_in(j);
      ++x_out;
    }

    for (int i=0; i<nb_P-6; i++) {
      for (int j=0; j<nb_fun; j++) {
        *x_out = (*P_in)(j,i);
        ++x_out;
      }
    }

    for (int j=0; j<nb_fun; j++) {
      *x_out = funf_in(j);
      ++x_out;
    }

    if (grad_of_time) {
      *x_out = T_in;
    }

  }


  void BasisFunction::convert_parameters_x2P(const double* x_in, boost::numeric::ublas::matrix<double> *P_out, double &T_out)
  {
    if ((defined_lc && free_fun_lc)) {
      std::cout << "error in BasisFunction, bad function convert_parameters_x2P has been called" << std::endl;
    }
    if (nb_fun!=(*P_out).size1()) {
      std::cout << "error in BasisFunction, P_out does not have nb_fun rows" << std::endl;
    }
    if (nb_P!=(*P_out).size2() && !defined_lc) {
      std::cout << "error in BasisFunction, P_out does not have nb_P columns" << std::endl;
    }
    if (nb_P-6!=(*P_out).size2() && defined_lc) {
      std::cout << "error in BasisFunction, P_out does not have nb_P-6 columns" << std::endl;
    }

    if (!defined_lc) {
      for (int i=0; i<nb_P; i++) {
        for (int j=0; j<nb_fun; j++) {
          (*P_out)(j,i) = *x_in/(P_scale*funs_scale(j));
          ++x_in;
        }
      }
    } else {
      for (int i=0; i<nb_P-6; i++) {
        for (int j=0; j<nb_fun; j++) {
          (*P_out)(j,i) = *x_in/(P_scale*funs_scale(j));
          ++x_in;
        }
      }
    }

    if (grad_of_time) {
      T_out = *x_in/T_scale;
    } else {
      T_out = _T_fixed;
    }

    //unscale_parameters(P_out,T_out);

  }


  void BasisFunction::convert_parameters_x2P(const double* x_in,
                                             boost::numeric::ublas::matrix<double> *P_out,
                                             boost::numeric::ublas::vector<double> &fun0_out,
                                             boost::numeric::ublas::vector<double> &funf_out, double &T_out)
  {
    if (!defined_lc || (defined_lc && !free_fun_lc)) {
      std::cout << "error in BasisFunction, bad function convert_parameters_x2P has been called" << std::endl;
    }
    if (nb_fun!=(*P_out).size1()) {
      std::cout << "error in BasisFunction, P_out does not have nb_fun rows" << std::endl;
    }
    if (nb_P-6!=(*P_out).size2()) {
      std::cout << "error in BasisFunction, P_out does not have nb_P-6 columns" << std::endl;
    }

    for (int j=0; j<nb_fun; j++) {
      fun0_out(j) = *x_in/(funif_scale*funs_scale(j));
      ++x_in;
    }

    for (int i=0; i<nb_P-6; i++) {
      for (int j=0; j<nb_fun; j++) {
        (*P_out)(j,i) = *x_in/(P_scale*funs_scale(j));
        ++x_in;
      }
    }

    for (int j=0; j<nb_fun; j++) {
      funf_out(j) = *x_in/(funif_scale*funs_scale(j));
      ++x_in;
    }

    if (grad_of_time) {
      T_out = *x_in/T_scale;
    } else {
      T_out = _T_fixed;
    }

    //unscale_parameters(P_out,fun0_out,funf_out,T_out);

  }


  void BasisFunction::change_P_parameter_in_x(boost::numeric::ublas::vector<double> &x_in_out, int i, int j, double P_ij_in)
  {
    if (i > nb_fun-1) {
      std::cout << "error in BasisFunction, i is more than nb_fun" << std::endl;
    }
    if (defined_lc && (j > nb_P-6 - 1)) {
      std::cout << "error in BasisFunction, j is more than nb_P-6" << std::endl;
    }
    if (!defined_lc && (j > nb_P - 1)) {
      std::cout << "error in BasisFunction, j is more than nb_P" << std::endl;
    }

    if (defined_lc && free_fun_lc) {
      x_in_out[nb_fun+i+j*nb_fun] = P_ij_in*(P_scale*funs_scale(i));
    } else {
      x_in_out[i+j*nb_fun] = P_ij_in*(P_scale*funs_scale(i));
    }

  }

  void BasisFunction::change_fun0_parameter_in_x(boost::numeric::ublas::vector<double> &x_in_out, int i, double fun0_i_in)
  {
    if (!defined_lc || (defined_lc && !free_fun_lc)) {
      std::cout << "error in BasisFunction, change_qi_parameter_in_x has been called while qi not defined or not free" << std::endl;
    }

    x_in_out[i] = fun0_i_in*(funif_scale*funs_scale(i));

  }

  void BasisFunction::change_funf_parameter_in_x(boost::numeric::ublas::vector<double> &x_in_out, int i, double funf_i_in)
  {
    if (!defined_lc || (defined_lc && !free_fun_lc)) {
      std::cout << "error in BasisFunction, change_qf_parameter_in_x has been called while qi not defined or not free" << std::endl;
    }

    x_in_out[(nb_P-6+1)*nb_fun+i] = funf_i_in*(funif_scale*funs_scale(i));
  }

  void BasisFunction::change_T_parameter_in_x(boost::numeric::ublas::vector<double> &x_in_out, double T_in)
  {
    if (grad_of_time) {
      x_in_out[nb_grad_uncomp-1] = T_in*T_scale;
    } else {
      std::cout << "warning in BasisFunction, change_T_parameter_in_x has been called while time is fixed" << std::endl;
    }

  }

  void BasisFunction::change_P_parameter_in_x(boost::numeric::ublas::vector<bool> &x_in_out, int i, int j, bool P_ij_in)
  {
    if (i > nb_fun-1) {
      std::cout << "error in BasisFunction, i is more than nb_fun" << std::endl;
    }
    if (defined_lc && (j > nb_P-6 - 1)) {
      std::cout << "error in BasisFunction, j is more than nb_P-6" << std::endl;
    }
    if (!defined_lc && (j > nb_P - 1)) {
      std::cout << "error in BasisFunction, j is more than nb_P" << std::endl;
    }

    if (defined_lc && free_fun_lc) {
      x_in_out[nb_fun+i+j*nb_fun] = P_ij_in;
    } else {
      x_in_out[i+j*nb_fun] = P_ij_in;
    }

  }

  void BasisFunction::change_fun0_parameter_in_x(boost::numeric::ublas::vector<bool> &x_in_out, int i, bool fun0_i_in)
  {
    if (!defined_lc || (defined_lc && !free_fun_lc)) {
      std::cout << "error in BasisFunction, change_qi_parameter_in_x has been called while qi not defined or not free" << std::endl;
    }

    x_in_out[i] = fun0_i_in;

  }

  void BasisFunction::change_funf_parameter_in_x(boost::numeric::ublas::vector<bool> &x_in_out, int i, bool funf_i_in)
  {
    if (!defined_lc || (defined_lc && !free_fun_lc)) {
      std::cout << "error in BasisFunction, change_qf_parameter_in_x has been called while qi not defined or not free" << std::endl;
    }

    x_in_out[(nb_P-6+1)*nb_fun+i] = funf_i_in;
  }

  void BasisFunction::change_T_parameter_in_x(boost::numeric::ublas::vector<bool> &x_in_out, bool T_in)
  {
    if (grad_of_time) {
      x_in_out[nb_grad_uncomp-1] = T_in;
    } else {
      std::cout << "warning in BasisFunction, change_T_parameter_in_x has been called while time is fixed" << std::endl;
    }

  }


  double BasisFunction::get_P_parameter_in_x(const boost::numeric::ublas::vector<double> &x_in, int i, int j) const
  {
    if (i > nb_fun-1) {
      std::cout << "error in BasisFunction, i is more than nb_fun" << std::endl;
    }
    if (defined_lc && (j > nb_P-6 - 1)) {
      std::cout << "error in BasisFunction, j is more than nb_P-6" << std::endl;
    }
    if (!defined_lc && (j > nb_P - 1)) {
      std::cout << "error in BasisFunction, j is more than nb_P" << std::endl;
    }

    if (defined_lc && free_fun_lc) {
      return x_in[nb_fun+i+j*nb_fun]/(P_scale*funs_scale(i));
    } else {
      return x_in[i+j*nb_fun]/(P_scale*funs_scale(i));
    }

  }

  double BasisFunction::get_fun0_parameter_in_x(const boost::numeric::ublas::vector<double> &x_in, int i) const
  {
    if (!defined_lc || (defined_lc && !free_fun_lc)) {
      std::cout << "error in BasisFunction, get_qi_parameter_in_x has been called while qi not defined or not free" << std::endl;
    }

    return x_in[i]/(funif_scale*funs_scale(i));

  }

  double BasisFunction::get_funf_parameter_in_x(const boost::numeric::ublas::vector<double> &x_in, int i) const
  {
    if (!defined_lc || (defined_lc && !free_fun_lc)) {
      std::cout << "error in BasisFunction, get_qf_parameter_in_x has been called while qi not defined or not free" << std::endl;
    }

    return x_in[(nb_P-6+1)*nb_fun+i]/(funif_scale*funs_scale(i));

  }

  double BasisFunction::get_T_parameter_in_x(const boost::numeric::ublas::vector<double> &x_in) const
  {
    if (grad_of_time) {
      return x_in[nb_grad_uncomp-1]/T_scale;
    } else {
      std::cout << "warning in BasisFunction, get_T_parameter_in_x has been called while time is fixed" << std::endl;
      return 0.;
    }

  }


  bool BasisFunction::get_P_parameter_in_x(const boost::numeric::ublas::vector<bool> &x_in, int i, int j) const
  {
    if (i > nb_fun-1) {
      std::cout << "error in BasisFunction, i is more than nb_fun" << std::endl;
    }
    if (defined_lc && (j > nb_P-6 - 1)) {
      std::cout << "error in BasisFunction, j is more than nb_P-6" << std::endl;
    }
    if (!defined_lc && (j > nb_P - 1)) {
      std::cout << "error in BasisFunction, j is more than nb_P" << std::endl;
    }

    if (defined_lc && free_fun_lc) {
      return x_in[nb_fun+i+j*nb_fun];
    } else {
      return x_in[i+j*nb_fun];
    }

  }

  bool BasisFunction::get_fun0_parameter_in_x(const boost::numeric::ublas::vector<bool> &x_in, int i) const
  {
    if (!defined_lc || (defined_lc && !free_fun_lc)) {
      std::cout << "error in BasisFunction, get_qi_parameter_in_x has been called while qi not defined or not free" << std::endl;
    }

    return x_in[i];

  }

  bool BasisFunction::get_funf_parameter_in_x(const boost::numeric::ublas::vector<bool> &x_in, int i) const
  {
    if (!defined_lc || (defined_lc && !free_fun_lc)) {
      std::cout << "error in BasisFunction, get_qf_parameter_in_x has been called while qi not defined or not free" << std::endl;
    }

    return x_in[(nb_P-6+1)*nb_fun+i];

  }

  bool BasisFunction::get_T_parameter_in_x(const boost::numeric::ublas::vector<bool> &x_in) const
  {
    if (grad_of_time) {
      return x_in[nb_grad_uncomp-1];
    } else {
      std::cout << "warning in BasisFunction, get_T_parameter_in_x has been called while time is fixed" << std::endl;
      return 0.;
    }

  }


  void BasisFunction::unscale_parameters(boost::numeric::ublas::matrix<double> *P, double &T)
  {
    if (defined_lc && free_fun_lc) {
      std::cout << "error in robot_motion::unscale_parameters, the other unscale_parameters function should be used" << std::endl;
    }

    if (!defined_lc) {
      for (int j=0; j<nb_fun; j++) {
        for (int i=0; i<nb_P; i++) {
          (*P)(j,i) /= P_scale*funs_scale(j);
        }
      }
    } else {
      for (int j=0; j<nb_fun; j++) {
        for (int i=0; i<nb_P-6; i++) {
          (*P)(j,i) /= P_scale*funs_scale(j);
        }
      }
    }

    T /= T_scale;


  }


  void BasisFunction::unscale_parameters(boost::numeric::ublas::matrix<double> *P,
                                         boost::numeric::ublas::vector<double> &fun0,
                                         boost::numeric::ublas::vector<double> &funf, double &T)
  {
    if ((defined_lc && !free_fun_lc) || !defined_lc) {
      std::cout << "error in robot_motion::unscale_parameters, the other unscale_parameters function should be used" << std::endl;
    }

    for (int j=0; j<nb_fun; j++) {
      for (int i=0; i<nb_P-6; i++) {
        (*P)(j,i) /= P_scale*funs_scale(j);
      }
    }

    for (int j=0; j<nb_fun; j++) {
      fun0(j) /= funif_scale*funs_scale(j);
      funf(j) /= funif_scale*funs_scale(j);
    }

    T /= T_scale;

  }


  void BasisFunction::scale_parameters(boost::numeric::ublas::matrix<double> *P, double &T)
  {
    if (defined_lc && free_fun_lc) {
      std::cout << "error in robot_motion::scale_parameters, the other scale_parameters function should be used" << std::endl;
    }

    if (!defined_lc) {
      for (int j=0; j<nb_fun; j++) {
        for (int i=0; i<nb_P; i++) {
          (*P)(j,i) *= P_scale*funs_scale(j);
        }
      }
    } else {
      for (int j=0; j<nb_fun; j++) {
        for (int i=0; i<nb_P-6; i++) {
          (*P)(j,i) *= P_scale*funs_scale(j);
        }
      }
    }

    T *= T_scale;

  }


  void BasisFunction::scale_parameters(boost::numeric::ublas::matrix<double> *P,
                                       boost::numeric::ublas::vector<double> &fun0,
                                       boost::numeric::ublas::vector<double> &funf, double &T)
  {
    if ((defined_lc && !free_fun_lc) || !defined_lc) {
      std::cout << "error in robot_motion::scale_parameters, the other scale_parameters function should be used" << std::endl;
    }

    for (int j=0; j<nb_fun; j++) {
      for (int i=0; i<nb_P-6; i++) {
        (*P)(j,i) *= P_scale*funs_scale(j);
      }
    }

    for (int j=0; j<nb_fun; j++) {
      fun0(j) *= funif_scale*funs_scale(j);
      funf(j) *= funif_scale*funs_scale(j);
    }

    T *= T_scale;

  }


  double BasisFunction::unscale_hess(double value_in, int x_dep1, int x_dep2)
  {
    double value_out;

    if (!defined_lc || !free_fun_lc) { // the limit conditions are not defined
      if (grad_of_time) {
        if (x_dep1<nb_grad_uncomp-1) {
          value_out = value_in * P_scale*funs_scale(x_dep1%nb_fun);
        } else if (x_dep1<nb_grad_uncomp) {
          value_out = value_in * T_scale;
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep1 is superior to the number of parameters" << std::endl;
        }
        if (x_dep2<nb_grad_uncomp-1) {
          value_out *= P_scale*funs_scale(x_dep2%nb_fun);
        } else if (x_dep2<nb_grad_uncomp) {
          value_out *= T_scale;
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep2 is superior to the number of parameters" << std::endl;
        }
      } else {
        if (x_dep1<nb_grad_uncomp) {
          value_out = value_in * P_scale*funs_scale(x_dep1%nb_fun);
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep1 is superior to the number of parameters" << std::endl;
        }
        if (x_dep2<nb_grad_uncomp) {
          value_out *= P_scale*funs_scale(x_dep2%nb_fun);
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep2 is superior to the number of parameters" << std::endl;
        }
      }

    } else { // the limit conditions are defined, limit conditions of fun are free
      if (grad_of_time) { // gradient of time is computed
        if (x_dep1<nb_fun) {
          value_out = value_in * funif_scale*funs_scale(x_dep1);
        } else if (x_dep1<(Nint-2)*nb_fun) {
          value_out = value_in * P_scale*funs_scale(x_dep1%nb_fun);
        } else if (x_dep1<nb_grad_uncomp-1) {
          value_out = value_in * funif_scale*funs_scale(x_dep1%nb_fun);
        } else if (x_dep1<nb_grad_uncomp) {
          value_out = value_in * T_scale;
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep1 is superior to the number of parameters" << std::endl;
        }
        if (x_dep2<nb_fun) {
          value_out *= funif_scale*funs_scale(x_dep2);
        } else if (x_dep2<(Nint-2)*nb_fun) {
          value_out *= P_scale*funs_scale(x_dep2%nb_fun);
        } else if (x_dep2<nb_grad_uncomp-1) {
          value_out *= funif_scale*funs_scale(x_dep2%nb_fun);
        } else if (x_dep2<nb_grad_uncomp) {
          value_out *= T_scale;
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep2 is superior to the number of parameters" << std::endl;
        }
      } else { // gradient of time is not computed
        if (x_dep1<nb_fun) {
          value_out = value_in * funif_scale*funs_scale(x_dep1);
        } else if (x_dep1<(Nint-2)*nb_fun) {
          value_out = value_in * P_scale*funs_scale(x_dep1%nb_fun);
        } else if (x_dep1<nb_grad_uncomp) {
          value_out = value_in * funif_scale*funs_scale(x_dep1%nb_fun);
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep1 is superior to the number of parameters" << std::endl;
        }
        if (x_dep2<nb_fun) {
          value_out *= funif_scale*funs_scale(x_dep2);
        } else if (x_dep2<(Nint-2)*nb_fun) {
          value_out *= P_scale*funs_scale(x_dep2%nb_fun);
        } else if (x_dep2<nb_grad_uncomp) {
          value_out *= funif_scale*funs_scale(x_dep2%nb_fun);
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep2 is superior to the number of parameters" << std::endl;
        }
      } // end gradient of time is not computed
    } // end limit conditions are defined

    return value_out;
  }

  double BasisFunction::scale_hess(double value_in, int x_dep1, int x_dep2)
  {
    double value_out;

    if (!defined_lc || !free_fun_lc) { // the limit conditions are not defined
      if (grad_of_time) {
        if (x_dep1<nb_grad_uncomp-1) {
          value_out = value_in / (P_scale*funs_scale(x_dep1%nb_fun));
        } else if (x_dep1<nb_grad_uncomp) {
          value_out = value_in / T_scale;
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep1 is superior to the number of parameters" << std::endl;
        }
        if (x_dep2<nb_grad_uncomp-1) {
          value_out /= P_scale*funs_scale(x_dep2%nb_fun);
        } else if (x_dep2<nb_grad_uncomp) {
          value_out /= T_scale;
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep2 is superior to the number of parameters" << std::endl;
        }
      } else {
        if (x_dep1<nb_grad_uncomp) {
          value_out = value_in / (P_scale*funs_scale(x_dep1%nb_fun));
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep1 is superior to the number of parameters" << std::endl;
        }
        if (x_dep2<nb_grad_uncomp) {
          value_out /= P_scale*funs_scale(x_dep2%nb_fun);
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep2 is superior to the number of parameters" << std::endl;
        }
      }

    } else { // the limit conditions are defined, limit conditions of fun are free
      if (grad_of_time) { // gradient of time is computed
        if (x_dep1<nb_fun) {
          value_out = value_in / (funif_scale*funs_scale(x_dep1));
        } else if (x_dep1<(Nint-2)*nb_fun) {
          value_out = value_in / (P_scale*funs_scale(x_dep1%nb_fun));
        } else if (x_dep1<nb_grad_uncomp-1) {
          value_out = value_in / (funif_scale*funs_scale(x_dep1%nb_fun));
        } else if (x_dep1<nb_grad_uncomp) {
          value_out = value_in / T_scale;
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep1 is superior to the number of parameters" << std::endl;
        }
        if (x_dep2<nb_fun) {
          value_out /= funif_scale*funs_scale(x_dep2);
        } else if (x_dep2<(Nint-2)*nb_fun) {
          value_out /= P_scale*funs_scale(x_dep2%nb_fun);
        } else if (x_dep2<nb_grad_uncomp-1) {
          value_out /= funif_scale*funs_scale(x_dep2%nb_fun);
        } else if (x_dep2<nb_grad_uncomp) {
          value_out /= T_scale;
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep2 is superior to the number of parameters" << std::endl;
        }
      } else { // gradient of time is not computed
        if (x_dep1<nb_fun) {
          value_out = value_in / (funif_scale*funs_scale(x_dep1));
        } else if (x_dep1<(Nint-2)*nb_fun) {
          value_out = value_in / (P_scale*funs_scale(x_dep1%nb_fun));
        } else if (x_dep1<nb_grad_uncomp) {
          value_out = value_in / (funif_scale*funs_scale(x_dep1%nb_fun));
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep1 is superior to the number of parameters" << std::endl;
        }
        if (x_dep2<nb_fun) {
          value_out /= funif_scale*funs_scale(x_dep2);
        } else if (x_dep2<(Nint-2)*nb_fun) {
          value_out /= P_scale*funs_scale(x_dep2%nb_fun);
        } else if (x_dep2<nb_grad_uncomp) {
          value_out /= funif_scale*funs_scale(x_dep2%nb_fun);
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep2 is superior to the number of parameters" << std::endl;
        }
      } // end gradient of time is not computed
    } // end limit conditions are defined

    return value_out;
  }


  double BasisFunction::unscale_grad(double value_in, int x_dep)
  {

    if (!defined_lc || !free_fun_lc) { // the limit conditions are not defined
      if (grad_of_time) {
        if (x_dep<nb_grad_uncomp-1) {
          return value_in * P_scale*funs_scale(x_dep%nb_fun);
        } else if (x_dep<nb_grad_uncomp) {
          return value_in * T_scale;
        } else {
          std::cout << "error in BasisFunction::unscale_grad, x_dep is superior to the number of parameters" << std::endl;
        }
      } else {
        if (x_dep<nb_grad_uncomp) {
          return value_in * P_scale*funs_scale(x_dep%nb_fun);
        } else {
          std::cout << "error in BasisFunction::unscale_grad, x_dep is superior to the number of parameters" << std::endl;
        }
      }

    } else { // the limit conditions are defined, limit conditions of fun are free
      if (grad_of_time) { // gradient of time is computed
        if (x_dep<nb_fun) {
          return value_in * funif_scale*funs_scale(x_dep);
        } else if (x_dep<(Nint-2)*nb_fun) {
          return value_in * P_scale*funs_scale(x_dep%nb_fun);
        } else if (x_dep<nb_grad_uncomp-1) {
          return value_in * funif_scale*funs_scale(x_dep%nb_fun);
        } else if (x_dep<nb_grad_uncomp) {
          return value_in * T_scale;
        } else {
          std::cout << "error in BasisFunction::unscale_grad, x_dep is superior to the number of parameters" << std::endl;
        }
      } else { // gradient of time is not computed
        if (x_dep<nb_fun) {
          return value_in * funif_scale*funs_scale(x_dep);
        } else if (x_dep<(Nint-2)*nb_fun) {
          return value_in * P_scale*funs_scale(x_dep%nb_fun);
        } else if (x_dep<nb_grad_uncomp) {
          return value_in * funif_scale*funs_scale(x_dep%nb_fun);
        } else {
          std::cout << "error in BasisFunction::unscale_grad, x_dep is superior to the number of parameters" << std::endl;
        }
      } // end gradient of time is not computed
    } // end limit conditions are defined

    return 0.;
  }

  double BasisFunction::scale_grad(double value_in, int x_dep)
  {
    if (!defined_lc || !free_fun_lc) { // the limit conditions are not defined
      if (grad_of_time) {
        if (x_dep<nb_grad_uncomp-1) {
          return value_in / (P_scale*funs_scale(x_dep%nb_fun));
        } else if (x_dep<nb_grad_uncomp) {
          return value_in / T_scale;
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep is superior to the number of parameters" << std::endl;
        }
      } else {
        if (x_dep<nb_grad_uncomp) {
          return value_in / (P_scale*funs_scale(x_dep%nb_fun));
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep is superior to the number of parameters" << std::endl;
        }
      }

    } else { // the limit conditions are defined, limit conditions of fun are free
      if (grad_of_time) { // gradient of time is computed
        if (x_dep<nb_fun) {
          return value_in / (funif_scale*funs_scale(x_dep));
        } else if (x_dep<(Nint-2)*nb_fun) {
          return value_in / (P_scale*funs_scale(x_dep%nb_fun));
        } else if (x_dep<nb_grad_uncomp-1) {
          return value_in / (funif_scale*funs_scale(x_dep%nb_fun));
        } else if (x_dep<nb_grad_uncomp) {
          return value_in / T_scale;
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep is superior to the number of parameters" << std::endl;
        }
      } else { // gradient of time is not computed
        if (x_dep<nb_fun) {
          return value_in / (funif_scale*funs_scale(x_dep));
        } else if (x_dep<(Nint-2)*nb_fun) {
          return value_in / (P_scale*funs_scale(x_dep%nb_fun));
        } else if (x_dep<nb_grad_uncomp) {
          return value_in / (funif_scale*funs_scale(x_dep%nb_fun));
        } else {
          std::cout << "error in BasisFunction::unscale_hess, x_dep is superior to the number of parameters" << std::endl;
        }
      } // end gradient of time is not computed
    } // end limit conditions are defined

    return 0.;
  }


  int BasisFunction::get_nb_grad_fun() //number of gradient components for a function in compressed form
  {
    return nb_grad_fun;
  }

  int BasisFunction::get_nb_grad_comp() //number of compressed gradient components for a variable depending of all functions
  {
    return nb_grad_comp;
  }

  int BasisFunction::get_nb_grad_uncomp() //number of gradient components of any uncompressed data computed from all functions
  {
    return nb_grad_uncomp;
  }

  void BasisFunction::get_rank_grad_fun2comp(boost::numeric::ublas::matrix<int> &rank_grad_fun2comp_in)
  {
    if (nb_fun!=rank_grad_fun2comp_in.size1()) {
      std::cout << "error in BasisFunction::get_rank_grad_fun2comp, rank_grad_fun2comp_in does not have nb_fun rows" << std::endl;
    }
    if (nb_grad_fun!=rank_grad_fun2comp_in.size2()) {
      std::cout << "error in BasisFunction::get_rank_grad_fun2comp, rank_grad_fun2comp_in does not have nb_grad_fun columns" << std::endl;
    }

    for (int j=0; j<nb_fun; j++) {
      for (int i=0; i<nb_grad_fun-grad_of_time; i++) {
        rank_grad_fun2comp_in(j,i) = j+i*nb_fun;
      }
      if (grad_of_time) {
        rank_grad_fun2comp_in(j,nb_grad_fun-1) = nb_grad_comp-1;
      }
    }
  }


  void BasisFunction::joint2grad_dep(boost::numeric::ublas::vector<int> &nb_fun_dep_in,
                                     boost::numeric::ublas::matrix<int> &fun_dep_in,
                                     boost::numeric::ublas::vector<int> &nb_grad_dep_out,
                                     boost::numeric::ublas::matrix<int> &grad_dep_out)
  {
    int nb_rows = (int)nb_fun_dep_in.size();

    if (nb_rows!=fun_dep_in.size1()) {
      std::cout << "error in BasisFunction::joint2grad_dep, nb_fun_dep_in and fun_dep_in does not have the same number of rows" << std::endl;
    }
    if (nb_rows!=nb_grad_dep_out.size()) {
      std::cout << "error in BasisFunction::joint2grad_dep, nb_fun_dep_in and nb_grad_dep_out does not have the same number of rows" << std::endl;
    }
    if (nb_rows!=grad_dep_out.size1()) {
      std::cout << "error in BasisFunction::joint2grad_dep, nb_fun_dep_in and grad_dep_out does not have the same number of rows" << std::endl;
    }
    if (nb_fun!=fun_dep_in.size2()) {
      std::cout << "error in BasisFunction::joint2grad_dep, fun_dep_in does not have nb_fun columns" << std::endl;
    }
    if (nb_grad_comp!=grad_dep_out.size2()) {
      std::cout << "error in BasisFunction::joint2grad_dep, grad_dep_out does not have nb_grad_comp columns" << std::endl;
    }

    for (int i=0; i<nb_rows; i++) {
      nb_grad_dep_out(i) = nb_fun_dep_in(i)*nb_P_per_fun;
      for (int j=0; j<nb_fun_dep_in(i); j++) {
        grad_dep_out(i,4*j)=fun_dep_in(i,j);
        grad_dep_out(i,4*j+1)=fun_dep_in(i,j)+nb_fun;
        grad_dep_out(i,4*j+2)=fun_dep_in(i,j)+2*nb_fun;
        grad_dep_out(i,4*j+3)=fun_dep_in(i,j)+3*nb_fun;
      }
      if (nb_fun_dep_in(i)>0 && grad_of_time) {
        nb_grad_dep_out(i) = nb_grad_dep_out(i)+1;
        grad_dep_out(i,nb_grad_dep_out(i)-1) = nb_grad_comp-1;
      }
    }
  }

  bool BasisFunction::get_if_grad_of_time()
  {
    return grad_of_time;
  }

  void BasisFunction::calc_N(boost::numeric::ublas::vector<double> *all_t, int nb_t) //compute the vector N
  {
    int i=0, j=0;

    N[0] = 0;

    //Calcul of N
    if (nb_t < 1)
      {
        std::cout << "Error 1: At least 1 sample of time is required!" << std::endl;
      }
    else
      {
        while (j < Nint)
          {
            j++;
            //        do
            //        {
            //            i++;
            //if ((i < nb_t) && ( (*all_t)[i] <= (j-1)*Tint )) {
            //	std::cout << "error in BasisFunction::calc_N, all_t must be monotone" << std::endl;
            //}
            //        } while ((i < nb_t) && ( (*all_t)[i] <= j*Tint ));
            while ((i < nb_t) && ( (*all_t)[i] <= j*Tint )) {
              i++;
              if ((i < nb_t) && ( (*all_t)[i] < (j-1)*Tint )) {
                std::cout << "error in BasisFunction::calc_N, all_t must be monotone" << std::endl;
              }
            }
            N[j] = i;
            //i--;
          }
      }
    N[j] = nb_t; //to take into account last point
  }

  int BasisFunction::intervalle_of_time(double t) //compute the intervalle of time t
  {
    int j=1;

    while ( t > j*Tint )
      {
        j++;
      }
    return j-1;
  }

  void BasisFunction::comp_dt(boost::numeric::ublas::vector<double> &t,
                              boost::numeric::ublas::matrix<double> &data,
                              boost::numeric::ublas::matrix<double> &data_dt, int nb_t)
  {
    if (nb_t> (int) data.size2()) {
      std::cout << "error in robot_motion::comp_dt, data has less columns than nb_t" << std::endl;
    }
    if (nb_t> (int) data_dt.size2()) {
      std::cout << "error in robot_motion::comp_dt, data_dt has less columns than nb_t" << std::endl;
    }
    if (data.size1()!=data_dt.size1()) {
      std::cout << "error in robot_motion::comp_dt, data and data_dt does not have the same number of rows" << std::endl;
    }

    int i,j;
    boost::numeric::ublas::matrix<double> data_dt_cd(data.size1(),nb_t); //data_dt obtained by central difference
    double err_abs, err_rel, err;
    double max_err;
    int i_max_err, j_max_err;

    max_err = 0;
    for (i=0; i< (int) data.size1(); i++) {
      data_dt_cd(i,0) = (data(i,1)-data(i,0))/(t(1)-t(0));
      err_abs = std::fabs(data_dt(i,0)-data_dt_cd(i,0));
      if (data_dt(i,0)==0) {
        err = err_abs;
      } else {
        err_rel = err_abs/std::fabs(data_dt(i,0));
        err = fmin(err_abs, err_rel);
      }
      if (err>0.005) {
        std::cout << "big error at (" << i << ",0) data_dt " << data_dt(i,0) << " data_dt_cd " << data_dt_cd(i,0) << std::endl;
      }
      if (err>max_err) {
        max_err = err;
        i_max_err = i;
        j_max_err = 0;
      }
      for (j=1; j<nb_t-1; j++) {
        data_dt_cd(i,j) = (data(i,j+1)-data(i,j-1))/(t(j+1)-t(j-1));
        err_abs = std::fabs(data_dt(i,j)-data_dt_cd(i,j));
        if (data_dt(i,j)==0) {
          err = err_abs;
        } else {
          err_rel = err_abs/std::fabs(data_dt(i,j));
          err = fmin(err_abs, err_rel);
        }
        if (err>0.005) {
          std::cout << "big error at (" << i << "," << j << ") data_dt " << data_dt(i,j) << " data_dt_cd " << data_dt_cd(i,j) << std::endl;
        }
        if (err>max_err) {
          max_err = err;
          i_max_err = i;
          j_max_err = j;
        }
      }
      data_dt_cd(i,nb_t-1) = (data(i,nb_t-1)-data(i,nb_t-2))/(t(nb_t-1)-t(nb_t-2));
      err_abs = std::fabs(data_dt(i,nb_t-1)-data_dt_cd(i,nb_t-1));
      if (data_dt(i,0)==0) {
        err = err_abs;
      } else {
        err_rel = err_abs/std::fabs(data_dt(i,nb_t-1));
        err = fmin(err_abs, err_rel);
      }
      if (err>0.005) {
        std::cout << "big error at (" << i << "," << nb_t-1 << ") data_dt " << data_dt(i,nb_t-1) << " data_dt_cd " << data_dt_cd(i,nb_t-1) << std::endl;
      }
      if (err>max_err) {
        max_err = err;
        i_max_err = i;
        j_max_err = nb_t-1;
      }
    }

    std::cout << "dimensions of data (" << data.size1() << "," << data.size2() << ")" << std::endl;
    std::cout << "error max " << max_err << " at (" << i_max_err << "," << j_max_err << ")" << std::endl;

    //plot
    double y_min, y_max;
    double *xval, *yval;
    xval = new double [nb_t];
    yval = new double [nb_t];

    //using namespace System;

    metafl("XWIN\0"); //define output
    setpag("da4l");   //define page style

    disini();  //initialise
    pagera();  //define a border
    complx();  //define the font
    axspos(450,1800); //define the axis position
    axslen(2200,1200); //define the axis length

    name("time","x");  //define name of x axis
    name("joint angles","y"); //define name of y axis

    labdig(1,"x"); //define number of digits in labels
    ticks(10,"xy"); //number of ticks in axis label

    titlin("Joint datas",1);

    y_min = 1e20;
    y_max = -1e20;
    //i=0;
    for (int i=0; i< (int) data.size1(); i++) {
      for (int j=0; j<nb_t; j++) {
        y_min = fmin(y_min,data_dt(i,j));
        y_min = fmin(y_min,data_dt_cd(i,j));
        y_max = fmax(y_max,data_dt(i,j));
        y_max = fmax(y_max,data_dt_cd(i,j));
      }
    }

    graf(0.,T,0.,T/nb_P,y_min,y_max,y_min,(y_max-y_min)/10);
    title();

    //color("red"); //define color for next plot
    chncrv("both");

    for (int j=0; j<nb_t; j++) {
      xval[j] = t(j);
    }

    //i=0;
    for (int i=0; i< (int) data.size1(); i++) {
      for (int j=0; j<nb_t; j++) {
        yval[j] = data_dt(i,j);
      }
      curve(xval,yval,nb_t); //plot curve
      for (int j=0; j<nb_t; j++) {
        yval[j] = data_dt_cd(i,j);
      }
      curve(xval,yval,nb_t); //plot curve
    }

    color("fore");
    dash();   // define dash style for next plots
    xaxgit(); //plot y=0 line

    endgrf();

    disfin();

    delete[] xval;
    delete[] yval;

  }

} // end of namespace roboptim.
