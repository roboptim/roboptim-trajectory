// Copyright (C) 2009 by Sylvain Miossec, Thomas Moulard, AIST, CNRS, INRIA.
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

#ifndef ROBOPTIM_TRAJECTORY_BASIS_FUNCTION_HH
# define ROBOPTIM_TRAJECTORY_BASIS_FUNCTION_HH

# include <boost/numeric/ublas/vector.hpp>
# include <boost/numeric/ublas/io.hpp>
# include <boost/numeric/ublas/matrix.hpp>
# include <boost/numeric/ublas/matrix_proxy.hpp>

//# include "real64/dislin.h"

namespace roboptim
{
  namespace ublas = boost::numeric::ublas;

  class BasisFunction
  {
  private:
    BasisFunction ();

    void define_sizes ();

  public:
    BasisFunction (int nb_P_per_fun_in);

    BasisFunction (int nb_fun_in, int nb_P_in, int nb_t_m_in, bool grad_of_time_in, int nb_P_per_fun_in);

    BasisFunction (int nb_fun_in, int nb_P_in, int nb_t_m_in, bool grad_of_time_in, bool defined_lc_in, bool free_fun_lc_in, int nb_P_per_fun_in);

    virtual ~BasisFunction ();

    void define_sizes (int nb_fun_in, int nb_P_in, int nb_t_m_in, bool grad_of_time_in);

    void define_sizes (int nb_fun_in, int nb_P_in, int nb_t_m_in, bool grad_of_time_in, bool defined_lc_in, bool free_fun_lc_in);

    void copy (const BasisFunction& bf);

    int get_nb_call_def_sizes (); //give the number of call to define_sizes, so that other class getting sizes from this class can track if sizes have changed

    void def_parameters (ublas::matrix<double> *P_in, double T_in);

    void def_parameters (ublas::matrix<double> *P_in,
                         const ublas::vector<double> &fun0,
                         const ublas::vector<double> &funf, double T_in);

    inline void check_T (double T_in)
    {
      if (T_in<0) {
        std::cout << "error in class basis_dunction, duration T must be positive, T = " << T_in << std::endl;
      }
    }

    void define_fixed_lc (
                         const ublas::vector<double> &fun0,
                         const ublas::vector<double> &funf);

    void get_fixed_lc (ublas::vector<double> &fun0,
                      ublas::vector<double> &funf) const;

    void define_fixed_time (double T_in);

    void get_fixed_time (double &T_out);

    double get_free_time ();

    bool get_defined_lc ();

    bool get_free_fun_lc ();

    int get_nb_fun ();

    int get_nb_P_free (); //number of free parameters for each joint (if limit conditions are defined it is less than nb_P)

    int get_Nint ();

    void save_parameters (const char* file_name, ublas::matrix<double> *P_in, double T_in);

    void save_parameters (
                         const char* file_name,
                         ublas::matrix<double> *P_in,
                         const ublas::vector<double> &fun0,
                         const ublas::vector<double> &funf, double T_in);

    void load_parameters (const char* file_name, ublas::matrix<double> *P_out, double &T_out);

    void load_parameters (const char* file_name,
                         ublas::matrix<double> *P_out,
                         ublas::vector<double> &fun0_out,
                         ublas::vector<double> &funf_out, double &T_out);

    // read all data of a motion, and save all necessary data inside basis_function class for direct trajectory computation
    void load_parameters_and_sizes (const char* file_name,
                                   ublas::matrix<double> *P_out,
                                   ublas::vector<double> &fun0_out,
                                   ublas::vector<double> &funf_out, double &T_out);

    // read all data of a motion, and save all necessary data inside basis_function class for direct trajectory computation
    void load_parameters_and_sizes (const char* file_name);

    void define_parameters_scale (double P_scale, double T_scale, double funif_scale,
                                 const ublas::vector<double> & funs_scale);

    void get_parameters_scale (double &P_scale, double &T_scale, double &funif_scale,
                              ublas::vector<double> &funs_scale);

    void get_parameters_scale (ublas::vector<double> &x_scale_out);

    void get_scale_for_auto_scale (int &nb_x_scale_out,
                                  ublas::vector<int> &nb_p_x_scale_p_out,
                                  ublas::matrix<int> &p_x_scale_p_out);

    void set_scale_from_auto_scale (ublas::vector<double> &p_x_scale_in);

    void convert_parameters_P2x (const ublas::matrix<double> *P_in, const double T_in, double* x_out);

    void convert_parameters_P2x (const ublas::matrix<double> *P_in,
                                const ublas::vector<double> &fun0,
                                const ublas::vector<double> &funf, const double T_in, double* x_out);

    void convert_parameters_P2x (const ublas::matrix<bool> *P_in, const bool T_in, bool* x_out);

    void convert_parameters_P2x (const ublas::matrix<bool> *P_in,
                                const ublas::vector<bool> &fun0,
                                const ublas::vector<bool> &funf, const bool T_in, bool* x_out);

    void convert_parameters_x2P (const double* x_in, ublas::matrix<double> *P_out, double &T_out);

    void convert_parameters_x2P (const double* x_in,
                                ublas::matrix<double> *P_out,
                                ublas::vector<double> &fun0_out,
                                ublas::vector<double> &funf_out, double &T_out);


    void change_P_parameter_in_x (ublas::vector<double> &x_in_out, int i, int j, double P_ij_in);

    void change_fun0_parameter_in_x (ublas::vector<double> &x_in_out, int i, double qi_i_in);

    void change_funf_parameter_in_x (ublas::vector<double> &x_in_out, int i, double qf_i_in);

    void change_T_parameter_in_x (ublas::vector<double> &x_in_out, double T_in);


    void change_P_parameter_in_x (ublas::vector<bool> &x_in_out, int i, int j, bool P_ij_in);

    void change_fun0_parameter_in_x (ublas::vector<bool> &x_in_out, int i, bool qi_i_in);

    void change_funf_parameter_in_x (ublas::vector<bool> &x_in_out, int i, bool qf_i_in);

    void change_T_parameter_in_x (ublas::vector<bool> &x_in_out, bool T_in);


    double get_P_parameter_in_x (const ublas::vector<double> &x_in, int i, int j) const;

    double get_fun0_parameter_in_x (const ublas::vector<double> &x_in, int i) const;

    double get_funf_parameter_in_x (const ublas::vector<double> &x_in, int i) const;

    double get_T_parameter_in_x (const ublas::vector<double> &x_in) const;


    bool get_P_parameter_in_x (const ublas::vector<bool> &x_in, int i, int j) const;

    bool get_fun0_parameter_in_x (const ublas::vector<bool> &x_in, int i)const;

    bool get_funf_parameter_in_x (const ublas::vector<bool> &x_in, int i) const;

    bool get_T_parameter_in_x (const ublas::vector<bool> &x_in)const ;


    void unscale_parameters (ublas::matrix<double> *P, double &T);

    void unscale_parameters (ublas::matrix<double> *P,
                            ublas::vector<double> &fun0,
                            ublas::vector<double> &funf, double &T);

    void scale_parameters (ublas::matrix<double> *P, double &T);

    void scale_parameters (ublas::matrix<double> *P,
                          ublas::vector<double> &fun0,
                          ublas::vector<double> &funf, double &T);

    double unscale_hess (double value_in, int x_dep1, int x_dep2);

    double scale_hess (double value_in, int x_dep1, int x_dep2);

    double unscale_grad (double value_in, int x_dep1);

    double scale_grad (double value_in, int x_dep1);

    int get_nb_grad_fun (); //number of gradient components for a function in compressed form
    int get_nb_grad_comp (); //number of compressed gradient components for a variable depending of all functions
    int get_nb_grad_uncomp (); //number of gradient components of any uncompressed data computed from all functions

    void get_rank_grad_fun2comp (ublas::matrix<int> &rank_grad_fun2comp_in);

    virtual void get_rank_nz_uncomp (ublas::vector<int> &time_points_in,
                                    ublas::vector<int> &funorcomp_type_in,
                                    ublas::vector<int> &nb_nz_per_row,
                                    ublas::matrix<int> &rank_grad_u2uncomp_out, int nb_t) = 0;

    void joint2grad_dep (ublas::vector<int> &nb_fun_dep_in,
                        ublas::matrix<int> &fun_dep_in,
                        ublas::vector<int> &nb_grad_dep_out,
                        ublas::matrix<int> &grad_dep_out);

    bool get_if_grad_of_time ();

    virtual void calc_fun (ublas::vector<double> *all_t, ublas::matrix<double> *fun, int nb_t) = 0;
    virtual void calc_dfun (ublas::vector<double> *all_t, ublas::matrix<double> *dfun, int nb_t) = 0;
    virtual void calc_ddfun (ublas::vector<double> *all_t, ublas::matrix<double> *ddfun, int nb_t) = 0;

    virtual void calc_fun_grad (ublas::vector<double> *all_t,
                               ublas::matrix<double> *fun,
                               ublas::matrix< ublas::vector<double> > *fun_grad, int nb_t) = 0;

    virtual void calc_dfun_grad (ublas::vector<double> *all_t,
                                ublas::matrix<double> *dfun,
                                ublas::matrix< ublas::vector<double> > *dfun_grad, int nb_t) = 0;
    virtual void calc_ddfun_grad (ublas::vector<double> *all_t,
                                 ublas::matrix<double> *ddfun,
                                 ublas::matrix< ublas::vector<double> > *ddfun_grad, int nb_t) = 0;

    virtual void uncompress_grad (ublas::vector<double> *all_t,
                                 ublas::matrix< ublas::vector<double> > *data_grad_comp,
                                 ublas::matrix< ublas::vector<double> > *data_grad_uncomp, int nb_t) = 0;

    void comp_dt (ublas::vector<double> &t,
                 ublas::matrix<double> &data,
                 ublas::matrix<double> &data_dt, int nb_t);


  protected:

    //variables
    int nb_fun;
    int nb_P;
    int nb_P_per_fun; //number of parameters of P from which the value of a function depend whatever the time (must be defined in the constructor of the child class)
    int nb_grad_fun; //number of compressed gradient components of each function returned by class basis_function
    int nb_grad_comp; //number of compressed gradient components for a variable depending of all functions
    int nb_grad_uncomp; //number of gradient components of any uncompressed data computed from all functions
    int nb_call_def_sizes; //number of call to define_sizes, so that other class getting sizes from this class can track if sizes have changed
    double T; //superior limit of time values
    double Tint; //duration of an intervalle
    int Nint; //number of intervalles
    int nb_t_m; //Number max of sample in the range of time
    ublas::matrix<double> *matP; //Matrix of splines coefficients for each actionnors
    ublas::vector<double> *all_tn; //Vector of samples of time, normalized
    ublas::vector<int> N; //N[j] number of values of time in [0,j*Tint]
    ublas::vector<double> fun0_fixed; //initial value of functions in the fixed initial/final case
    ublas::vector<double> funf_fixed; //initial value of functions in the fixed initial/final case
    double _T_fixed; //saved value of fixed time

    bool grad_of_time; //if the gradient of time is computed
    bool defined_lc; //1 : the limit conditions of fun are given and limit conditions of dfun and ddfun are zero
    bool free_fun_lc;  //the limit conditions of fun are free (1) or fixed (0)
    // (the difference is for the gradient and where the limit values of fun are given)
    bool given_lc; //in the case of fixed limit conditions, if the limit conditions have been defined
    bool given_T; //in the case of fixed time, if the time have been defined

    //parameters for the scaling
    double funif_scale;
    double P_scale;
    double T_scale;
    ublas::vector<double> funs_scale;

    virtual void def_Nint() = 0; //must define Nint from nb_P

    virtual void solve_P_0(const ublas::vector<double> &fun0) = 0;

    virtual void solve_P_f(const ublas::vector<double> &funf) = 0;

    void calc_N(ublas::vector<double> *all_t, int nb_t); //compute the vector N

    int intervalle_of_time(double t); //give the intervalle of time t
  };

} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_BASIS_FUNCTION_HH
