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

#ifndef ROBOPTIM_TRAJECTORY_BSPLINE_HH
# define ROBOPTIM_TRAJECTORY_BSPLINE_HH

# include <boost/numeric/ublas/vector.hpp>
# include <boost/numeric/ublas/io.hpp>
# include <boost/numeric/ublas/matrix.hpp>
# include <boost/numeric/ublas/matrix_proxy.hpp>

# include <roboptim-trajectory/basis-function.hh>

namespace roboptim
{
  namespace ublas = boost::numeric::ublas;

  class BSpline : public BasisFunction
  {
  public:
    BSpline ();
    BSpline (int nb_fun_in, int nb_P_in, int nb_t_m_in, bool grad_of_time_in);
    BSpline (int nb_fun_in, int nb_P_in, int nb_t_m_in, bool grad_of_time_in, bool defined_lc_in, bool free_fun_lc_in);

    virtual ~BSpline ();

    void calc_fun (ublas::vector<double> *all_t, ublas::matrix<double> *fun, int nb_t);
    void calc_dfun (ublas::vector<double> *all_t, ublas::matrix<double> *dfun, int nb_t);
    void calc_ddfun (ublas::vector<double> *all_t, ublas::matrix<double> *ddfun, int nb_t);

    void calc_fun (double t, ublas::vector<double> *fun);
    void calc_dfun (double t, ublas::vector<double> *dfun);
    void calc_ddfun (double t, ublas::vector<double> *ddfun);

    void calc_fun_grad (ublas::vector<double> *all_t,
                       ublas::matrix<double> *fun,
                       ublas::matrix< ublas::vector<double> > *fun_grad, int nb_t);
    void calc_dfun_grad (ublas::vector<double> *all_t,
                        ublas::matrix<double> *dfun,
                        ublas::matrix< ublas::vector<double> > *dfun_grad, int nb_t);
    void calc_ddfun_grad (ublas::vector<double> *all_t,
                         ublas::matrix<double> *ddfun,
                         ublas::matrix< ublas::vector<double> > *ddfun_grad, int nb_t);

    void get_rank_nz_uncomp (ublas::vector<int> &time_points_in,
                            ublas::vector<int> &funorcomp_type_in,
                            ublas::vector<int> &nb_nz_per_row_out,
                            ublas::matrix<int> &rank_grad_nz_out, int nb_t);

    void grad_comp2sparse (ublas::vector<int> &time_points_in,
                          ublas::vector<int> &funorcomp_type_in,
                          ublas::matrix<double> &data_grad, int nb_t);

    void uncompress_grad (ublas::vector<double> *all_t,
                         ublas::matrix< ublas::vector<double> > *data_grad_comp,
                         ublas::matrix< ublas::vector<double> > *data_grad_uncomp, int nb_t);

    void uncompress_grad_int (ublas::vector<double> *all_t,
                             ublas::vector<double> *data,
                             ublas::matrix<double> *data_grad_comp,
                             ublas::vector<double> *data_grad_uncomp_sum, int nb_t, double Tfinal);

    template<class Type>
    void uncompress_grad (double t,
                         ublas::matrix< Type > *data_grad_comp,
                         ublas::matrix< Type > *data_grad_uncomp)
    {
      if ((*data_grad_uncomp).size1()!=(*data_grad_comp).size1()) {
        std::cout << "error in bspline::uncompress_grad, data_grad_comp and data_grad_uncomp does not have the same number of rows" << std::endl;
      }
      if (nb_grad_comp!=(*data_grad_comp).size2()) {
        std::cout << "error in bspline::uncompress_grad, data_grad_comp does not have nb_grad_comp components of gradient" << std::endl;
      }
      if (nb_grad_uncomp!=(*data_grad_uncomp).size2()) {
        std::cout << "error in bspline::uncompress_grad, data_grad_uncomp does not have nb_grad_uncomp components of gradient" << std::endl;
        std::cout << "nb_grad_uncomp " << nb_grad_uncomp << std::endl;
        std::cout << "number of components of gradient of data_grad_uncomp " << (*data_grad_uncomp).size2() << std::endl;
      }

      unsigned int i;
      int k, kk;

      k = intervalle_of_time(t);


      if (!defined_lc) { // the limit conditions are not defined

        if (grad_of_time) {
          for (i=0; i<(*data_grad_comp).size1(); i++) {
            for (kk=0; kk<k*nb_fun; kk++) {
              (*data_grad_uncomp)(i,kk) = 0;
            }
            for (kk=0; kk<nb_grad_comp-1; kk++) {
              (*data_grad_uncomp)(i,kk+k*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
            }
            for (kk=k*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
              (*data_grad_uncomp)(i,kk) = 0;
            }
            (*data_grad_uncomp)(i,nb_grad_uncomp-1) = (*data_grad_comp)(i,nb_grad_comp-1)/T_scale;
          }
        } else {
          for (i=0; i<(*data_grad_comp).size1(); i++) {
            for (kk=0; kk<k*nb_fun; kk++) {
              (*data_grad_uncomp)(i,kk) = 0;
            }
            for (kk=0; kk<nb_grad_comp; kk++) {
              (*data_grad_uncomp)(i,kk+k*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
            }
            for (kk=k*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) {
              (*data_grad_uncomp)(i,kk) = 0;
            }
          }
        }

        //---------------------------------------------------------------------------
      } else if (defined_lc & !free_fun_lc) { // the limit conditions are completely defined

        if (grad_of_time) { // gradient of time is computed

          if (k<3)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(k+1)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = (*data_grad_comp)(i,kk+(3-k)*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=(k+1)*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  (*data_grad_uncomp)(i,nb_grad_uncomp-1) = (*data_grad_comp)(i,nb_grad_comp-1)/T_scale;
                }
            }
          else if (k<Nint-3)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(k-3)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  for (kk=0; kk<nb_grad_comp-1; kk++) {
                    (*data_grad_uncomp)(i,kk+(k-3)*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=(k-3)*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  (*data_grad_uncomp)(i,nb_grad_uncomp-1) = (*data_grad_comp)(i,nb_grad_comp-1)/T_scale;
                }
            }
          else
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(k-3)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  for (kk=0; kk<(Nint-k)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(k-3)*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  (*data_grad_uncomp)(i,nb_grad_uncomp-1) = (*data_grad_comp)(i,nb_grad_comp-1)/T_scale;
                }
            }
          //---------------------------------------------------------------------------
        } else { // gradient of time is not computed
          if (k<3)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(k+1)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = (*data_grad_comp)(i,kk+(3-k)*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=(k+1)*nb_fun; kk<nb_grad_uncomp; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                }
            }
          else if (k<Nint-3)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(k-3)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  for (kk=0; kk<nb_grad_comp; kk++) {
                    (*data_grad_uncomp)(i,kk+(k-3)*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=(k-3)*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                }
            }
          else
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(k-3)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  for (kk=0; kk<(Nint-k)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(k-3)*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
                  }
                }
            }
        } // end gradient of time is not computed

        //---------------------------------------------------------------------------
      } else { // the limit conditions are defined, limit conditions of fun are free
        calc_P_fun_0();
        calc_P_fun_f();

        if (grad_of_time) { // gradient of time is computed
          if (k==0)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = ((*data_grad_comp)(i,kk)*P0_fun
                                                 + (*data_grad_comp)(i,kk+nb_fun)*P1_fun
                                                 + (*data_grad_comp)(i,kk+2*nb_fun)*P2_fun)/(funif_scale*funs_scale(kk%nb_fun));
                    (*data_grad_uncomp)(i,kk+nb_fun) = (*data_grad_comp)(i,kk+3*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=2*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  (*data_grad_uncomp)(i,nb_grad_uncomp-1) = (*data_grad_comp)(i,nb_grad_comp-1)/T_scale;
                }
            }
          else if (k==1)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = ((*data_grad_comp)(i,kk)*P1_fun
                                                 + (*data_grad_comp)(i,kk+nb_fun)*P2_fun)/(funif_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=0; kk<2*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+nb_fun) = (*data_grad_comp)(i,kk+2*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=3*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  (*data_grad_uncomp)(i,nb_grad_uncomp-1) = (*data_grad_comp)(i,nb_grad_comp-1)/T_scale;
                }
            }
          else if (k==2)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = (*data_grad_comp)(i,kk)*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=0; kk<3*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+nb_fun) = (*data_grad_comp)(i,kk+nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=4*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  (*data_grad_uncomp)(i,nb_grad_uncomp-1) = (*data_grad_comp)(i,nb_grad_comp-1)/T_scale;
                }
            }
          else if (k<Nint-3)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(k-2)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  for (kk=0; kk<nb_grad_comp-1; kk++) {
                    (*data_grad_uncomp)(i,kk+(k-2)*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=(k-2)*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  (*data_grad_uncomp)(i,nb_grad_uncomp-1) = (*data_grad_comp)(i,nb_grad_comp-1)/T_scale;
                }
            }
          else if (k==Nint-3)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  for (kk=0; kk<3*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(Nint-5)*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(Nint-2)*nb_fun) = (*data_grad_comp)(i,kk+3*nb_fun)*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                  }
                  (*data_grad_uncomp)(i,nb_grad_uncomp-1) = (*data_grad_comp)(i,nb_grad_comp-1)/T_scale;
                }
            }
          else if (k==Nint-2)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  for (kk=0; kk<2*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(Nint-4)*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(Nint-2)*nb_fun) = ((*data_grad_comp)(i,kk+2*nb_fun)*Pendm2_fun
                                                                 + (*data_grad_comp)(i,kk+3*nb_fun)*Pendm1_fun)/(funif_scale*funs_scale(kk%nb_fun));
                  }
                  (*data_grad_uncomp)(i,nb_grad_uncomp-1) = (*data_grad_comp)(i,nb_grad_comp-1)/T_scale;
                }
            }
          else if (k==Nint-1)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(Nint-3)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(Nint-3)*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(Nint-2)*nb_fun) = ((*data_grad_comp)(i,kk+1*nb_fun)*Pendm2_fun
                                                                 + (*data_grad_comp)(i,kk+2*nb_fun)*Pendm1_fun
                                                                 + (*data_grad_comp)(i,kk+3*nb_fun)*Pend_fun)/(funif_scale*funs_scale(kk%nb_fun));
                  }
                  (*data_grad_uncomp)(i,nb_grad_uncomp-1) = (*data_grad_comp)(i,nb_grad_comp-1)/T_scale;
                }
            }
          //---------------------------------------------------------------------------
        } else { // gradient of time is not computed
          if (k==0)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = ((*data_grad_comp)(i,kk)*P0_fun
                                                 + (*data_grad_comp)(i,kk+nb_fun)*P1_fun
                                                 + (*data_grad_comp)(i,kk+2*nb_fun)*P2_fun)/(funif_scale*funs_scale(kk%nb_fun));
                    (*data_grad_uncomp)(i,kk+nb_fun) = (*data_grad_comp)(i,kk+3*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=2*nb_fun; kk<nb_grad_uncomp; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                }
            }
          else if (k==1)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = ((*data_grad_comp)(i,kk)*P1_fun
                                                 + (*data_grad_comp)(i,kk+nb_fun)*P2_fun)/(funif_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=0; kk<2*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+nb_fun) = (*data_grad_comp)(i,kk+2*nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=3*nb_fun; kk<nb_grad_uncomp; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                }
            }
          else if (k==2)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = (*data_grad_comp)(i,kk)*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=0; kk<3*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+nb_fun) = (*data_grad_comp)(i,kk+nb_fun)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=4*nb_fun; kk<nb_grad_uncomp; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                }
            }
          else if (k<Nint-3)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(k-2)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  for (kk=0; kk<nb_grad_comp; kk++) {
                    (*data_grad_uncomp)(i,kk+(k-2)*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=(k-2)*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                }
            }
          else if (k==Nint-3)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  for (kk=0; kk<3*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(Nint-5)*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(Nint-2)*nb_fun) = (*data_grad_comp)(i,kk+3*nb_fun)*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                  }
                }
            }
          else if (k==Nint-2)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  for (kk=0; kk<2*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(Nint-4)*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(Nint-2)*nb_fun) = ((*data_grad_comp)(i,kk+2*nb_fun)*Pendm2_fun
                                                                 + (*data_grad_comp)(i,kk+3*nb_fun)*Pendm1_fun)/(funif_scale*funs_scale(kk%nb_fun));
                  }
                }
            }
          else if (k==Nint-1)
            {
              for (i=0; i<(*data_grad_comp).size1(); i++)
                {
                  for (kk=0; kk<(Nint-3)*nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk) = 0;
                  }
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(Nint-3)*nb_fun) = (*data_grad_comp)(i,kk)/(P_scale*funs_scale(kk%nb_fun));
                  }
                  for (kk=0; kk<nb_fun; kk++) {
                    (*data_grad_uncomp)(i,kk+(Nint-2)*nb_fun) = ((*data_grad_comp)(i,kk+1*nb_fun)*Pendm2_fun
                                                                 + (*data_grad_comp)(i,kk+2*nb_fun)*Pendm1_fun
                                                                 + (*data_grad_comp)(i,kk+3*nb_fun)*Pend_fun)/(funif_scale*funs_scale(kk%nb_fun));
                  }
                }
            }
        } // end gradient of time is not computed
      } // end limit conditions are defined


    }


    template<class Type>
    void uncompress_grad_it (ublas::vector<double> *all_t,
                            typename ublas::matrix< Type >::iterator1 &data_grad_comp_rit,
                            typename ublas::matrix< Type >::iterator1 &data_grad_uncomp_rit, int nb_t)
    {
      if (nb_t!=data_grad_comp_rit().size1()) {
        std::cout << "error in bspline::uncompress_grad, data_grad_comp does not have nb_t rows" << std::endl;
      }
      if (nb_t!=data_grad_uncomp_rit().size1()) {
        std::cout << "error in bspline::uncompress_grad, data_grad_uncomp does not have nb_t rows" << std::endl;
      }
      if (nb_grad_comp!=data_grad_comp_rit().size2()) {
        std::cout << "error in bspline::uncompress_grad, data_grad_comp does not have nb_grad_comp components of gradient" << std::endl;
      }
      if (nb_grad_uncomp!=data_grad_uncomp_rit().size2()) {
        std::cout << "error in bspline::uncompress_grad, data_grad_uncomp does not have nb_grad_uncomp components of gradient" << std::endl;
        std::cout << "nb_grad_uncomp " << nb_grad_uncomp << std::endl;
        std::cout << "number of components of gradient of data_grad_uncomp " << data_grad_uncomp_rit().size2() << std::endl;
      }

      int j, k, kk;
      typename ublas::matrix< Type >::iterator2 data_grad_comp_cit;
      typename ublas::matrix< Type >::iterator2 data_grad_uncomp_cit;

      //computation of vector N such that N[j] is the number of time sample in [0,j*Tint]
      calc_N(all_t, nb_t);


      if (!defined_lc) { // the limit conditions are not defined

        if (grad_of_time) {
          for (k=0 ; k<Nint ; k++)
            {
              for (j=N[k] ; j<N[k+1] ; j++)
                {
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  for (kk=0; kk<k*nb_fun; kk++) {
                    *data_grad_uncomp_cit = 0;
                    ++data_grad_uncomp_cit;
                  }
                  for (kk=0; kk<nb_grad_comp-1; kk++) {
                    *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                    ++data_grad_comp_cit;
                    ++data_grad_uncomp_cit;
                  }
                  for (kk=k*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
                    *data_grad_uncomp_cit = 0;
                    ++data_grad_uncomp_cit;
                  }
                  *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
            }
        } else {
          for (k=0 ; k<Nint ; k++)
            {
              for (j=N[k] ; j<N[k+1] ; j++)
                {
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  for (kk=0; kk<k*nb_fun; kk++) {
                    *data_grad_uncomp_cit = 0;
                    ++data_grad_uncomp_cit;
                  }
                  for (kk=0; kk<nb_grad_comp; kk++) {
                    *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                    ++data_grad_comp_cit;
                    ++data_grad_uncomp_cit;
                  }
                  for (kk=k*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) {
                    *data_grad_uncomp_cit = 0;
                    ++data_grad_uncomp_cit;
                  }
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
            }
        }

        //---------------------------------------------------------------------------
      } else if (defined_lc & !free_fun_lc) { // the limit conditions are completely defined

        if (grad_of_time) { // gradient of time is computed
          for (j=N[0] ; j<N[1] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              data_grad_comp_cit += 3*nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=nb_fun; kk<nb_grad_uncomp-1; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[1] ; j<N[2] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              data_grad_comp_cit += 2*nb_fun;
              for (kk=0; kk<2*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=2*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[2] ; j<N[3] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              data_grad_comp_cit += nb_fun;
              for (kk=0; kk<3*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=3*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (k=3 ; k<Nint-3 ; k++)
            {
              for (j=N[k] ; j<N[k+1] ; j++)
                {
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  for (kk=0; kk<(k-3)*nb_fun; kk++) {
                    *data_grad_uncomp_cit = 0;
                    ++data_grad_uncomp_cit;
                  }
                  for (kk=0; kk<nb_grad_comp-1; kk++) {
                    *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                    ++data_grad_comp_cit;
                    ++data_grad_uncomp_cit;
                  }
                  for (kk=(k-3)*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
                    *data_grad_uncomp_cit = 0;
                    ++data_grad_uncomp_cit;
                  }
                  *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
            }
          for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<(Nint-6)*nb_fun; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<3*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_comp_cit += nb_fun;
              *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<2*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_comp_cit += 2*nb_fun;
              *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[Nint-1] ; j<N[Nint] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_comp_cit += 3*nb_fun;
              *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          //---------------------------------------------------------------------------
        } else { // gradient of time is not computed
          for (j=N[0] ; j<N[1] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              data_grad_comp_cit += 3*nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=nb_fun; kk<nb_grad_uncomp; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[1] ; j<N[2] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              data_grad_comp_cit += 2*nb_fun;
              for (kk=0; kk<2*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=2*nb_fun; kk<nb_grad_uncomp; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[2] ; j<N[3] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              data_grad_comp_cit += nb_fun;
              for (kk=0; kk<3*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=3*nb_fun; kk<nb_grad_uncomp; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (k=3 ; k<Nint-3 ; k++)
            {
              for (j=N[k] ; j<N[k+1] ; j++)
                {
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  for (kk=0; kk<(k-3)*nb_fun; kk++) {
                    *data_grad_uncomp_cit = 0;
                    ++data_grad_uncomp_cit;
                  }
                  for (kk=0; kk<nb_grad_comp; kk++) {
                    *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                    ++data_grad_comp_cit;
                    ++data_grad_uncomp_cit;
                  }
                  for (kk=(k-3)*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) {
                    *data_grad_uncomp_cit = 0;
                    ++data_grad_uncomp_cit;
                  }
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
            }
          for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<(Nint-6)*nb_fun; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<3*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<2*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[Nint-1] ; j<N[Nint] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
        } // end gradient of time is not computed

        //---------------------------------------------------------------------------
      } else { // the limit conditions are defined, limit conditions of fun are free
        calc_P_fun_0();
        calc_P_fun_f();

        if (grad_of_time) { // gradient of time is computed
          for (j=N[0] ; j<N[1] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit*P0_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_uncomp_cit -= nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit += *data_grad_comp_cit*P1_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_uncomp_cit -= nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit += *data_grad_comp_cit*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=2*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[1] ; j<N[2] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit*P1_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_uncomp_cit -= nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit += *data_grad_comp_cit*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<2*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=3*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[2] ; j<N[3] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<3*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=4*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (k=3 ; k<Nint-3 ; k++)
            {
              for (j=N[k] ; j<N[k+1] ; j++)
                {
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  for (kk=0; kk<(k-2)*nb_fun; kk++) {
                    *data_grad_uncomp_cit = 0;
                    ++data_grad_uncomp_cit;
                  }
                  for (kk=0; kk<nb_grad_comp-1; kk++) {
                    *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                    ++data_grad_comp_cit;
                    ++data_grad_uncomp_cit;
                  }
                  for (kk=(k-2)*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
                    *data_grad_uncomp_cit = 0;
                    ++data_grad_uncomp_cit;
                  }
                  *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
            }
          for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<3*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<2*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_uncomp_cit -= nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit += *data_grad_comp_cit*Pendm1_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[Nint-1] ; j<N[Nint] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<(Nint-3)*nb_fun; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_uncomp_cit -= nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit += *data_grad_comp_cit*Pendm1_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_uncomp_cit -= nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit += *data_grad_comp_cit*Pend_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          //---------------------------------------------------------------------------
        } else { // gradient of time is not computed
          for (j=N[0] ; j<N[1] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit*P0_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_uncomp_cit -= nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit += *data_grad_comp_cit*P1_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_uncomp_cit -= nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit += *data_grad_comp_cit*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=2*nb_fun; kk<nb_grad_uncomp; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[1] ; j<N[2] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit*P1_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_uncomp_cit -= nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit += *data_grad_comp_cit*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<2*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=3*nb_fun; kk<nb_grad_uncomp; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[2] ; j<N[3] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<3*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=4*nb_fun; kk<nb_grad_uncomp; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (k=3 ; k<Nint-3 ; k++)
            {
              for (j=N[k] ; j<N[k+1] ; j++)
                {
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  for (kk=0; kk<(k-2)*nb_fun; kk++) {
                    *data_grad_uncomp_cit = 0;
                    ++data_grad_uncomp_cit;
                  }
                  for (kk=0; kk<nb_grad_comp; kk++) {
                    *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                    ++data_grad_comp_cit;
                    ++data_grad_uncomp_cit;
                  }
                  for (kk=(k-2)*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) {
                    *data_grad_uncomp_cit = 0;
                    ++data_grad_uncomp_cit;
                  }
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
            }
          for (j=N[Nint-3] ; j<N[Nint-2] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<3*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[Nint-2] ; j<N[Nint-1] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<2*nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_uncomp_cit -= nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit += *data_grad_comp_cit*Pendm1_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
          for (j=N[Nint-1] ; j<N[Nint] ; j++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<(Nint-3)*nb_fun; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_uncomp_cit -= nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit += *data_grad_comp_cit*Pendm1_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              data_grad_uncomp_cit -= nb_fun;
              for (kk=0; kk<nb_fun; kk++) {
                *data_grad_uncomp_cit += *data_grad_comp_cit*Pend_fun/(funif_scale*funs_scale(kk%nb_fun));
                ++data_grad_comp_cit;
                ++data_grad_uncomp_cit;
              }
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
        } // end gradient of time is not computed
      } // end limit conditions are defined

    }


    template<class Type>
    void uncompress_grad(ublas::vector<double> *all_t,
                         ublas::matrix< Type > *data_grad_comp,
                         ublas::matrix< Type > *data_grad_uncomp, int nb_t)
    {
      typename ublas::matrix< Type >::iterator1 data_grad_comp_rit=(*data_grad_comp).begin1();
      typename ublas::matrix< Type >::iterator1 data_grad_uncomp_rit=(*data_grad_uncomp).begin1();
      uncompress_grad_it<Type>(all_t, data_grad_comp_rit, data_grad_uncomp_rit, nb_t);
    }

    template<class Type>
    void uncompress_grad_it(ublas::vector<int> &time_points,
                            ublas::vector<int> &funorcomp_type,
                            typename ublas::matrix< Type >::iterator1 &data_grad_comp_rit,
                            typename ublas::matrix< Type >::iterator1 &data_grad_uncomp_rit, int nb_t)
    {

      if (time_points.size()!=funorcomp_type.size()) {
        std::cout << "error in basis_function::uncompress_grad, funorcomp_type does not have as much rows as there is time_points_in" << std::endl;
      }
      if (time_points.size()!=data_grad_comp_rit().size1()) {
        std::cout << "error in basis_function::uncompress_grad, data_grad_comp does not have as much rows as there is time_points_in" << std::endl;
      }
      if (nb_grad_comp!=data_grad_comp_rit().size2()) {
        std::cout << "error in basis_function::uncompress_grad, data_grad_comp does not have nb_grad_comp columns" << std::endl;
      }
      if (time_points.size()!=data_grad_uncomp_rit().size1()) {
        std::cout << "error in basis_function::uncompress_grad, data_grad_uncomp does not have as much rows as there is time_points_in" << std::endl;
      }
      if (nb_grad_uncomp!=data_grad_uncomp_rit().size2()) {
        std::cout << "error in basis_function::uncompress_grad, data_grad_uncomp does not have nb_grad_comp columns" << std::endl;
      }


      unsigned int i;
      int k, kk;
      typename ublas::matrix< Type >::iterator2 data_grad_comp_cit;
      typename ublas::matrix< Type >::iterator2 data_grad_uncomp_cit;

      //determination of the intervalle corresponding to a time_point
      ublas::vector<int> intervalles(time_points.size());
      for (i=0; i<time_points.size(); i++) {
        if (time_points(i)==0) {
          intervalles(i) = 0;
        } else {
          intervalles(i) = (time_points(i)*Nint-1)/(nb_t-1);
        }
        //std::cout << Nint << " " << intervalles(i) << std::endl;
      }

      if (!defined_lc) { // the limit conditions are not defined

        if (grad_of_time) {
          for (i=0 ; i<time_points.size() ; i++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<intervalles(i)*nb_fun; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              if (funorcomp_type(i)<0)
                {
                  for (kk=0; kk<nb_grad_comp-1; kk++) {
                    *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                    ++data_grad_comp_cit;
                    ++data_grad_uncomp_cit;
                  }
                }
              else //for a joint angle
                {
                  double fun_scale = funs_scale(funorcomp_type(i));
                  for (k=0; k<4; k++) {
                    for (kk=0; kk<funorcomp_type(i); kk++) {
                      *data_grad_uncomp_cit = 0;
                      ++data_grad_uncomp_cit;
                    }
                    *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                    ++data_grad_uncomp_cit;
                    ++data_grad_comp_cit;
                    for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                      *data_grad_uncomp_cit = 0;
                      ++data_grad_uncomp_cit;
                    }
                  }
                }
              for (kk=intervalles(i)*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
        } else {
          for (i=0 ; i<time_points.size() ; i++)
            {
              data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
              data_grad_comp_cit = data_grad_comp_rit.begin();
              for (kk=0; kk<intervalles(i)*nb_fun; kk++) { //we put first zeros
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              if (funorcomp_type(i)<0)
                {
                  for (kk=0; kk<nb_grad_comp; kk++) {
                    *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                    ++data_grad_comp_cit;
                    ++data_grad_uncomp_cit;
                  }
                }
              else //for a joint angle
                {
                  double fun_scale = funs_scale(funorcomp_type(i));
                  for (k=0; k<4; k++) {
                    for (kk=0; kk<funorcomp_type(i); kk++) {
                      *data_grad_uncomp_cit = 0;
                      ++data_grad_uncomp_cit;
                    }
                    *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                    ++data_grad_uncomp_cit;
                    ++data_grad_comp_cit;
                    for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                      *data_grad_uncomp_cit = 0;
                      ++data_grad_uncomp_cit;
                    }
                  }
                }
              for (kk=intervalles(i)*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) { //we put last zeros
                *data_grad_uncomp_cit = 0;
                ++data_grad_uncomp_cit;
              }
              ++data_grad_comp_rit;
              ++data_grad_uncomp_rit;
            }
        }

        //---------------------------------------------------------------------------
      } else if (defined_lc & !free_fun_lc) { // the limit conditions are completely defined

        if (grad_of_time) { // gradient of time is computed
          for (i=0 ; i<time_points.size() ; i++)
            {
              if (funorcomp_type(i)<0)
                {
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  if (intervalles(i) == 0)
                    {
                      data_grad_comp_cit += 3*nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=nb_fun; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 1)
                    {
                      data_grad_comp_cit += 2*nb_fun;
                      for (kk=0; kk<2*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=2*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 2)
                    {
                      data_grad_comp_cit += nb_fun;
                      for (kk=0; kk<3*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=3*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) < Nint-3)
                    {
                      for (kk=0; kk<(intervalles(i)-3)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_grad_comp-1; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=(intervalles(i)-3)*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-3)
                    {
                      for (kk=0; kk<(Nint-6)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<3*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_comp_cit += nb_fun;
                    }
                  else if (intervalles(i) == Nint-2)
                    {
                      for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<2*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_comp_cit += 2*nb_fun;
                    }
                  else if (intervalles(i) == Nint-1)
                    {
                      for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_comp_cit += 3*nb_fun;
                    }
                  *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
              else //for a joint angle
                {
                  double fun_scale = funs_scale(funorcomp_type(i));
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  if (intervalles(i) == 0)
                    {
                      data_grad_comp_cit += 3;
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                      ++data_grad_uncomp_cit;
                      ++data_grad_comp_cit;
                      for (kk=funorcomp_type(i)+1; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 1)
                    {
                      data_grad_comp_cit += 2;
                      for (k=0; k<2; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      for (kk=2*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 2)
                    {
                      ++data_grad_comp_cit;
                      for (k=0; k<3; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      for (kk=3*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) < Nint-3)
                    {
                      for (kk=0; kk<(intervalles(i)-3)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (k=0; k<4; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      for (kk=(intervalles(i)-3)*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-3)
                    {
                      for (kk=0; kk<(Nint-6)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (k=0; k<3; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      ++data_grad_comp_cit;
                    }
                  else if (intervalles(i) == Nint-2)
                    {
                      for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (k=0; k<2; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      data_grad_comp_cit += 2;
                    }
                  else if (intervalles(i) == Nint-1)
                    {
                      for (kk=0; kk<(Nint-4)*nb_fun+funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                      ++data_grad_uncomp_cit;
                      ++data_grad_comp_cit;
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_comp_cit += 3;
                    }
                  *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
            }
          //---------------------------------------------------------------------------
        } else { // gradient of time is not computed
          for (i=0 ; i<time_points.size() ; i++)
            {
              if (funorcomp_type(i)<0)
                {
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  if (intervalles(i) == 0)
                    {
                      data_grad_comp_cit += 3*nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=nb_fun; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 1)
                    {
                      data_grad_comp_cit += 2*nb_fun;
                      for (kk=0; kk<2*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=2*nb_fun; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 2)
                    {
                      data_grad_comp_cit += nb_fun;
                      for (kk=0; kk<3*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=3*nb_fun; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) < Nint-3)
                    {
                      for (kk=0; kk<(intervalles(i)-3)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_grad_comp-1; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=(intervalles(i)-3)*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-3)
                    {
                      for (kk=0; kk<(Nint-6)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<3*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-2)
                    {
                      for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<2*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-1)
                    {
                      for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
              else //for a joint angle
                {
                  double fun_scale = funs_scale(funorcomp_type(i));
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  if (intervalles(i) == 0)
                    {
                      data_grad_comp_cit += 3;
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                      ++data_grad_uncomp_cit;
                      ++data_grad_comp_cit;
                      for (kk=funorcomp_type(i)+1; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 1)
                    {
                      data_grad_comp_cit += 2;
                      for (k=0; k<2; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      for (kk=2*nb_fun; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 2)
                    {
                      ++data_grad_comp_cit;
                      for (k=0; k<3; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      for (kk=3*nb_fun; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) < Nint-3)
                    {
                      for (kk=0; kk<(intervalles(i)-3)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (k=0; k<4; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      for (kk=(intervalles(i)-3)*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-3)
                    {
                      for (kk=0; kk<(Nint-6)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (k=0; k<3; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                    }
                  else if (intervalles(i) == Nint-2)
                    {
                      for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (k=0; k<2; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                    }
                  else if (intervalles(i) == Nint-1)
                    {
                      for (kk=0; kk<(Nint-4)*nb_fun+funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                      ++data_grad_uncomp_cit;
                      ++data_grad_comp_cit;
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
            }
        } // end gradient of time is not computed

        //---------------------------------------------------------------------------
      } else { // the limit conditions are defined, limit conditions of fun are free
        calc_P_fun_0();
        calc_P_fun_f();

        if (grad_of_time) { // gradient of time is computed
          for (i=0 ; i<time_points.size() ; i++)
            {
              if (funorcomp_type(i)<0)
                {
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  if (intervalles(i) == 0)
                    {
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit*P0_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_uncomp_cit -= nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit += *data_grad_comp_cit*P1_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_uncomp_cit -= nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit += *data_grad_comp_cit*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=2*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 1)
                    {
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit*P1_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_uncomp_cit -= nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit += *data_grad_comp_cit*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<2*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=3*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 2)
                    {
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<3*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=4*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) < Nint-3)
                    {
                      for (kk=0; kk<(intervalles(i)-2)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_grad_comp-1; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=(intervalles(i)-2)*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-3)
                    {
                      for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<3*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-2)
                    {
                      for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<2*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_uncomp_cit -= nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit += *data_grad_comp_cit*Pendm1_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-1)
                    {
                      for (kk=0; kk<(Nint-3)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_uncomp_cit -= nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit += *data_grad_comp_cit*Pendm1_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_uncomp_cit -= nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit += *data_grad_comp_cit*Pend_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
              else //for a joint angle
                {
                  double fun_scale = funs_scale(funorcomp_type(i));
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  if (intervalles(i) == 0)
                    {
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit*P0_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      *data_grad_uncomp_cit += *data_grad_comp_cit*P1_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      *data_grad_uncomp_cit += *data_grad_comp_cit*P2_fun/(funif_scale*fun_scale);
                      ++data_grad_uncomp_cit;
                      ++data_grad_comp_cit;
                      for (kk=0; kk<nb_fun-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                      ++data_grad_comp_cit;
                      ++data_grad_uncomp_cit;
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=2*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 1)
                    {
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit*P1_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      *data_grad_uncomp_cit += *data_grad_comp_cit*P2_fun/(funif_scale*fun_scale);
                      ++data_grad_uncomp_cit;
                      ++data_grad_comp_cit;
                      for (k=0; k<2; k++) {
                        for (kk=0; kk<nb_fun-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=3*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 2)
                    {
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit*P2_fun/(funif_scale*fun_scale);
                      ++data_grad_uncomp_cit;
                      ++data_grad_comp_cit;
                      for (k=0; k<3; k++) {
                        for (kk=0; kk<nb_fun-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=4*nb_fun; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) < Nint-3)
                    {
                      for (kk=0; kk<(intervalles(i)-2)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (k=0; k<4; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      for (kk=(intervalles(i)-2)*nb_fun+nb_grad_comp-1; kk<nb_grad_uncomp-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-3)
                    {
                      for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (k=0; k<3; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      ++data_grad_uncomp_cit;
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-2)
                    {
                      for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (k=0; k<2; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      *data_grad_uncomp_cit += *data_grad_comp_cit*Pendm1_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      ++data_grad_uncomp_cit;
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-1)
                    {
                      for (kk=0; kk<(Nint-3)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                      ++data_grad_uncomp_cit;
                      ++data_grad_comp_cit;
                      for (kk=0; kk<nb_fun-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      *data_grad_uncomp_cit += *data_grad_comp_cit*Pendm1_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      *data_grad_uncomp_cit += *data_grad_comp_cit*Pend_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      ++data_grad_uncomp_cit;
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  *data_grad_uncomp_cit = *data_grad_comp_cit/T_scale;
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
            }
          //---------------------------------------------------------------------------
        } else { // gradient of time is not computed
          for (i=0 ; i<time_points.size() ; i++)
            {
              if (funorcomp_type(i)<0)
                {
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  if (intervalles(i) == 0)
                    {
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit*P0_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_uncomp_cit -= nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit += *data_grad_comp_cit*P1_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_uncomp_cit -= nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit += *data_grad_comp_cit*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=2*nb_fun; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 1)
                    {
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit*P1_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_uncomp_cit -= nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit += *data_grad_comp_cit*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<2*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=3*nb_fun; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 2)
                    {
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit*P2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<3*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=4*nb_fun; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) < Nint-3)
                    {
                      for (kk=0; kk<(intervalles(i)-2)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_grad_comp-1; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=(intervalles(i)-2)*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-3)
                    {
                      for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<3*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-2)
                    {
                      for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<2*nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_uncomp_cit -= nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit += *data_grad_comp_cit*Pendm1_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-1)
                    {
                      for (kk=0; kk<(Nint-3)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_uncomp_cit -= nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit += *data_grad_comp_cit*Pendm1_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      data_grad_uncomp_cit -= nb_fun;
                      for (kk=0; kk<nb_fun; kk++) {
                        *data_grad_uncomp_cit += *data_grad_comp_cit*Pend_fun/(funif_scale*funs_scale(kk%nb_fun));
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
              else //for a joint angle
                {
                  double fun_scale = funs_scale(funorcomp_type(i));
                  data_grad_uncomp_cit = data_grad_uncomp_rit.begin();
                  data_grad_comp_cit = data_grad_comp_rit.begin();
                  if (intervalles(i) == 0)
                    {
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit*P0_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      *data_grad_uncomp_cit += *data_grad_comp_cit*P1_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      *data_grad_uncomp_cit += *data_grad_comp_cit*P2_fun/(funif_scale*fun_scale);
                      ++data_grad_uncomp_cit;
                      ++data_grad_comp_cit;
                      for (kk=0; kk<nb_fun-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                      ++data_grad_comp_cit;
                      ++data_grad_uncomp_cit;
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=2*nb_fun; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 1)
                    {
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit*P1_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      *data_grad_uncomp_cit += *data_grad_comp_cit*P2_fun/(funif_scale*fun_scale);
                      ++data_grad_uncomp_cit;
                      ++data_grad_comp_cit;
                      for (k=0; k<2; k++) {
                        for (kk=0; kk<nb_fun-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=3*nb_fun; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == 2)
                    {
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit*P2_fun/(funif_scale*fun_scale);
                      ++data_grad_uncomp_cit;
                      ++data_grad_comp_cit;
                      for (k=0; k<3; k++) {
                        for (kk=0; kk<nb_fun-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_comp_cit;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=4*nb_fun; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) < Nint-3)
                    {
                      for (kk=0; kk<(intervalles(i)-2)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (k=0; k<4; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      for (kk=(intervalles(i)-2)*nb_fun+nb_grad_comp; kk<nb_grad_uncomp; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-3)
                    {
                      for (kk=0; kk<(Nint-5)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (k=0; k<3; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      ++data_grad_uncomp_cit;
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-2)
                    {
                      for (kk=0; kk<(Nint-4)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (k=0; k<2; k++) {
                        for (kk=0; kk<funorcomp_type(i); kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                        *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                        ++data_grad_uncomp_cit;
                        ++data_grad_comp_cit;
                        for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                          *data_grad_uncomp_cit = 0;
                          ++data_grad_uncomp_cit;
                        }
                      }
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      *data_grad_uncomp_cit += *data_grad_comp_cit*Pendm1_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      ++data_grad_uncomp_cit;
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  else if (intervalles(i) == Nint-1)
                    {
                      for (kk=0; kk<(Nint-3)*nb_fun; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      for (kk=0; kk<funorcomp_type(i); kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit/(P_scale*fun_scale);
                      ++data_grad_uncomp_cit;
                      ++data_grad_comp_cit;
                      for (kk=0; kk<nb_fun-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                      *data_grad_uncomp_cit = *data_grad_comp_cit*Pendm2_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      *data_grad_uncomp_cit += *data_grad_comp_cit*Pendm1_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      *data_grad_uncomp_cit += *data_grad_comp_cit*Pend_fun/(funif_scale*fun_scale);
                      ++data_grad_comp_cit;
                      ++data_grad_uncomp_cit;
                      for (kk=0; kk<nb_fun-funorcomp_type(i)-1; kk++) {
                        *data_grad_uncomp_cit = 0;
                        ++data_grad_uncomp_cit;
                      }
                    }
                  ++data_grad_comp_rit;
                  ++data_grad_uncomp_rit;
                }
            }
        } // end gradient of time is not computed
      } // end limit conditions are defined

    }


    template<class Type>
    void uncompress_grad(ublas::vector<int> &time_points,
                         ublas::vector<int> &funorcomp_type,
                         ublas::matrix< Type > *data_grad_comp,
                         ublas::matrix< Type > *data_grad_uncomp, int nb_t)
    {
      typename ublas::matrix< Type >::iterator1 data_grad_comp_rit=(*data_grad_comp).begin1();
      typename ublas::matrix< Type >::iterator1 data_grad_uncomp_rit=(*data_grad_uncomp).begin1();
      uncompress_grad_it<Type>(time_points, funorcomp_type, data_grad_comp_rit, data_grad_uncomp_rit, nb_t);
    }


  protected:

    void def_Nint(); //define Nint from nb_P

    void solve_P_0(const ublas::vector<double> &fun0);

    void solve_P_f(const ublas::vector<double> &funf);

  private:

    //variables

    double f11v, f12v, f13v;
    double f21v, f22v;
    double f31v;
    double fi1v, fi2v, fi3v, fi4v;

    double df11v, df12v, df13v;
    double df21v, df22v;
    double df31v;
    double dfi1v, dfi2v, dfi3v, dfi4v;

    double d2f11v, d2f12v, d2f13v;
    double d2f21v, d2f22v;
    double d2f31v;
    double d2fi1v, d2fi2v, d2fi3v, d2fi4v;

    double P0_fun, P1_fun, P2_fun; // coefficient between P0, P1, P2 and fun at initial time
    double Pend_fun, Pendm1_fun, Pendm2_fun; // coefficient between Pend, Pend-1, Pend-2 and fun at final time

    //member functions

    //functions for spline3
    //base functions of the first phase
    inline double f11 (double tpow3, double tpow2, double t);
    inline double f12 (double tpow3, double tpow2, double t);
    inline double f13 (double tpow3, double tpow2, double t);

    //base functions of the second phase
    inline double f21 (double tpow3, double tpow2, double t);
    inline double f22 (double tpow3, double tpow2, double t);

    //base functions of the third phase
    inline double f31 (double tpow3, double tpow2, double t);

    //base functions of the phase (3<i<Nint-2)
    inline double fi1 (double tpow3, double tpow2, double t);
    inline double fi2 (double tpow3, double tpow2, double t);
    inline double fi3 (double tpow3, double tpow2, double t);
    inline double fi4 (double tpow3, double tpow2, double t);
    //end of functions for spline3

    //function for dspline3
    //base functions of the first phase
    inline double df11 (double tpow2, double t);
    inline double df12 (double tpow2, double t);
    inline double df13 (double tpow2, double t);

    //base functions of the second phase
    inline double df21 (double tpow2, double t);
    inline double df22 (double tpow2, double t);

    //base functions of the third phase
    inline double df31 (double tpow2, double t);

    //base functions of the phase (3<i<Nint-2)
    inline double dfi1 (double tpow2, double t);
    inline double dfi2 (double tpow2, double t);
    inline double dfi3 (double tpow2, double t);
    inline double dfi4 (double tpow2, double t);
    //end of functions for dspline3

    //function for d2spline3
    //base functions of the first phase
    inline double d2f11 (double t);
    inline double d2f12 (double t);
    inline double d2f13 (double t);

    //base functions of the second phase
    inline double d2f21 (double t);
    inline double d2f22 (double t);

    //base functions of the third phase
    inline double d2f31 (double t);

    //base functions of the phase (3<i<Nint-2)
    inline double d2fi1 (double t);
    inline double d2fi2 (double t);
    inline double d2fi3 (double t);
    inline double d2fi4 (double t);
    //end of functions for d2spline3

    void calc_P_fun_0 (); //compute P0_fun, P1_fun, P2_fun

    void calc_P_fun_f (); //compute Pend_fun, Pendm1_fun, Pendm2_fun
  };
} // end of namespace roboptim

#endif //! ROBOPTIM_TRAJECTORY_BSPLINE_HH
