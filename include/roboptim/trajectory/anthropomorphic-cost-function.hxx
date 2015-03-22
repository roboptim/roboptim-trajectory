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

#ifndef ROBOPTIM_TRAJECTORY_ANTHROPOMORPHIC_COST_FUNCTION_HXX
# define ROBOPTIM_TRAJECTORY_ANTHROPOMORPHIC_COST_FUNCTION_HXX
# include <roboptim/trajectory/frontal-speed.hh>
# include <roboptim/trajectory/orthogonal-speed.hh>

namespace roboptim
{
namespace trajectory
{
  template <typename T>
  AnthropomorphicCostFunction<T>::AnthropomorphicCostFunction
  (const T& trajectory,
   const typename boost::optional<vector_t> alpha,
   const typename boost::optional<value_type> ksi1,
   const typename boost::optional<value_type> ksi2)
    : DerivableFunction (trajectory.parameters ().size (), 1,
			 "cost function"),
      trajectory_ (trajectory),
      alpha_ (alpha ? *alpha : defaultAlpha ()),
      ksi1_ (ksi1 ? *ksi1 : defaultKsi1 ()),
      ksi2_ (ksi2 ? *ksi2 : defaultKsi2 ()),
      alpha3_ ()
  {
    assert (alpha_.size () == 4 && "Invalid alpha array.");

    vector_t p = trajectory.parameters ();

    //FIXME: +1 as we're doing free time trajectory.
    //should be generic.
    const value_type dx = p[p.size () - 3] - p[1 + 0];
    const value_type dy = p[p.size () - 2] - p[1 + 1];
    const value_type dtheta = p[p.size () - 1] - p[1 + 2];

    alpha3_ = this->alpha3 (std::fabs (dtheta), dx * dx + dy * dy);
  }

  template <typename T>
  AnthropomorphicCostFunction<T>::~AnthropomorphicCostFunction ()
  {}

  namespace detail
  {
    template <typename T>
    struct ComputeIntegral
    {
      ComputeIntegral (const T& traj,
		       Function::const_vector_ref alpha,
		       const double& alpha3,
		       double& res)
	: traj_ (traj),
	  alpha_ (alpha),
	  alpha3_ (alpha3),
	  res_ (res)
      {}

      void operator () (const double& t)
      {
	static Function::vector_t t_ (1);
	FrontalSpeed frontalSpeed;
	OrthogonalSpeed orthogonalSpeed;
	t_[0] = t;
	const Function::value_type u1 = frontalSpeed.gradient (traj_.state (t, 1))[0];
	const Function::value_type u2 = traj_.derivative (t, 2)[2];
	const Function::value_type u3 = orthogonalSpeed.gradient (traj_.state (t, 1))[0];
	res_ +=
	  alpha_[0]
	  + alpha_[1] * u1 * u1
	  + alpha_[2] * u2 * u2
	  + alpha3_ * u3 * u3;
      }

    private:
      const T& traj_;
      Function::const_vector_ref alpha_;
      const double& alpha3_;
      double& res_;
    };
  } // end of namespace detail.

  template <typename T>
  void
  AnthropomorphicCostFunction<T>::impl_compute (result_ref res,
						const_argument_ref p)
    const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
      Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    res.setZero();

    static T updatedTrajectory = trajectory_;
    updatedTrajectory.setParameters (p);

    const size_type nbDiscretizationPoints = 50;
    detail::ComputeIntegral<T> ci (updatedTrajectory, alpha_, alpha3_, res[0]);
    foreach (updatedTrajectory.timeRange (), nbDiscretizationPoints, ci);

    res[0] *= updatedTrajectory.length () / nbDiscretizationPoints;
  }

  template <typename T>
  void
  AnthropomorphicCostFunction<T>::impl_gradient (gradient_ref grad,
						 const_argument_ref p,
						 size_type)
    const
  {
#ifndef ROBOPTIM_DO_NOT_CHECK_ALLOCATION
      Eigen::internal::set_is_malloc_allowed (true);
#endif //! ROBOPTIM_DO_NOT_CHECK_ALLOCATION

    GenericFiniteDifferenceGradient<EigenMatrixDense> fdfunction (*this);
    fdfunction.gradient (grad, p, 0);
  }

  template <typename T>
  typename AnthropomorphicCostFunction<T>::value_type
  AnthropomorphicCostFunction<T>::alpha3 (value_type deltaTheta,
					  value_type dsquare) const
  {
    const value_type ksi1 = M_PI / 18.;
    const value_type ksi2 = .5;

    return alpha_[3] * (1. + (deltaTheta / ksi1)) * (1 + (dsquare / ksi2));
  }

} // end of namespace trajectory.
} // end of namespace roboptim.

#endif //! ROBOPTIM_TRAJECTORY_ANTHROPOMORPHIC_COST_FUNCTION_HXX
