#ifndef ROBOPTIM_TRAJECTORY_CONSTRAINTS_OVER_SPLINES_HH
# define ROBOPTIM_TRAJECTORY_CONSTRAINTS_OVER_SPLINES_HH

# include <vector>

# include <boost/shared_ptr.hpp>

# include <roboptim/core/differentiable-function.hh>

namespace roboptim
{
  namespace trajectory
  {
    /// \brief Constraint function on a spline's interval
    ///
    /// \tparam T Matrix type
    /// \tparam S Spline type
    template <typename T, typename S>
    class ConstraintsOverSplines
    : public roboptim::GenericDifferentiableFunction<T>
    {
    public:
      ROBOPTIM_DIFFERENTIABLE_FUNCTION_FWD_TYPEDEFS_
        (roboptim::GenericDifferentiableFunction<T>);

      typedef typename S::polynomial_t polynomial_t;
      typedef Function::interval_t interval_t;
      typedef std::vector<polynomial_t> polynomials_t;
      typedef S spline_t;
      typedef boost::shared_ptr<spline_t> splinePtr_t;
      typedef std::vector<boost::shared_ptr<S> > splines_t;

      /// \brief Constructor.
      ///
      /// \param splines vector of splines.
      /// \param splineIdx index of the spline to constrain.
      /// \param order derivation order of the spline to constain.
      /// \param startingPoint starting point of the constraint.
      /// \param inputSize input size of the problem.
      ConstraintsOverSplines (const splines_t& splines, size_t splineIdx,
                              unsigned int order, value_type startingPoint,
                              size_type inputSize);

    protected:

      /// \brief Generate the name of the constraint.
      /// \param splines vector of splines.
      /// \param splineIdx index of the spline to constrain.
      /// \param order derivation order of the spline to constain.
      /// \param startingPoint starting point of the constraint.
      static std::string GenerateName (const splines_t& splines,
                                       size_t splineIdx,
                                       unsigned int order,
                                       value_type startingPoint);

      /// \brief Compute the range of the time interval.
      /// \param splines vector of splines.
      /// \param splineIdx index of the spline to constrain.
      /// \param startingPoint starting point of the constraint.
      static interval_t ComputeInterval (const splines_t& splines,
                                         size_t splineIdx,
                                         value_type startingPoint);

      /// \brief Determine the interval index of the constraint.
      /// \param splines vector of splines.
      /// \param splineIdx index of the spline to constrain.
      /// \param startingPoint starting point of the constraint.
      static size_t ComputeIntervalIdx (const splines_t& splines,
                                        size_t splineIdx,
                                        value_type startingPoint);

      /// \brief Determine the first index in the argument vector, assuming all
      /// joint control points are contiguous.
      /// \param splines vector of splines.
      /// \param splineIdx index of the spline to constrain.
      static size_type ComputeStartIdx (const splines_t& splines,
                                        size_t splineIdx);

      polynomial_t toPoly (const_argument_ref x) const;

      virtual void impl_compute (result_ref result,
                                 const_argument_ref x) const;

      virtual void impl_gradient (gradient_ref grad,
                                  const_argument_ref x,
                                  size_type i) const;

      /// \brief Order of the spline (e.g. 3 for cubic).
      size_t order_;

      /// \brief Time range of the constraint.
      typename S::interval_t interval_;

      /// \brief Index of the influencing control points.
      size_type startingIndex_;

      /// \brief Basis polynomials used for gradient computation.
      polynomials_t basisPolynomials_;

      /// \brief Constrained spline.
      mutable spline_t spline_;

      /// \brief Index of the constrained time interval.
      size_t intervalIdx_;
    };
  }
}

#include <roboptim/trajectory/constraints-over-splines.hxx>

#endif //! ROBOPTIM_TRAJECTORY_CONSTRAINTS_OVER_SPLINES_HH
