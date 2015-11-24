#ifndef ROBOPTIM_TRAJECTORY_JERK_OVER_SPLINES_FACTORY_HH
# define ROBOPTIM_TRAJECTORY_JERK_OVER_SPLINES_FACTORY_HH

# include <vector>
# include <boost/shared_ptr.hpp>

# include <roboptim/core/numeric-quadratic-function.hh>

namespace roboptim
{
  namespace trajectory
  {
    /// \brief Factory generating a jerk cost function over a vector of
    /// splines.
    ///
    /// \tparam S spline type.
    /// \tparam T matrix type.
    template <typename S, typename T>
    class JerkOverSplinesFactory
    {
    public:
      typedef S spline_t;
      typedef boost::shared_ptr<spline_t> splinePtr_t;

      typedef GenericNumericQuadraticFunction<T> numericQuadraticFunction_t;
      typedef boost::shared_ptr<numericQuadraticFunction_t>
        numericQuadraticFunctionPtr_t;
      typedef typename numericQuadraticFunction_t::size_type size_type;
      typedef typename numericQuadraticFunction_t::value_type value_type;
      typedef typename numericQuadraticFunction_t::interval_t interval_t;
      typedef typename numericQuadraticFunction_t::vector_t vector_t;
      typedef typename numericQuadraticFunction_t::matrix_t matrix_t;

    public:
      /// \brief Constructor.
      /// \param splines vector of splines.
      /// \param range time interval on which the jerk is computed.
      /// \param scaling scaling factor.
      JerkOverSplinesFactory (const std::vector<splinePtr_t>& splines,
                              const interval_t& range,
                              value_type scaling = -1.);

      /// \brief Updates the time range on which we do the optimization.
      void updateRange (const interval_t& range);

      /// \brief Retrieves the cost function.
      numericQuadraticFunctionPtr_t getJerk();

    protected:

      /// \brief Try to find a "good" scaling factor.
      value_type chooseScaling () const;

    private:
      /// \brief Shared pointers to the splines of the problem
      const std::vector<splinePtr_t>& splines_;

      /// \brief Time range of the optimization.
      interval_t range_;

      /// \brief Input size of the cost function.
      size_type inputSize_;

      /// \brief Matrix used for the NumericFunction.
      matrix_t A_;

      /// \brief Vector used for the NumericFunction.
      vector_t B_;

      /// \brief Scaling of the Jerk.
      value_type scaling_;

      /// \brief Shared pointer to the cost function.
      numericQuadraticFunctionPtr_t f_;
    };
  }
}

# include <roboptim/trajectory/jerk-over-splines-factory.hxx>

#endif //! ROBOPTIM_TRAJECTORY_JERK_OVER_SPLINES_FACTORY_HH
