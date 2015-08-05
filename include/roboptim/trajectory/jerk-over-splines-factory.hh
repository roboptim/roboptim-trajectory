#ifndef ROBOPTIM_JERK_OVER_SPLINES_FACTORY_HH
#define ROBOPTIM_JERK_OVER_SPLINES_FACTORY_HH

#include <roboptim/core/numeric-quadratic-function.hh>
namespace roboptim
{
  namespace trajectory
  {
    /// \brief Factory of the Jerk cost function
    ///
    /// \tparam S Spline type
    /// \tparam T Matrix type
    template <typename S, typename T>
    class JerkOverSplinesFactory
    {
      typedef typename S::value_type value_type;
      typedef typename S::interval_t interval_t;
      typedef typename GenericNumericLinearFunction<T>::vector_t vector_t;
      typedef typename GenericNumericLinearFunction<T>::matrix_t matrix_t;

    public:
      JerkOverSplinesFactory(const std::vector<boost::shared_ptr<S> >& splines, interval_t range);

      /// \brief Updates the time range on which we do the optimization
      void updateRange(interval_t range);

      /// \brief Retrieves the cost function
      boost::shared_ptr<GenericNumericQuadraticFunction<T> > getJerk()
      {
        return f_;
      }

    private:
      /// \brief Shared pointers to the splines of the problem
      const std::vector<boost::shared_ptr<S> >& splines_;

      /// \brief Time range of the oprimization
      interval_t range_;

      /// \brief Inputsize of the cost function
      int inputsize_;

      /// \brief Matrix used for the NumericQuadraticFunction
      matrix_t A_;

      /// \brief Vector used for the NumericQuadraticFunction
      vector_t B_;

      /// \brief Scaling of the Jerk
      value_type scale_;

      /// \brief Shared pointer to the cost function
      boost::shared_ptr<GenericNumericQuadraticFunction<T> > f_;
    };
  }
}

#include <roboptim/trajectory/jerk-over-splines-factory.hxx>

#endif //!ROBOPTIM_JERK_OVER_SPLINES_FACTORY_HH
