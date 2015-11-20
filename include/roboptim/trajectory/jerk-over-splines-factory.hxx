#ifndef ROBOPTIM_TRAJECTORY_JERK_OVER_SPLINES_FACTORY_HXX
# define ROBOPTIM_TRAJECTORY_JERK_OVER_SPLINES_FACTORY_HXX

namespace roboptim
{
  namespace trajectory
  {
    template <typename S, typename T>
    JerkOverSplinesFactory<S, T>::JerkOverSplinesFactory
    (const std::vector<splinePtr_t>& splines, const interval_t& range)
      : splines_ (splines),
        range_ (range),
        inputSize_ (static_cast<size_type> (splines.size ())
                    * static_cast<size_type> (splines[0]->getNumberControlPoints ())),
        A_ (inputSize_, inputSize_),
        B_ (inputSize_),
        scale_ (1e-6),
        f_ ()
    {
      assert (splines[0]->dimension () == 3);
      ROBOPTIM_DEBUG_ONLY
        (typename spline_t::size_type nbp_ = splines[0]->getNumberControlPoints();
         for (size_t i = 1; i < splines.size (); ++i)
           assert (nbp_ == splines[i]->getNumberControlPoints());
        );

      A_.setZero ();
      B_.setZero ();

      assert (range.first < range.second);
      size_t start = static_cast<size_t> (splines[0]->interval (range.first));
      size_t end = static_cast<size_t> (splines[0]->interval (range.second));

      size_t N = static_cast<size_t> (splines[0]->dimension ());
      size_t nbParameters = static_cast<size_t> (splines[0]->getNumberControlPoints());

      for (size_t n = 0; n < splines.size (); ++n)
        for (size_t i = start; i < end + 1; ++i)
          for (size_t j = 0; j < N+1; ++j)
            for (size_t k = 0; k < N+1; ++k)
            {
              size_type ii = static_cast<size_type> (i);
              A_.coeffRef (static_cast<size_type> (i + n*nbParameters - j),
                           static_cast<size_type> (i + n*nbParameters - k))
                += splines_[n]->basisPolynomials ()[i-j][j].derivative (0,3) *
                   splines_[n]->basisPolynomials ()[i-k][k].derivative (0,3) *
                   (std::min (splines_[n]->knotVector ()[ii+1], range.second)
                    - std::max (splines_[n]->knotVector ()[ii], range.first));
            }

      A_ *= scale_;
      f_ = boost::make_shared<numericQuadraticFunction_t> (A_, B_);
    }

    template <typename S, typename T>
    void JerkOverSplinesFactory<S, T>::updateRange (const interval_t& range)
    {
      assert (range.first < range.second);

      size_t newStart = static_cast<size_t> (splines_[0]->interval (range.first));
      size_t newEnd = static_cast<size_t> (splines_[0]->interval (range.second));
      size_t N = static_cast<size_t> (splines_[0]->dimension ());
      size_t nbParameters = static_cast<size_t> (splines_[0]->getNumberControlPoints ());
      A_.setZero ();

      for (size_t n = 0; n < splines_.size (); ++n)
        for (size_t i = newStart; i < newEnd + 1; ++i)
          for (size_t j = 0; j < N+1; ++j)
            for (size_t k = 0; k < N+1; ++k)
            {
              size_type ii = static_cast<size_type> (i);
              A_.coeffRef (static_cast<size_type> (i + n*nbParameters - j),
                           static_cast<size_type> (i + n*nbParameters - k))
                += splines_[n]->basisPolynomials ()[i-j][j].derivative (0,3) *
                   splines_[n]->basisPolynomials ()[i-k][k].derivative (0,3) *
                   (std::min (splines_[n]->knotVector ()[ii+1], range.second)
                    - std::max (splines_[n]->knotVector ()[ii], range.first));
            }

      range_ = range;
      A_ *= scale_;
      f_ = boost::make_shared<numericQuadraticFunction_t> (A_, B_);
    }

    template <typename S, typename T>
    typename JerkOverSplinesFactory<S, T>::numericQuadraticFunctionPtr_t
    JerkOverSplinesFactory<S, T>::getJerk ()
    {
      return f_;
    }
  }
}

#endif //! ROBOPTIM_TRAJECTORY_JERK_OVER_SPLINES_FACTORY_HXX
