#ifndef ROBOPTIM_JERK_OVER_SPLINES_FACTORY_HXX
#define ROBOPTIM_JERK_OVER_SPLINES_FACTORY_HXX

namespace roboptim
{
  namespace trajectory
  {
    template <typename S, typename T>
    JerkOverSplinesFactory<S, T>::JerkOverSplinesFactory(const std::vector<boost::shared_ptr<S> >& splines, interval_t range)
      : splines_ (splines),
	range_ (range),
	inputsize_(static_cast<int>(splines.size()*(static_cast<unsigned long>(splines[0]->getNumberControlPoints())))),
	A_ (inputsize_, inputsize_),
	B_ (inputsize_),
  scale_ (1e-6),
	f_ ()
    {
      assert (splines[0]->dimension() == 3);
      ROBOPTIM_DEBUG_ONLY(
                          long nbp_ = splines[0]->getNumberControlPoints();
                          for (unsigned long i = 1; i < splines.size(); ++i)
			    assert (nbp_ == splines[i]->getNumberControlPoints());
			  );
      A_.setZero();
      B_.setZero();
      int start = static_cast<int>(splines[0]->interval(range.first));
      int end = static_cast<int>(splines[0]->interval(range.second));
      assert (range.first < range.second);
      int N = splines[0]->dimension();
      int nbParameters = static_cast<int>(splines[0]->getNumberControlPoints());
      for (unsigned long n = 0; n < splines.size(); ++n)
        for (int i = start; i < end + 1; ++i)
          for (int j = N; j >= 0; --j)
            for (int k = N; k >= 0; --k)
              A_.coeffRef(i+static_cast<int>(n)*nbParameters-j, i+static_cast<int>(n)*nbParameters-k) +=
                splines_[n]->basisPolynomials()[static_cast<unsigned long>(i-j)][static_cast<unsigned long>(j)].derivative(0,3) *
                splines_[n]->basisPolynomials()[static_cast<unsigned long>(i-k)][static_cast<unsigned long>(k)].derivative(0,3) *
                (std::min(splines_[n]->knotVector()[i+1], range.second) - std::max(splines_[n]->knotVector()[i], range.first));
      A_ *= scale_;
      f_ = boost::make_shared<GenericNumericQuadraticFunction<T> >(A_, B_);
    }

    template <typename S, typename T>
    void JerkOverSplinesFactory<S, T>::updateRange(interval_t range)
    {
      assert (range.first < range.second);
      int newStart = static_cast<int>(splines_[0]->interval(range.first));
      int newEnd = static_cast<int>(splines_[0]->interval(range.second));
      int N = splines_[0]->dimension();
      int nbParameters = static_cast<int>(splines_[0]->getNumberControlPoints());
      A_.setZero();
      for (unsigned long n = 0; n < splines_.size(); ++n)
        for (int i = newStart; i < newEnd + 1; ++i)
          for (int j = N; j >= 0; --j)
            for (int k = N; k >= 0; --k)
              A_.coeffRef(i+static_cast<int>(n)*nbParameters-j, i+static_cast<int>(n)*nbParameters-k) +=
                splines_[n]->basisPolynomials()[static_cast<unsigned long>(i-j)][static_cast<unsigned long>(j)].derivative(0,3) *
                splines_[n]->basisPolynomials()[static_cast<unsigned long>(i-k)][static_cast<unsigned long>(k)].derivative(0,3) *
                (std::min(splines_[n]->knotVector()[i+1], range.second) - std::max(splines_[n]->knotVector()[i], range.first));
      range_ = range;
      A_ *= scale_;
      f_ = boost::make_shared<GenericNumericQuadraticFunction<T> >(A_, B_);
    }
  }
}

#endif //!ROBOPTIM_JERK_OVER_SPLINES_FACTORY_HXX
