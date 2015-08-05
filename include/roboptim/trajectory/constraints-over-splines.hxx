#ifndef ROBOPTIM_TRAJECTORY_CONSTRAINTS_OVER_SPLINES_HXX
#define ROBOPTIM_TRAJECTORY_CONSTRAINTS_OVER_SPLINES_HXX

namespace roboptim
{
  namespace trajectory
  {
    template <typename T, typename S>
    ConstraintsOverSplines<T, S>::ConstraintsOverSplines(const splines_t& splines, unsigned long splineNumber, unsigned int order, value_type startingPoint, int inputSize)
    : GenericDifferentiableFunction<T>(inputSize, 2),
    dimension_ (static_cast<unsigned long>(splines[splineNumber]->dimension())),
    startingIndex_ (0),
    coefs_ ()
    {
      for (unsigned long i = 0; i < splineNumber; ++i)
        startingIndex_ += static_cast<unsigned long>(splines[i]->getNumberControlPoints());
      unsigned long intervalNumber_ = static_cast<unsigned long>(splines[splineNumber]->interval(startingPoint)) - dimension_;
      startingIndex_ += intervalNumber_;
      interval_ = this->makeInterval(startingPoint, splines[splineNumber]->knotVector()[static_cast<long>(intervalNumber_ + dimension_ + 1)]);
      for (unsigned long j = 0; j < dimension_ + 1; ++j)
      {
        polynomial_t p;
        switch (order)
        {
          case 1 :
          p = splines[splineNumber]->basisPolynomials()[intervalNumber_ + j][dimension_ - j].template derivative<1>();
            break;
          case 2 :
          p = splines[splineNumber]->basisPolynomials()[intervalNumber_ + j][dimension_ - j].template derivative<2>();
            break;
          default :
          p = splines[splineNumber]->basisPolynomials()[intervalNumber_ + j][dimension_ - j];
            break;
        }
        coefs_.push_back(p);
      }
    }

    template <typename T, typename S>
    void ConstraintsOverSplines<T, S>::update (const_argument_ref x) const
    {
      p_.coefs().setZero();
      for (unsigned long i = 0; i < dimension_ + 1; ++i)
        p_ += coefs_[i] * x[static_cast<long>(i)];
    }

    template <typename T, typename S>
    void ConstraintsOverSplines<T, S>::impl_compute (result_ref result, const_argument_ref x) const
    {
      this->update(x.middleRows(static_cast<long>(startingIndex_), static_cast<long>(dimension_ + 1)));
      result[0] = p_.min(interval_).second;
      result[1] = p_.max(interval_).second;
    }

    template <typename T, typename S>
    void ConstraintsOverSplines<T, S>::impl_gradient (gradient_ref grad, const_argument_ref x, size_type functionId) const
    {
      this->update(x.middleRows(static_cast<long>(startingIndex_), static_cast<long>(dimension_ + 1)));
      if (functionId == 0)
      {
        value_type tmin_ = p_.min(interval_).first;
        for (unsigned long i = 0; i < dimension_ + 1; ++i)
          grad.coeffRef(static_cast<int>(startingIndex_+i)) = coefs_[i](tmin_);
      }
      else
      {
        value_type tmax_ = p_.max(interval_).first;
        for (unsigned long i = 0; i < dimension_ + 1; ++i)
          grad.coeffRef(static_cast<int>(startingIndex_+i)) = coefs_[i](tmax_);
      }
    }
  }
}

#endif //! ROBOPTIM_TRAJECTORY_CONSTRAINTS_OVER_SPLINES_HXX
