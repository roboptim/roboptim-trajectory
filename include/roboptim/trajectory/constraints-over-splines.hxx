#ifndef ROBOPTIM_TRAJECTORY_CONSTRAINTS_OVER_SPLINES_HXX
# define ROBOPTIM_TRAJECTORY_CONSTRAINTS_OVER_SPLINES_HXX

# include <boost/format.hpp>

namespace roboptim
{
  namespace trajectory
  {
    template <typename T, typename S>
    ConstraintsOverSplines<T, S>::ConstraintsOverSplines
    (const splines_t& splines, size_t splineIdx, unsigned int order,
     value_type startingPoint, size_type inputSize)
    : GenericDifferentiableFunction<T>
        (inputSize, 2, GenerateName (splines, splineIdx,
                                     order, startingPoint)),
      dimension_ (static_cast<size_t> (splines[splineIdx]->dimension ())),
      interval_ (ComputeInterval (splines, splineIdx, startingPoint)),
      startingIndex_ (ComputeStartIdx (splines, splineIdx)),
      p_ (),
      coefs_ ()
    {
      size_t interval_idx = ComputeIntervalIdx (splines, splineIdx, startingPoint);

      for (size_t j = 0; j < dimension_ + 1; ++j)
      {
        polynomial_t p;
        switch (order)
        {
          case 1 :
          p = splines[splineIdx]->basisPolynomials()[interval_idx + j][dimension_ - j].template derivative<1>();
            break;
          case 2 :
          p = splines[splineIdx]->basisPolynomials()[interval_idx + j][dimension_ - j].template derivative<2>();
            break;
          default :
          p = splines[splineIdx]->basisPolynomials()[interval_idx + j][dimension_ - j];
            break;
        }
        coefs_.push_back(p);
      }
    }

    template <typename T, typename S>
    std::string ConstraintsOverSplines<T, S>::GenerateName
    (const splines_t& splines, size_t splineIdx, unsigned int order,
     value_type startingPoint)
    {
      std::string name;
      interval_t range = ComputeInterval (splines, splineIdx, startingPoint);

      if (order == 0)
        name = (boost::format ("f(t) / t ∈ [%.2f,%.2f] (%s)")
            % range.first % range.second
            % splines[splineIdx]->getName ()).str ();
      else
        name = (boost::format ("f^(%i)(t) / t ∈ [%.2f,%.2f] (%s)")
            % order % range.first % range.second
            % splines[splineIdx]->getName ()).str ();

      return name;
    }

    template <typename T, typename S>
    typename ConstraintsOverSplines<T, S>::size_type
    ConstraintsOverSplines<T, S>::ComputeStartIdx (const splines_t& splines,
                                                   size_t splineIdx)
    {
      size_type start_idx = 0;
      for (size_t i = 0; i < splineIdx; ++i)
        start_idx += static_cast<size_type>
                       (splines[i]->getNumberControlPoints ());
      return start_idx;
    }

    template <typename T, typename S>
    size_t ConstraintsOverSplines<T, S>::ComputeIntervalIdx
    (const splines_t& splines, size_t splineIdx, value_type startingPoint)
    {
      size_t dimension = static_cast<size_t> (splines[splineIdx]->dimension ());
      size_t interval = static_cast<size_t>
        (splines[splineIdx]->interval (startingPoint));

      return interval - dimension;
    }

    template <typename T, typename S>
    typename ConstraintsOverSplines<T, S>::interval_t
    ConstraintsOverSplines<T, S>::ComputeInterval
    (const splines_t& splines, size_t splineIdx, value_type startingPoint)
    {
      size_t dimension = static_cast<size_t> (splines[splineIdx]->dimension ());
      size_t interval_idx = ComputeIntervalIdx (splines, splineIdx,
                                                startingPoint);

      size_type idx = static_cast<size_type> (interval_idx + dimension + 1);

      return S::makeInterval (startingPoint,
                              splines[splineIdx]->knotVector ()[idx]);
    }

    template <typename T, typename S>
    void ConstraintsOverSplines<T, S>::update (const_argument_ref x) const
    {
      p_.coefs().setZero();
      for (size_t i = 0; i < dimension_ + 1; ++i)
        p_ += coefs_[i] * x[static_cast<size_type> (i)];
    }

    template <typename T, typename S>
    void ConstraintsOverSplines<T, S>::impl_compute (result_ref result, const_argument_ref x) const
    {
      this->update (x.middleRows (static_cast<long> (startingIndex_),
                                  static_cast<long> (dimension_ + 1)));
      result[0] = p_.min(interval_).second;
      result[1] = p_.max(interval_).second;
    }

    template <typename T, typename S>
    void ConstraintsOverSplines<T, S>::impl_gradient
    (gradient_ref grad, const_argument_ref x, size_type functionId) const
    {
      update (x.middleRows (startingIndex_,
                            static_cast<size_type> (dimension_ + 1)));

      if (functionId == 0)
      {
        value_type tmin_ = p_.min(interval_).first;
        for (size_t i = 0; i < dimension_ + 1; ++i)
          grad.coeffRef (startingIndex_ + static_cast<size_type> (i))
            = coefs_[i](tmin_);
      }
      else
      {
        value_type tmax_ = p_.max(interval_).first;
        for (size_t i = 0; i < dimension_ + 1; ++i)
          grad.coeffRef (startingIndex_ + static_cast<size_type> (i))
            = coefs_[i](tmax_);
      }
    }
  }
}

#endif //! ROBOPTIM_TRAJECTORY_CONSTRAINTS_OVER_SPLINES_HXX
