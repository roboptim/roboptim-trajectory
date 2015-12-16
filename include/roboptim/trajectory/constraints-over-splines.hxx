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
      order_ (static_cast<size_t> (splines[splineIdx]->order ())),
      interval_ (ComputeInterval (splines, splineIdx, startingPoint)),
      startingIndex_ (ComputeStartIdx (splines, splineIdx)),
      basisPolynomials_ (),
      spline_ (*splines[splineIdx]),
      intervalIdx_ (ComputeIntervalIdx (splines, splineIdx, startingPoint))
    {
      for (size_t j = 0; j < order_ + 1; ++j)
      {
        polynomial_t p;
        switch (order)
        {
          case 1:
          {
            p = splines[splineIdx]->basisPolynomials ()
                [intervalIdx_ + j][order_ - j].template derivative<1>();
            break;
          }

          case 2:
          {
            p = splines[splineIdx]->basisPolynomials ()
                [intervalIdx_ + j][order_ - j].template derivative<2>();
            break;
          }

          default:
          case 0:
          {
            p = splines[splineIdx]->basisPolynomials ()
                [intervalIdx_ + j][order_ - j];
            break;
          }
        }

        basisPolynomials_.push_back (p);
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
        name = (boost::format ("f(t) / t ∈ [%g,%g] (%s)")
            % range.first % range.second
            % splines[splineIdx]->getName ()).str ();
      else
        name = (boost::format ("f^(%i)(t) / t ∈ [%g,%g] (%s)")
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
      size_t order = static_cast<size_t> (splines[splineIdx]->order ());
      size_t interval = static_cast<size_t>
        (splines[splineIdx]->interval (startingPoint));

      return interval - order;
    }

    template <typename T, typename S>
    typename ConstraintsOverSplines<T, S>::interval_t
    ConstraintsOverSplines<T, S>::ComputeInterval
    (const splines_t& splines, size_t splineIdx, value_type startingPoint)
    {
      size_t order = static_cast<size_t> (splines[splineIdx]->order ());
      size_t interval_idx = ComputeIntervalIdx (splines, splineIdx,
                                                startingPoint);

      size_type idx = static_cast<size_type> (interval_idx + order + 1);

      return S::makeInterval (startingPoint,
                              splines[splineIdx]->knotVector ()[idx]);
    }

    template <typename T, typename S>
    void ConstraintsOverSplines<T, S>::impl_compute
    (result_ref result, const_argument_ref x) const
    {
      polynomial_t p = toPoly (x);

      result[0] = p.min (interval_).second;
      result[1] = p.max (interval_).second;
    }

    template <typename T, typename S>
    typename ConstraintsOverSplines<T, S>::polynomial_t
    ConstraintsOverSplines<T, S>::toPoly (const_argument_ref x) const
    {
      typename spline_t::const_argument_ref
        param = x.segment (static_cast<size_type> (startingIndex_),
                           spline_.parameters ().size ());
      const size_type n = static_cast<size_type> (spline_.outputSize ());
      const size_type offset = n - 1;
      const size_type k = static_cast<size_type> (intervalIdx_ + order_);

      polynomial_t p;
      for (size_type i = 0; i <= static_cast<size_type> (order_); ++i)
      {
        p += basisPolynomials_[order_ - i] * param (n * (k - i) + offset);
      }

      return p;
    }

    template <typename T, typename S>
    void ConstraintsOverSplines<T, S>::impl_gradient
    (gradient_ref grad, const_argument_ref x, size_type functionId) const
    {
      polynomial_t p = toPoly (x);

      if (functionId == 0)
      {
        value_type tmin_ = p.min (interval_).first;
        for (size_t i = 0; i < order_ + 1; ++i)
          grad.coeffRef (startingIndex_
                         + static_cast<size_type> (intervalIdx_)
                         + static_cast<size_type> (i))
            = basisPolynomials_[i](tmin_);
      }
      else
      {
        value_type tmax_ = p.max (interval_).first;
        for (size_t i = 0; i < order_ + 1; ++i)
          grad.coeffRef (startingIndex_
                         + static_cast<size_type> (intervalIdx_)
                         + static_cast<size_type> (i))
            = basisPolynomials_[i](tmax_);
      }
    }
  }
}

#endif //! ROBOPTIM_TRAJECTORY_CONSTRAINTS_OVER_SPLINES_HXX
