// Copyright (C) 2009, 2010 by Thomas Moulard, AIST, CNRS, INRIA.
// Copyright (C) 2013 by Alexander Werner, DLR.
// Copyright (C) 2014 by Benjamin Chrétien, CNRS-LIRMM.
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

#ifndef ROBOPTIM_TRAJECTORY_B_SPLINE_HH
# define ROBOPTIM_TRAJECTORY_B_SPLINE_HH
# include <roboptim/trajectory/sys.hh>
# include <roboptim/trajectory/deprecated.hh>

# include <roboptim/trajectory/trajectory.hh>
# include <roboptim/trajectory/polynomial.hh>
# include <roboptim/trajectory/fwd.hh>

# include <map>


namespace roboptim
{
namespace trajectory
{
  /// \addtogroup roboptim_function
  /// @{
  /// B-Spline trajectory.
  /// Implement a B-Spline as a trajectory based in doc/quintic-b-spline.tex
  /// The template parameter N determines the order of the spline.
  template <int N>
  class BSpline : public Trajectory<N>
  {
  public:
    /// \brief Parent type and imports.
    ROBOPTIM_NTIMES_DERIVABLE_FUNCTION_FWD_TYPEDEFS_
      (Trajectory<N>);

    typedef typename parent_t::interval_t interval_t;

    typedef Polynomial<N> polynomial_t;
    typedef Monomial<N> monomial_t;
    typedef std::map<
      int,
      polynomial_t,
      std::less<int>,
      Eigen::aligned_allocator<std::pair<const int, polynomial_t> >
      > cox_map;
    typedef typename cox_map::iterator cox_map_itr_t;

    typedef std::vector <
      polynomial_t,
      Eigen::aligned_allocator<polynomial_t>
      > basisPolynomials_t;
    typedef std::vector <basisPolynomials_t>
    basisPolynomialsVector_t;

    /// \brief Instantiate a B-Spline from its definition.
    ///
    /// \param timeRange spline time range: $\f$[t_3,t_n]\f$
    /// \param dimension spline dimension: \f$n\f$
    /// \param parameters vector of parameters defining control points
    /// \param name function title
    /// \param clamped whether the spline should be clamped
    ///
    /// Number of control points is inferred from dimension of parameter
    /// vector.
    BSpline (const interval_t& timeRange, size_type dimension,
             const vector_t& parameters,
             const std::string name = "B-Spline",
             bool clamped = false);

    /// \brief Instantiate a B-Spline with parameters and knot points.
    ///
    /// \param timeRange spline time range: $\f$[t_3,t_n]\f$
    /// \param dimension spline dimension: \f$n\f$
    /// \param parameters vector of parameters
    /// \param knots vector of control points
    /// \param name function title
    /// The number of knot points must be the number of
    /// parameters + N + 1.
    /// In the knot vector, N knots at the beginning
    /// must lie before the start of the spline.
    /// The rest of the knot point must lie before the
    /// end of the spline interval.
    BSpline (const interval_t& tr, size_type dimension,
             const vector_t& parameters,
             const_vector_ref knots,
             std::string name = "B-Spline");

    /// \brief Copy constructor.
    /// \param spline spline that will be copied.
    BSpline (const BSpline<N>& spline);

    /// \brief Compute the derivative of the B-spline for a given order.
    ///
    /// Note: here we return a spline of the same degree for simplicity (we
    /// keep the same knot vector), but the basis polynomials are M degrees
    /// lower.
    ///
    /// \tparam M derivation order (M ≥ 0).
    ///
    /// \return derived spline.
    template <int M>
    BSpline<N> derivative () const;

    /// \brief Virtual destructor.
    virtual ~BSpline () {};

    /// \brief Modify spline parameters.
    virtual void setParameters (const vector_t&);

    virtual jacobian_t variationConfigWrtParam (value_type t) const;
    virtual jacobian_t variationDerivWrtParam (value_type t, size_type order)
      const;
    virtual value_type singularPointAtRank (size_type rank) const;
    virtual vector_t derivBeforeSingularPoint (size_type rank, size_type order)
      const;
    virtual vector_t derivAfterSingularPoint (size_type rank, size_type order)
      const;

    ROBOPTIM_IMPLEMENT_CLONE (BSpline<N>)

    virtual Trajectory<N>* resize (interval_t timeRange) const;

    /// \brief Display the function on the specified output stream.
    ///
    /// \param o output stream used for display
    /// \return output stream
    virtual std::ostream& print (std::ostream& o) const;

    jacobian_t variationConfigWrtParam (StableTimePoint tp) const;
    jacobian_t variationDerivWrtParam (StableTimePoint tp, size_type order)
      const;

    value_type Dt () const ROBOPTIM_TRAJECTORY_DEPRECATED;

    /// \brief Return the knot vector of the spline.
    /// \return knot vector of the spline (const).
    const vector_t& knotVector () const;

    size_type interval (value_type t) const;

    /// \brief Get the number of control points of the spline.
    /// \return Number of control points of the spline.
    size_type getNumberControlPoints() const
    {
      return nbp_;
    }

    /// \brief Return the polynomial expression of the B-spline on each
    /// time interval.
    void toPolynomials (basisPolynomials_t& res) const;

    /// \brief Constant getter for the basis polynomials of the B-spline.
    ///
    /// Note: computeBasisPolynomials() needs to be called beforehand (which is
    /// done in the BSpline constructor).
    ///
    /// \return constant reference to the basis polynomials.
    const basisPolynomialsVector_t&
    basisPolynomials () const
    {
      return basisPolynomials_;
    }

    using Trajectory<N>::derivative;

    /// \brief Retrieve the dimension of the BSpline
    int dimension () const;

  protected:

    /// \brief Enum for special constructors.
    enum ConstructionMode
      {
        /// \brief Automatically generate basis functions and knot vector.
        NORMAL = 0,
        /// \brief Do not generate basis functions or knot vector
        /// automatically.
        UNINITIALIZED = 1
      };

    /// \brief Special constructor used internally for specific
    /// cases (e.g. spline derivation).
    ///
    /// \param timeRange time range of the spline.
    /// \param dimension output dimension.
    /// \param parameters control point vector.
    /// \param name name of the spline.
    /// \param clamped whether the spline is clamped.
    BSpline (ConstructionMode mode,
             const interval_t& timeRange, size_type dimension,
             const vector_t& parameters,
             const std::string name = "B-Spline",
             bool clamped = false);

    using Trajectory<N>::impl_compute;

    void impl_compute (result_ref, value_type) const;
    void impl_derivative (derivative_ref g, value_type x, size_type order)
      const;
    void impl_derivative (derivative_ref g, StableTimePoint, size_type order)
      const;

    vector_t basisFunctions (value_type t, size_type order) const
      ROBOPTIM_TRAJECTORY_DEPRECATED;

    /// \brief Compute the basis polynomials.
    void computeBasisPolynomials ();

    /// \brief Derive basis polynomials, but express the result as
    /// polynomials of BSpline<N> instead of BSpline<N-M>.
    /// \tparam M derivation order.
    template <int M>
    basisPolynomialsVector_t deriveBasisPolynomials () const;

    /// \brief Order of the B-Spline
    static const size_type order_ = N;

    /// \brief Generate base polynomial set
    /// For the basic spline formula noted in the pdf from the docs section,
    /// this implements (3), the recursion formula. It returns the results as
    /// factors of b_j,0 where j is the index of the returned map.
    /// \param j current knot index
    /// \param n current basis function order
    /// \return a std::map.
    /// Call this with j as the index of the knot for which spline segment your
    /// want to generate the basis function for, and n as the order of your
    /// spline.
    cox_map cox_de_boor (size_type j, size_type n) const;

    /// \brief Initialize the knot vector.
    /// \param clamped whether the spline is clamped.
    void initializeKnots (bool clamped);

  private:

    /// \brief Number of control points.
    size_type nbp_;

    /// \brief Vector of knots.
    vector_t knots_;

    /// \brief basisPolynomials_[i][j] = B_{i,i+j}
    basisPolynomialsVector_t basisPolynomials_;

    /// \brief Whether the B-spline is uniform.
    /// Note: used for faster interval lookup.
    bool uniform_;

  protected:
    /// \brief Pointer to B-spline logger (see log4cxx documentation).
    static log4cxx::LoggerPtr logger;
  };

  /// \brief LOG4CXX logger for B-splines.
  /// \tparam N order of the B-Spline.
  /// \return logger.
  template <int N>
  log4cxx::LoggerPtr BSpline<N>::logger
  (log4cxx::Logger::getLogger ("roboptim.trajectory"));

  /// @}

} // end of namespace trajectory.
} // end of namespace roboptim.

# include <roboptim/trajectory/b-spline.hxx>
#endif //! ROBOPTIM_TRAJECTORY_B_SPLINE_HH
