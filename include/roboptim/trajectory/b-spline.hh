// Copyright (C) 2009, 2010 by Thomas Moulard, AIST, CNRS, INRIA.
// Copyright (C) 2013 by Alexander Werner, DLR.
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
  /// \addtogroup roboptim_function
  /// @{
  /// B-Spline trajectory.
  /// Implement a B-Spline as a trajectory based in doc/quintic-b-spline.tex
  /// The template parameter N determines the order of the spline.
  template<int N>
  class BSpline : public Trajectory<N>
  {
  public:
		typedef Trajectory<N> parent_t;
		typedef typename parent_t::interval_t interval_t;
		typedef typename parent_t::size_type size_type;
		typedef typename parent_t::value_type value_type;
		typedef typename parent_t::vector_t vector_t;
		typedef typename parent_t::matrix_t matrix_t;
		typedef typename parent_t::result_t result_t;
		typedef typename parent_t::gradient_t gradient_t;
		typedef typename parent_t::jacobian_t jacobian_t;
    /// \brief Instantiate a B-Spline from its definition.
    ///
    /// \param timeRange spline time range: $\f$[t_3,t_n]\f$
    /// \param dimension spline dimension: \f$n\f$
    /// \param parameters vector of parameters defining control points
    /// \param name function title
    ///
    /// Number of control points is inferred from dimension of dimenion of
    /// parameter vector.
    BSpline (interval_t timeRange, size_type dimension,
		  const vector_t& parameters,
		  const std::string name = "B-Spline") throw ();

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
    BSpline (interval_t tr, size_type dimension,
  			      const vector_t& parameters,
  			      const vector_t& knots,
  			      std::string name = "B-Spline") throw ();

    /// \brief Copy constructor.
    /// \param spline spline that will be copied
    BSpline (const BSpline<N>& spline) throw ();

    virtual ~BSpline () throw (){};

    /// \brief Modify spline parameters.
    virtual void setParameters (const vector_t&) throw ();

    virtual jacobian_t variationConfigWrtParam (double t) const throw ();
    virtual jacobian_t variationDerivWrtParam (double t, size_type order)
      const throw ();
    virtual value_type singularPointAtRank (size_type rank) const;
    virtual vector_t derivBeforeSingularPoint (size_type rank, size_type order)
      const;
    virtual vector_t derivAfterSingularPoint (size_type rank, size_type order)
      const;

    ROBOPTIM_IMPLEMENT_CLONE(BSpline<N>)

    virtual Trajectory<N>* resize (interval_t timeRange) const throw ();

    /// \brief Display the function on the specified output stream.
    ///
    /// \param o output stream used for display
    /// \return output stream
    virtual std::ostream& print (std::ostream& o) const throw ();

    jacobian_t variationConfigWrtParam (StableTimePoint tp) const throw ();
    jacobian_t variationDerivWrtParam (StableTimePoint tp, size_type order)
      const throw ();

    value_type Dt () const ROBOPTIM_TRAJECTORY_DEPRECATED;

    vector_t const & knots() const;

  protected:
    using Trajectory<N>::impl_compute;
    void impl_compute (result_t&, double) const throw ();
    void impl_derivative (gradient_t& g, double x, size_type order)
      const throw ();
    void impl_derivative (gradient_t& g, StableTimePoint, size_type order)
      const throw ();
    size_type interval (value_type t) const;

    vector_t basisFunctions (value_type t, size_type order) const
      ROBOPTIM_TRAJECTORY_DEPRECATED;
    void computeBasisPolynomials ();
    /// order of the B-Spline
    static const size_type order_ = N;

    typedef Polynomial<N> polynomial_t;
    typedef Monomial<N> monomial_t;
    typedef std::map<int,polynomial_t> cox_map;
    typedef typename cox_map::iterator cox_map_itr_t;
    cox_map cox_de_boor(size_type j, size_type n) const;
  private:
    /// \brief Number of control points.
    size_type nbp_;
    /// Vector of knots
    vector_t knots_;
    /// basisPolynomials_[i][j] = B_{i,i+j}
    std::vector <std::vector <Polynomial<N> > > basisPolynomials_;
    /// For backward compatibility only
    bool uniform_;
  };

  /// @}

} // end of namespace roboptim.

# include <roboptim/trajectory/b-spline.hxx>
#endif //! ROBOPTIM_TRAJECTORY_B_SPLINE_HH
