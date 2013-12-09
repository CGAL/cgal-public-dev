#ifndef CGAL_BEZIER_TRAITS_CURVE_2_H
#define CGAL_BEZIER_TRAITS_CURVE_2_H

/*! \file
 * Header file for the _Bezier_curve_2 class.
 */

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Arr_geometry_traits/de_Casteljau_2.h>

#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Get_arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_d_1.h>
#include <CGAL/Fraction_traits.h>

#include <algorithm>
#include <deque>
#include <vector>
#include <list>
#include <ostream>

namespace CGAL {

/*! \class _Bezier_curve_2
 * Representation of a Bezier curve, specified by (n+1) control points
 * p_0, ... , p_n that define the curve (X(t), Y(t)) for 0 <= t <= 1,
 * where X(t) and Y(t) are polynomials of degree n.
 *
 * The class is templated with three parameters: 
 * Rat_kernel A geometric kernel, where Alg_kernel::FT is the number type
 *            for the coordinates of control points (and subsequently also for
 *            the polynomial coefficients). This number type must be the same
 *            as Nt_traits::Rational.
 * Alg_kernel A geometric kernel, where Alg_kernel::FT is a number type
 *            for representing algebraic numbers. This number type must be the
 *            same as Nt_traits::Algebraic.
 * Nt_traits A traits class that defines the Integer, Rational and Algebraic
 *           number-types, as well as the Polynomial class (a polynomial with
 *           integer coefficients) and enables performing various operations
 *           on objects of these types.
 * Bounding_traits A traits class for filtering the exact computations.
 */

// Forward declaration:
template <class Algebraic_kernel_d_1_, class RatKernel_, class BoundingTraits_>
class _Bezier_curve_2;

template <class Algebraic_kernel_d_1_, class RatKernel_, class BoundingTraits_>
class _Bezier_curve_2_rep
{
  friend class _Bezier_curve_2<Algebraic_kernel_d_1_, RatKernel_, BoundingTraits_>;

public:

  typedef RatKernel_                              Rat_kernel;
  typedef Algebraic_kernel_d_1_			  Algebraic_kernel_d_1;
//  typedef AlgKernel_                              Alg_kernel;
//  typedef NtTraits_                               Nt_traits;
  typedef BoundingTraits_                         Bounding_traits;


  typedef typename Algebraic_kernel_d_1::Bound						Bound;
  typedef typename Algebraic_kernel_d_1::Algebraic_real_1				Algebraic_real_1;
  typedef typename Algebraic_kernel_d_1::Coefficient					Coefficient

  //typedef typename Nt_traits::Rational            Rational;
  typedef typename Get_arithmetic_kernel<Coefficient>::Arithmetic_kernel::Rational	Rational;

  typedef typename Rat_kernel::Point_2            					Rat_point_2;
  typedef typename std::pair<Algebraic_real_1, Algebraic_real_1>			Alg_point_2;
  //typedef typename Alg_kernel::Point_2            Alg_point_2;

  // TODO Maybe use the Coefficient type for Integer later
  typedef typename Get_arithmetic_kernel<Coefficient>::Arithmetic_kernel::Integer	Integer;
  //typedef typename Nt_traits::Integer             Integer;
  //typedef typename Nt_traits::Algebraic           Algebraic;
  //typedef typename Nt_traits::Polynomial          Polynomial;
  typedef typename Algebraic_kernel_d_1::Polynomial_1					Polynomial_1;

  typedef typename CGAL::Polynomial_traits_d<Polynomial_1>				Polynomial_traits_1;


private:

  typedef std::deque<Rat_point_2>                 Control_point_vec;

  Control_point_vec   _ctrl_pts;      /*!< The control points (we prefer deque
                                           to a vector, as it enables
                                           push_front()). */
  Bbox_2              _bbox;          /*!< A bounding box for the curve. */
  bool                _no_self_inter; /*!< Whether the curve surely has  no
                                           self intersections. */

  /// \name Lazily-evaluated values of the polynomials B(t) = (X(t), Y(t)).
  //@{
 
  // X(t) is given by *p_polyX(t) / _normX:
  mutable Polynomial_1        *p_polyX;       // The polynomial for x.
  mutable Integer             *p_normX;       // Normalizing factor for y.

  // Y(t) is given by _polyY(t) / _normY:
  mutable Polynomial_1        *p_polyY;       // The polynomial for y.
  mutable Integer             *p_normY;       // Normalizing factor for y.
  //@}

public:

  /*! Default constructor. */
  _Bezier_curve_2_rep () : 
    _no_self_inter (true),
    p_polyX(NULL),
    p_normX(NULL),
    p_polyY(NULL),
    p_normY(NULL)
  {}

  /*! Copy constructor (isn't really used). */
  _Bezier_curve_2_rep (const _Bezier_curve_2_rep& other) :
    _ctrl_pts(other._ctrl_pts),
    _bbox(other._bbox),
    _no_self_inter(other._no_self_inter),
    p_polyX(NULL),
    p_normX(NULL),
    p_polyY(NULL),
    p_normY(NULL)
  {
    if (other.p_polyX != NULL)
      p_polyX = new Polynomial_1 (*(other.p_polyX));
    if (other.p_polyY != NULL)
      p_polyY = new Polynomial_1 (*(other.p_polyY));
    if (other.p_normX != NULL)
      p_normX = new Integer(*(other.p_normX));
    if (other.p_normY != NULL)
      p_normY = new Integer(*(other.p_normY));
  }

  /*!
   * Constructor from a given range of control points.
   * \param pts_begin An iterator pointing to the first point in the range.
   * \param pts_end A past-the-end iterator for the range.
   * \pre The value-type of the input iterator must be Rat_kernel::Point_2.
   *      It is forbidden to specify two identical consecutive control points.
   */


  template <class InputIterator>
  _Bezier_curve_2_rep (InputIterator pts_begin, InputIterator pts_end) :
    p_polyX(NULL),
    p_normX(NULL),
    p_polyY(NULL),
    p_normY(NULL)
  {
    // Copy the control points and compute their bounding box.
    const int   pts_size = std::distance (pts_begin, pts_end);
    double      x, y;
    double      x_min = 0, x_max = 0;
    double      y_min = 0, y_max = 0;
    int         k;

    CGAL_precondition_msg (pts_size > 1, "There must be at least 2 control points.");

    _ctrl_pts.resize (pts_size);
    
    for (k = 0; pts_begin != pts_end; ++pts_begin, k++)
    {
//SL: Acccording to the fact that all operations are based on polynomials
//    duplicated control points can be allowed.
//      // Make sure that we do not have two identical consecutive control
//      // points.
//      CGAL_precondition_msg
//          (k == 0 || ! equal (*pts_begin, _ctrl_pts[k - 1]),
//           "Two consecutive control points must not be identical.");

      // Copy the current control point.
      _ctrl_pts[k] = *pts_begin;

      // Update the bounding box, if necessary.
      // pts_begin references to Rat_point_2 type
      x = CGAL::to_double (pts_begin->x());
      y = CGAL::to_double (pts_begin->y());

      if (k == 0)
      {
        x_min = x_max = x;
        y_min = y_max = y;
      }
      else
      {
        if (x < x_min)
          x_min = x;
        else if (x > x_max)
          x_max = x;

        if (y < y_min)
          y_min = y;
        else if (y > y_max)
          y_max = y;
      }
    }

    // Construct the bounding box.
    _bbox = Bbox_2 (x_min, y_min, x_max, y_max);

    // Use the bounding traits to determine whether the curve surely has
    // not self intersections.
    Bounding_traits     bound_tr;

    _no_self_inter = ! bound_tr.may_have_self_intersections (_ctrl_pts);
  }

  /*! Destructor. */
  ~_Bezier_curve_2_rep ()
  {
    if (p_polyX != NULL) 
      delete p_polyX;
    if (p_normX != NULL) 
      delete p_normX;
    if (p_polyY != NULL) 
      delete p_polyY;
    if (p_normY != NULL) 
      delete p_normY;
  }

  /// \name Access the polynomials (lazily evaluated).
  //@{

  /*! Check if the polynomials are already constructed. */
  bool has_polynomials () const
  {
    return (p_polyX != NULL && p_normX != NULL &&
            p_polyY != NULL && p_normY != NULL);
  }

  /*! Get the polynomial X(t). */
  const Polynomial_1& x_polynomial () const 
  {
    if (p_polyX == NULL)
      _construct_polynomials ();
    
    return (*p_polyX);
  }

  /*! Get the normalizing factor for X(t). */
  const Integer& x_norm () const 
  {
    if (p_normX == NULL)
      _construct_polynomials ();
    
    return (*p_normX);
  }

  /*! Get the polynomial Y(t). */
  const Polynomial_1& y_polynomial () const 
  {
    if (p_polyY == NULL)
      _construct_polynomials ();

    return (*p_polyY);
  }

  /*! Get the normalizing factor for Y(t). */
  const Integer& y_norm () const 
  {
    if (p_normY == NULL)
      _construct_polynomials ();
    
    return (*p_normY);
  }
  //@}

private:

  /*!
   * Construct the representation of X(t) and Y(t).
   * The function is declared as "const" as it changes only mutable members.
   */
  void _construct_polynomials () const;

  /*!
   * Compute the value of n! / (j! k! (n-k-j)!).
   */
  Integer _choose (int n, int j, int k) const;

};

template <class Algebraic_kernel_d_1_, class RatKernel_, class BoundingTraits_>
class _Bezier_curve_2 :
  public Handle_for<_Bezier_curve_2_rep<Algebraic_kernel_d_1_, RatKernel_, BoundingTraits_> >
{
public:

  typedef RatKernel_                              Rat_kernel;
  typedef Algebraic_kernel_d_1_			  Algebraic_kernel_d_1;
//  typedef AlgKernel_                              Alg_kernel;
//  typedef NtTraits_                               Nt_traits;
  typedef BoundingTraits_                         Bounding_traits;

  typedef _Bezier_curve_2<Algebraic_kernel_d_1, Rat_kernel, Bounding_traits>        Self;

private:

  typedef _Bezier_curve_2_rep<Algebraic_kernel_d_1, Rat_kernel, Bounding_traits>    Bcv_rep;
  typedef Handle_for<Bcv_rep>                     Bcv_handle;

  typedef typename Bcv_rep::Control_point_vec     Control_pt_vec;

public:

  typedef typename Bcv_rep::Rat_point_2           Rat_point_2;
  typedef typename Bcv_rep::Alg_point_2           Alg_point_2;
  
  typedef typename Bcv_rep::Integer               Integer;
  typedef typename Bcv_rep::Rational              Rational;
  typedef typename Bcv_rep::Bound	          Bound;
  typedef typename Bcv_rep::Algebraic_real_1      Algebraic_real_1;
  typedef typename Bcv_rep::Polynomial_1          Polynomial_1;
  typedef typename Bcv_rep::Polynomial_traits_1	  Polynomial_traits_1;

  typedef typename Control_pt_vec::const_iterator Control_point_iterator;

public:

  /*!
   * Default constructor.
   */
  _Bezier_curve_2 () :
    Bcv_handle (Bcv_rep())
  {}

  /*!
   * Copy constructor.
   */
  _Bezier_curve_2 (const Self& bc) :
    Bcv_handle (bc)
  {}

  /*!
   * Constructor from a given range of control points.
   * \param pts_begin An iterator pointing to the first point in the range.
   * \param pts_end A past-the-end iterator for the range.
   * \pre The value-type of the input iterator must be Rat_kernel::Point_2.
   *      It is forbidden to specify two identical consecutive control points.
   */
  template <class InputIterator>
  _Bezier_curve_2 (InputIterator pts_begin, InputIterator pts_end) :
    Bcv_handle (Bcv_rep (pts_begin, pts_end))
  {}

  /*!
   * Assignment operator.
   */
  Self& operator= (const Self& bc)
  {
    if (this == &bc || this->identical (bc))
      return (*this);

    Bcv_handle::operator= (bc);
    return (*this);
  }
  
  /*!
   * Get a unique curve ID (based on the actual representation pointer).
   */
  size_t id () const
  {
    return (reinterpret_cast<size_t> (this->ptr()));
  }

  /*!
   * Get the polynomial for the x-coordinates of the curve.
   */
  const Polynomial_1& x_polynomial () const
  {
    return (this->_rep().x_polynomial());
  }

  /*!
   * Get the normalizing factor for the x-coordinates.
   */
  const Integer& x_norm () const 
  {
    return (this->_rep().x_norm());
  }

  /*!
   * Get the polynomial for the y-coordinates of the curve.
   */
  const Polynomial_1& y_polynomial () const 
  {
    return (this->_rep().y_polynomial());
  }

  /*!
   * Get the normalizing factor for the y-coordinates.
   */
  const Integer& y_norm () const 
  {
    return (this->_rep().y_norm());
  }

  /*!
   * Get the number of control points inducing the Bezier curve.
   */
  unsigned int number_of_control_points () const
  {
    return (this->_rep()._ctrl_pts.size());
  }

  /*!
   * Get the i'th control point.
   * \pre i must be between 0 and n - 1, where n is the number of control
   *      points.
   */
  const Rat_point_2& control_point (unsigned int i) const
  {
    CGAL_precondition (i < number_of_control_points());

    return ((this->_rep()._ctrl_pts)[i]);
  }

  /*!
   * Get an interator for the first control point.
   */
  Control_point_iterator control_points_begin () const
  {
    return (this->_rep()._ctrl_pts.begin());
  }

  /*!
   * Get a past-the-end interator for control points.
   */
  Control_point_iterator control_points_end () const
  {
    return (this->_rep()._ctrl_pts.end());
  }

  /*!
   * Check if both curve handles refer to the same object.
   */
  bool is_same (const Self& bc) const
  {
    return (this->identical (bc));
  }

  /*!
   * Compute a point of the Bezier curve given a rational t-value.
   * \param t The given t-value.
   * \return The point B(t).
   */
// TODO Bound data type might not work in this case
  Rat_point_2 operator() (const Bound& t) const;
  
  /*!
   * Compute a point of the Bezier curve given an algebraic t-value.
   * \param t The given t-value.
   * \return The point B(t).
   */
  Alg_point_2 operator() (const Algebraic_real_1& t) const;
 
  /*!
   * Sample a portion of the curve (for drawing purposes, etc.).
   * \param t_start The t-value to start with.
   * \param t_end The t-value to end at.
   * \param n_samples The required number of samples.
   * \param oi Output: An output iterator for the samples. The value-type of
   *                   this iterator must be std::pair<double, double>.
   * \return A past-the-end iterator for the samples.
   */
  template <class OutputIterator>
  OutputIterator sample (const double& t_start, const double& t_end,
                         unsigned int n_samples,
                         OutputIterator oi) const 
  {
    // Convert the coordinates of the control points to doubles.
    typedef Simple_cartesian<double>                        App_kernel;
    typedef App_kernel::Point_2                             App_point_2;

    const unsigned int             n_pts = number_of_control_points();
    std::vector<App_point_2>       app_ctrl_pts (n_pts);
    unsigned int                   k;

    // .x() and .y() called for rational_point type
    for (k = 0; k < n_pts; k++)
    {
      const Rat_point_2&   pt = control_point(k);

      app_ctrl_pts[k] = App_point_2 (CGAL::to_double (pt.x()), CGAL::to_double (pt.y()));
    }

    // Sample the approximated curve.
    const unsigned int   n = (n_samples >= 2) ? n_samples : 2; 
    const double         delta_t = (t_end - t_start) / (n - 1);
    App_point_2          p;

    for (k = 0; k < n; k++)
    {
      p = point_on_Bezier_curve_2 (app_ctrl_pts.begin(), app_ctrl_pts.end(),
                                   t_start + k * delta_t);

      *oi = std::make_pair (p.x(), p.y());
      ++oi;
    }

    return (oi);
  }

  /*!
   * Compute all parameter values t such that the x-coordinate of B(t) is x0.
   * Note that the function does not return only values between 0 and 1, so
   * the output t-values may belong to the imaginary continuation of the curve.
   * \param x0 The given x-coordinate.
   * \param oi Output: An output iterator for the t-values.
   * \return A past-the-end iterator.
   */
  template <class OutputIterator>
  OutputIterator get_t_at_x (const Bound& x0,
                             OutputIterator oi) const
  {
    return (_solve_t_values (this->_rep().x_polynomial(), 
                             this->_rep().x_norm(), 
                             x0,
                             oi));
  }

  /*!
   * Compute all parameter values t such that the y-coordinate of B(t) is y0.
   * Note that the function does not return only values between 0 and 1, so
   * the output t-values may belong to the imaginary continuation of the curve.
   * \param y0 The given y-coordinate.
   * \param oi Output: An output iterator for the t-values.
   * \return A past-the-end iterator.
   */
  template <class OutputIterator>
  OutputIterator get_t_at_y (const Bound& y0,
                             OutputIterator oi) const
  {
    return (_solve_t_values (this->_rep().y_polynomial(), 
                             this->_rep().y_norm(), y0,
                             oi));
  }

  /*!
   * Check if the two curves have the same support.
   */
  bool has_same_support (const Self& bc) const;

  /*!
   * Get the bounding box of the curve.
   */
  const Bbox_2& bbox () const
  {
    return (this->_rep()._bbox);
  }

  /*!
   * Check if the curve contains not self intersections.
   * Note that there may not be any self intersections even if the
   * function returns true (but not vice versa).
   */
  bool has_no_self_intersections () const
  {
    return (this->_rep()._no_self_inter);
  }

private:

  // Get the representation.
  inline const Bcv_rep& _rep () const
  {
    return (*(this->ptr()));
  }

  inline Bcv_rep& _rep ()
  {
    return (*(this->ptr()));
  }

  /*!
   * Compute all parameter values t, such that P(t) = val.
   * \param poly The polynomial.
   * \param norm Its normalizing factor.
   * \param val The required value.
   * \param oi Output: An output iterator for the t-values.
   * \return A past-the-end iterator.
   */
  template <class OutputIterator>
  OutputIterator _solve_t_values (const Polynomial_1& poly,
                                  const Integer& norm,
                                  const Bound& val,
                                  OutputIterator oi) const
  {
    // Construct the polynomial P(t) - val = 0:

    typename Polynomial_traits_1::Degree		degree;
    typename Polynomial_traits_1::Get_coefficient	get_coefficient;
    typename Polynomial_traits_1::Construct_polynomial	construct_polynomial;


    // TODO Check if fraction_trait should be on rational type

    typedef CGAL::Fraction_traits<Bound>	FT_bound;
    typedef FT_bound::Numerator_type		NT;
    typedef FT_bound::Denominator_type		DT;

    //Nt_traits             nt_traits;
    // TODO Check the return type for the NULLTAG and conversion

    //NT numer; DT denom;
    Integer numer, denom;
    FT_bound::Decompose() (val, numer, denom);
    //const Integer         numer = nt_traits.numerator (val);
    //const Integer         denom = nt_traits.denominator (val);

    //const int             deg = nt_traits.degree (poly);
    const int             deg = degree (poly);
    if (deg <=0 ) 
    	return oi;

    Integer              *coeffs = new Integer [deg + 1];
    int                   k;

    for (k = 1; k <= deg; k++)
    {
      //coeffs[k] = nt_traits.get_coefficient (poly, k) * denom;
      coeffs[k] = get_coefficient (poly, k) * denom;
    }
    //coeffs[0] = nt_traits.get_coefficient (poly, 0) * denom - numer * norm;
    coeffs[0] = get_coefficient (poly, 0) * denom - numer * norm;

    // Solve the polynomial and obtain the t-values.
    typename Algebraic_kernel_d algebraic_kernel_d_1;
    typename Algebraic_kernel_d_1::Solve_1 solve_1 = algebraic_kernel_d_1.solve_1_object();
    typename std::list<Algebraic_real_1, Algebraic_real_1> root_pair;

    solve_1 (construct_polynomial (coeffs, coeffs + deg), std::back_inserter(root_pair));

    typename std::list<Algebraic_real_1, Algebraic_real_1>::iterator rit = root_pair.begin();
    while(rpit != root_pair.end())
    {
    	*oi++ = rpit->first;
	rpit++;
    }

    //OutputIterator end = oi;

    //OutputIterator  end = nt_traits.compute_polynomial_roots (nt_traits.construct_polynomial (coeffs, deg), oi);

    delete[] coeffs;
    return oi;
    //return (end);
  }
};

/*!
 * Exporter for Bezier curves.
 */
template <class Algebraic_kernel_d_1_, class Rat_kernel, class Bounding_traits>
std::ostream& 
operator<< (std::ostream& os, 
            const _Bezier_curve_2<Algebraic_kernel_d_1_, Rat_kernel, Bounding_traits> & bc)
{
  const unsigned int  n = bc.number_of_control_points();
  unsigned int        k;

  os << n;
  for (k = 0; k < n; k++)
    os << "  " << bc.control_point(k);

  return (os);
}

/*!
 * Importer for Bezier curves.
 */
template <class Algebraic_kernel_d_1_, class Rat_kernel, class Bounding_traits>
std::istream& 
operator>> (std::istream& is, 
            _Bezier_curve_2<Algebraic_kernel_d_1_, Rat_kernel, Bounding_traits> & bc)
{
  // Read the number of control points.
  unsigned int  n;

  is >> n;

  // Read the control points.
  std::vector<typename Rat_kernel::Point_2>   ctrl_pts (n);
  unsigned int                                k;

  for (k = 0; k < n; k++)
    is >> ctrl_pts[k];

  // Set the Bezier curve.
  bc = _Bezier_curve_2<Algebraic_kernel_d_1_, Rat_kernel, Bounding_traits> (ctrl_pts.begin(), ctrl_pts.end());

  return (is);
}


// ---------------------------------------------------------------------------
// Construct the representation of X(t) and Y(t).
//
template <class Algebraic_kernel_d_1_, class RatKer, class BndTrt>
void _Bezier_curve_2_rep<Algebraic_kernel_d_1_, RatKer, BndTrt>::_construct_polynomials () const
{
  const int        n = _ctrl_pts.size() - 1;
  Bound        *coeffsX = new Bound [n + 1];
  Bound        *coeffsY = new Bound [n + 1];
  const Bound   rat_zero = Bound (0);
  int              j, k;
  
  CGAL_precondition_msg (n > 0,
                         "There must be at least 2 control points.");
  
  for (j = 0; j <= n; j++)
    coeffsX[j] = coeffsY[j] = rat_zero;
  
  // Compute the rational coefficients, given by the formula:
  //
  //                     n
  //                   *****
  //                   *   *      / n \   k        n-k
  //   (X(t), Y(t)) =   *    p_k (     ) t  (1 - t)
  //                   *   *      \ k /
  //                   *****
  //                    k=0
  //
  Bound                  px, py;
  Integer                n_over_k_j;
  bool                   even_exp;
  
  typename Control_point_vec::const_iterator pts_begin = _ctrl_pts.begin();  
  typename Control_point_vec::const_iterator pts_end = _ctrl_pts.end();  

  for (k = 0; pts_begin != pts_end; ++pts_begin, k++)
  {
    px = pts_begin->x();
    py = pts_begin->y();
    
    // By simplifying (1 - t)^(n-k) we obtain that the k'th expression of
    // the above sum is given by:
    //
    //     n-k
    //    *****
    //    *   *     j            n!         j+k
    //     *    (-1)  p_k ---------------- t
    //    *   *            j! k! (n-k-j)!
    //    *****
    //     j=0
    //
    even_exp = true; 
    for (j = 0; j <= n - k; j++)
    {
      n_over_k_j = _choose (n, k, j);
      
      if (even_exp)
      {
        // We should add the current values to the coefficients of the
        // monomial t^(n_j).
        coeffsX[j + k] += px * n_over_k_j;
        coeffsY[j + k] += py * n_over_k_j;
      }
      else
      {
        // We should subtract the current values from the coefficients of the
        // monomial t^(n_j).
        coeffsX[j + k] -= px * n_over_k_j;
        coeffsY[j + k] -= py * n_over_k_j;
      }
      
      // As we increment j, negate the "even" flag for the exponent (n-j).
      even_exp = !even_exp;
    } // loop on j.
  } // loop on k.
  
  // Convert the rational polynomials to polynomials with rational
  // coefficients (plus normalizing factors).
  //Nt_traits        nt_traits;
  p_polyX = new Polynomial_1;
  p_normX = new Integer;
  p_polyY = new Polynomial_1;
  p_normY = new Integer;
  
// TODO Check for this kind of polynomial construction
// TODO Need the new functionality
  nt_traits.construct_polynomial (coeffsX, n,
                                  *p_polyX, *p_normX);
  delete[] coeffsX;
  
  nt_traits.construct_polynomial (coeffsY, n,
                                  *p_polyY, *p_normY);
  delete[] coeffsY;

//  CGAL_assertion (nt_traits.degree (*p_polyX) >= 0);
//  CGAL_assertion (nt_traits.degree (*p_polyY) >= 0);

  return;
}

// ---------------------------------------------------------------------------
// Compute the value of n! / (j! k! (n-k-j)!).
//
template <class Algebraic_kernel_d_1_, class RatKer, class BndTrt>
typename _Bezier_curve_2_rep<Algebraic_kernel_d_1_, RatKer, BndTrt>::Integer
_Bezier_curve_2_rep<Algebraic_kernel_d_1_, RatKer, BndTrt>::_choose (int n, int j, int k) const
{
  Integer   reduced_fact = 1;
  Integer   j_fact = 1, k_fact = 1;
  int       i;
  
  for (i = n - k - j + 1; i <= n; i++)
    reduced_fact *= Integer (i);
  
  for (i = 2; i <= j; i++)
    j_fact *= Integer (i);
  
  for (i = 2; i <= k; i++)
    k_fact *= Integer (i);

  //TODO Check the div functionality
  return (CGAL::div (reduced_fact, (j_fact * k_fact)));
}

// ---------------------------------------------------------------------------
// Compute a point on the Bezier curve B(t) given a rational t-value.
//
template <class Algebraic_kernel_d_1_, class RatKer, class BndTrt>
typename _Bezier_curve_2<Algebraic_kernel_d_1_, RatKer, BndTrt>::Rat_point_2
_Bezier_curve_2<Algebraic_kernel_d_1_, RatKer, BndTrt>::operator() (const Bound& t) const 
{
  // Check for extermal t values (either 0 or 1).
  const CGAL::Sign   sign_t = CGAL::sign (t);

  CGAL_precondition (sign_t != NEGATIVE);

  if (sign_t == ZERO)
  {
    // If t is 0, simply return the first control point.
    return (this->_rep()._ctrl_pts[0]);
  }


  typename Algebraic_kernel_d_1 alg_kernel_d_1;
  typename Algebraic_kernel_d_1::Compare_1 compare_1 = alg_kernel_d_1.compare_1_object();
  typename Algebraic_kernel_d_1::Construct_algebraic_real_1 construct_algebraic_real_1 = alg_kernel_d_1.construct_algebraic_real_1_object();


  const Algebraic_real_1& temp_t = construct_algebraic_real_1 (t);
  Comparison_result  res = compare(temp_t, 1);
  //Comparison_result  res = CGAL::compare (t, Rational(1));

  CGAL_precondition (res != LARGER);

  if (res == EQUAL)
  {
    // If t is 1, simply return the first control point.
    return (this->_rep()._ctrl_pts[this->_rep()._ctrl_pts.size() - 1]);
  }

  // Evaluate the point for 0 < t < 1.
  if (! this->_rep().has_polynomials())
  {
    // Use de Casteljau's algorithm to evaluate the polynomial at t.
    return (point_on_Bezier_curve_2 (this->_rep()._ctrl_pts.begin(),
                                     this->_rep()._ctrl_pts.end(),
                                     t));
  }

  // Compute the x and y coordinates using the X(t), Y(t) polynomials.
  Bound           x, y;

  typename Polynomial_traits_1::Evaluate	evaluate;

  //Nt_traits          nt_traits;
    
  //x = nt_traits.evaluate_at (this->_rep().x_polynomial(), t) / Rational (this->_rep().x_norm(), 1);

  // TODO Construct uses conversion using CGAL::to_double
  // TODO Same problem as that of Bezier_traits_cache for computation of the exact value of a polynomial
  x = Bound (CGAL::to_double( evaluate (this->_rep().x_polynomial(), t)) / CGAL::to_double(this->_rep().x_norm()) );
  //y = nt_traits.evaluate_at (this->_rep().y_polynomial(), t) / Rational (this->_rep().y_norm(), 1);
  y = Bound (CGAL::to_double( evaluate (this->_rep().y_polynomial(), t)) / CGAL::to_double(this->_rep().y_norm()) );
  
  // Return the point.
  // TODO Maybe need to convert to Algebraic_real_1
  // OR Maybe use the Rational_traits and use make_rational(x,y)
  return (Rat_point_2 (x, y));
}



// ---------------------------------------------------------------------------
// Compute a point on the Bezier curve B(t) given an algebraic t-value.
//
template <class Algebraic_kernel_d_1_, class RatKer, class BndTrt>
typename _Bezier_curve_2<Algebraic_kernel_d_1_, RatKer, BndTrt>::Alg_point_2
_Bezier_curve_2<Algebraic_kernel_d_1_, RatKer, BndTrt>::operator() (const Algebraic_real_1& t) const
{
  // Check for extermal t values (either 0 or 1).
  //Nt_traits          nt_traits;
  const CGAL::Sign   sign_t = CGAL::sign (t);

  CGAL_precondition (sign_t != NEGATIVE);

  if (sign_t == ZERO)
  {
    // If t is 0, simply return the first control point.
    const Rat_point_2&  p_0 = this->_rep()._ctrl_pts[0];
    
    // TODO maybe need to create the Algebraic_real_1 data type using CGAL::to_double
    //return (Alg_point_2 (nt_traits.convert (p_0.x()), nt_traits.convert (p_0.y())));
    return (std::make_pair ( p_0.x(), p_0.y() ));
  }

  typename Algebraic_kernel_d_1 alg_kernel_d_1;
  typename Algebraic_kernel_d_1::Compare_1 compare_1 = alg_kernel_d_1.compare_1_object();
  typename Algebraic_kernel_d_1::Construct_algebraic_real_1 construct_algebraic_real_1 = alg_kernel_d_1.construct_algebraic_real_1_object();
  Comparison_result	res = compare_1(t, 1);
  //Comparison_result  res = CGAL::compare (t, Algebraic(1));

  CGAL_precondition (res != LARGER);

  if (res == EQUAL)
  {
    // If t is 1, simply return the first control point.
    const Rat_point_2&  p_n = 
      this->_rep()._ctrl_pts[this->_rep()._ctrl_pts.size() - 1];
    
    //return (Alg_point_2 (nt_traits.convert (p_n.x()), nt_traits.convert (p_n.y())));
    return ( std::make_pair( p_n.x(), p_n.y() ) );
  }

  // The t-value is between 0 and 1: Compute the x and y coordinates.


  // TODO Check for the convert function
  typename Polynomial_traits_1::Evaluate evaluate;
  //const Algebraic    x = nt_traits.evaluate_at (this->_rep().x_polynomial(), t)/ nt_traits.convert (this->_rep().x_norm());
  const Algebraic_real_1    x = construct_algebraic_real_1 ( evaluate (this->_rep().x_polynomial(), CGAL::to_double(t)) / (this->_rep().x_norm()) );
  //const Algebraic    y = nt_traits.evaluate_at (this->_rep().y_polynomial(), t)/ nt_traits.convert (this->_rep().y_norm());
  const Algebraic_real_1    y = construct_algebraic_real_1 ( evaluate (this->_rep().y_polynomial(), CGAL::to_double(t)) / (this->_rep().y_norm()) );

  // Return the Alg_point_2 type which is std::pair<Algebraic_real_1, Algebraic_real_1>
  //return (Alg_point_2 (x, y));
  return ( std::make_pair (x,y) );
}

// ---------------------------------------------------------------------------
// Check if the two curves have the same support.
//
template <class Algebraic_kernel_d_1_, class RatKer, class BndTrt>
bool _Bezier_curve_2<Algebraic_kernel_d_1_, RatKer, BndTrt>::has_same_support
        (const Self& bc) const
{
  // If one curve is of degree d1 and the other of degree d2, there can be
  // at most d1*d2 intersection points between them.
  const int      deg1 = number_of_control_points() - 1;
  const int      deg2 = bc.number_of_control_points() - 1;
  const int      n_samples = deg1*deg2;
  Rat_point_2    p1;
  int            k;
  
  for (k = 0; k <= n_samples; k++)
  {
    // Compute p1 = B1(k/n_samples), where B1 is (*this) curve.
    if (k == 0)
      p1 = (this->_rep()._ctrl_pts[0]);
    else if (k == 1)
      p1 = (this->_rep()._ctrl_pts[this->_rep()._ctrl_pts.size() - 1]);
    else
      p1 = this->operator() (Bound (k, n_samples));
    
    // Get all t-values such that the x-coordinate of B2(t) equals x1,
    // and check if there exists a t-value such that the y-coordinate of
    // b2(t) equals the y-coordinate of p1.
    std::list<Algebraic_real_1>                           t_vals;
    typename std::list<Algebraic_real_1>::const_iterator  t_iter;

    typename Algebraic_kernel_d_1 alg_kernel_d_1;
    typename Algebraic_kernel_d_1::Construct_algebraic_real_1 construct_algebraic_real_1 = alg_kernel_d_1.construct_algebraic_real_1_object();
    typename Algebraic_kernel_d_1::Compare_1 compare_1 = alg_kernel_d_1.compare_1_object();

    //Nt_traits                             nt_traits;
    // TODO Conversion of Rat_point_2 with CGAL::to_double
    const Algebraic_real_1&               y1 = construct_algebraic_real_1( CGAL::to_double(p1.y()) );
    //const Algebraic_real_1&               y1 = nt_traits.convert (p1.y());

    bool                                  eq_y = false;
    
    bc.get_t_at_x (p1.x(), std::back_inserter(t_vals));
    
    for (t_iter = t_vals.begin(); t_iter != t_vals.end(); ++t_iter)
    {

      const Alg_point_2&  p2 = bc (*t_iter);
      
      //if (CGAL::compare (y1, p2.y()) == CGAL::EQUAL)
      if (compare_1 (y1, p2.second) == CGAL::EQUAL)
      {
        eq_y = true;
        break;
      }
    }
    
    // If we found a point on B1 which is not of B2, the two curves do not
    // have the same support.
    if (! eq_y)
      return (false);
  }
  
  // If we reached here, we found (d1*d2 + 1) common points of B1 and B2.
  // This means they have the same support.
  return (true);
}

} //namespace CGAL

#endif
