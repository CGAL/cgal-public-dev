// Copyright (c) 2005-2007 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id:
// 
//
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>
//
//

/*! \file
  This is the file containing a traits class to create anistropic voronoi 
  diagram using the lower envelope code.
  The traits class is based on the CKvA.
 */


#ifndef CGAL_ANISOTROPIC_VORONOI_TRAITS_2_H
#define CGAL_ANISOTROPIC_VORONOI_TRAITS_2_H

#include <CGAL/basic.h>

#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Envelope_voronoi_traits_2/Algebraic_apollonius_traits_base_2.h>

namespace CGAL {

//! \todo Move to Algebraic_structure_traits
/*@{*/
template <typename T>
bool is_zero(const boost::numeric::interval<T> &num)
{
  return (num == boost::numeric::interval<T>(0));
}

template <typename T>
bool is_one(const boost::numeric::interval<T> &num)
{
  return (num == boost::numeric::interval<T>(1));
}

/*@}*/


namespace Envelope_voronoi_traits_2
{
  /*!
   * \class The class represents an anistropic voronoi diagram site.
   *        The distance function from an anistropic site is:
   *        
   *        $d(x) = (x - p)^t * M * (x - p) - D$
   *
   *        Where M is symmetric positive definite matrix.
   *
   *        The class members are:
   *        _center The center of the anistropic site (p in the above formula).
   *        _A, _B, and _C - The values of the anistropic matrix:
   *               / _A _B \
   *        M =    \ _B _C /
   *        _D - The D in the above formula.
   *        
   */  
  template <typename Rational>
  class _Anisotropic_site
  {
  public:
    typedef std::pair< Rational, Rational> Rat_point_2;
    
  protected:
    Rat_point_2 _center;
    Rational _A;
    Rational _B;
    Rational _C;
    Rational _D;
    
  public:
    _Anisotropic_site(const Rat_point_2& center, 
                      const Rational &A, 
                      const Rational &B, 
                      const Rational &C,
                      const Rational &D)
      : _center(center), _A(A), _B(B), _C(C), _D(D)
      {}
    
    const Rat_point_2& center() const
      {
        return _center;
      }
    
    const Rational& A() const
      {
        return _A;
      }
    
    const Rational& B() const
      {
        return _B;
      }
    
    const Rational& C() const
      {
        return _C;
      }
    
    const Rational& D() const
      {
        return _D;
      }
    
  };
  
} // Envelope_voronoi_2

/*!
 * \class The traits class
 */
template <class CurvedKernel_>
class Anisotropic_voronoi_traits_2 : 
: public public Algebraic_apollonius_traits_base_2< 
             CurvedKernel_, 
             Anisotropic_voronoi_traits_2< CurvedKernel_ >,
             typename Envelope_voronoi_traits_2::_Anisotropic_site< 
               typename CurvedKernel_::Curve_kernel_2::         
                 Coordinate_1::Rational
             >,
             CGAL::Field_tag
           >
{
 public:
  typedef CurvedKernel_                                     Curved_kernel_2;
  
  typedef typename Curved_kernel_2::Curve_kernel_2          Algebraic_kernel_2;
  typedef typename Algebraic_kernel_2::Coefficient          Coefficient;
  typedef typename Algebraic_kernel_2::
  Coordinate_1::Rational                                    Rational;
  
  typedef typename CGAL::
  Get_arithmetic_kernel<Coefficient>::
  Arithmetic_kernel::Field_with_sqrt                        Field_with_sqrt;
  
  typedef boost::numeric::interval<Field_with_sqrt>         Interval;
  
public:

  typedef Envelope_voronoi_traits_2::
    _Anisotropic_site<Rational>                             Site_2;

  typedef Anisotropic_voronoi_traits_2<Curved_kernel_2>     Self;

  typedef Algebraic_apollonius_traits_base_2< 
    Curved_kernel_2, Self, Site_2, CGAL::Field_tag >        Base;

  typedef typename Base::Point_2                            Point_2;
  typedef typename Base::X_monotone_curve_2                 X_monotone_curve_2;
  typedef typename Base::Curve_2                            Curve_2;
  typedef typename Base::Multiplicity                       Multiplicity;

  typedef typename Base::Surface_3                          Surface_3;
  typedef typename Base::Xy_monotone_surface_3            Xy_monotone_surface_3;
  

protected:
  
  //! Compute the distance from a point to a site.
  /*! Compute the distance from a point to a site
    \param px The x-coordinate of the point.
    \param py The y-coordinate of the point.
    \param site The site we compute the distance to.
    \return The distance from the point to this site. (+ output argument).
  */
  template <typename NT>
  static NT distance(const NT &px,
                     const NT &py, const Site_2& site)
    {
      
      /*!
       * Computes the distance from the anistropic site.
       * The distance from an anistropic diagram is defined as:
       * 
       * $d(x) = (x - p)^t * M * (x - p) - D$
       *
       * \param p The point to calculate the distance from.
       */
      NT x = px - site.center().first;
      nt y = py - site.center().second;
      
      NT temp_x = A * x + B * y;
      NT temp_y = B * x + C * y;
        
      NT dist = temp_x * x + temp_y * y - D;
      return dist;
    }

public:

  //! Make_xy_monotone_3 functor required by the EnvelopeTraits_3 concept.
  /*! Make_xy_monotone_3
    The traits class treats each segment as a set of 3~different Voronoi sites:
    the two vertices of the segment and the open set that is the interior of
    the segment. (A ray is treated as 2~different sites.)
    The handling of all predicates, specifically predicates related to 
    bisector construction, are made much simpler with the cost of an enlarged
    set of sites.
    The split of a site into several sites is implemented in this functor.
  */
  class Make_xy_monotone_3
  {
  public:
    template <class OutputIterator>
    OutputIterator operator() (Surface_3 S,
                               bool is_lower, 
                               OutputIterator oi) const
      {
        *oi++ = S;
        return oi;
      }
  };

  Make_xy_monotone_3
  make_xy_monotone_3_object () const
    {
      return Make_xy_monotone_3();
    }

  //! Construct_projected_intersections_2 class constructs the bisector curve 
  //  of two sites.
  /*! After the splitting we only need to consider three types of bisectors,
    namely, the bisector between two points, the bisector between two
    lines (where a line can also relate to the interior of a segment or a
    ray) that is a two-line bisector, and the bisector between a point and
    a line that is a full parabola. 
  */
  class Construct_projected_intersections_2
  {
    typedef typename Site_2::Point_2                       KPoint_2;
    
  protected:
    Self *m_traits;
    
  public:
    
    Construct_projected_intersections_2(Self *traits)
      : m_traits(traits)
      {}
    
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s1,
                              const Xy_monotone_surface_3& s2,
                              OutputIterator o) const
      {
        // Construct the bisector according to the equation.
        Rational& x1 = s1.center().x();
        Rational& y1 = s1.center().y();
        Rational& A1 = s1.A();
        Rational& B1 = s1.B();
        Rational& C1 = s1.C();
        Rational& D1 = s1.D();
        
        Rational& x2 = s2.center().x();
        Rational& y2 = s2.center().y();
        Rational& A2 = s2.A();
        Rational& B2 = s2.B();
        Rational& C2 = s2.C();
        Rational& D2 = s2.D();
        
        
        Rational& r = A1 - A2;
        Rational& s = C1 - C2;
        Rational& t = 2 * (B1 - B2);
        Rational& u = -2*(A1*x1 + B1*y1) + 2*(A2*x2 + B2*y2);
        Rational& v = -2*(B1*x1 + C1*y1) + 2*(B2*x2 + C2*y2);
        Rational& w = (A1*x1*x1 + 2*B1*x1*y1 + C1*y1*y1 - D1) - 
          (A2*x2*x2 + 2*B2*x2*y2 + C2*y2*y2 - D2);
        
        
        // construct the degenerate hyperbola.
        return construct_bisector_quadratic(r, s, t, u, v, w, o);
      }
    
  protected:
    
    //! construct_bisector_quadratic function construct a quadratic curve
    // to represent a bisector of two Voronoi sites.
    /*! The function constructs the quadratic of the form:
      $rx^2 + sy^2 + txy + ux + vy + w = 0$.
      It uses the traits class function to construct the conic, and then 
      converts the types to the correct type required by the concept.
      \return An output iterator whose value type is CGAL::Object which contain
      pair<X_monotone_curve_2,Multiplicity> or Point_2.
    */
    template <class OutputIterator>
    OutputIterator construct_bisector_quadratic (
      const Rational& r,
      const Rational& s,
      const Rational& t,
      const Rational& u,
      const Rational& v,
      const Rational& w,
      OutputIterator o) const
      {
        typedef std::list<X_monotone_curve_2> Xcurve_list;

        std::list<X_monotone_curve_2> x_conix;
        m_traits->construct_conic(r, s, t, u, v, w, 
                                  std::back_inserter(x_conix));

        convert_to_bisector(x_conix.begin(), x_conix.end(), o);

        return o;
      }

    /*! Converts a sequence of X_monotone_curve_2 to a sequence of CGAL::Object
      which contain pair<X_monotone_curve_2,Multiplicity>.
      The multiplicity used is always 1 (this is the multiplicity of the 
      bisector in this function).
      \return An output iterator containing the sequence.
    */
    template <class InputIterator, class OutputIterator>
    OutputIterator convert_to_bisector (InputIterator begin, InputIterator end,
                                         OutputIterator o) const
      {
        while (begin != end)
        {
          *o++ = make_object(std::make_pair(*begin++,
                                            Multiplicity(0)));
        }
        
        return o;
      }

  };

  Construct_projected_intersections_2 
  construct_projected_intersections_2_object()
    {
      return Construct_projected_intersections_2(this);
    }

  
  class Construct_projected_boundary_2
  {
  protected:
    Self *m_traits;
    
  public:
    
    Construct_projected_boundary_2(Self *traits)
      : m_traits(traits)
      {}
    
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s,
                              OutputIterator o) const
      {
        // Anisotropic sites are defined over the whole plane.
        return o;
      }
  };
  
  Construct_projected_boundary_2 
  construct_projected_boundary_2_object()
    {
      return Construct_projected_boundary_2(this);
    }
  
  
  
protected:
  //! construct_conic
  /*! The function constructs X_monotone_curve_2 objects of CKvA that 
    repesent the conic: r*x^2 + s*y^2 + t*x*y + u*x + v*y + w = 0
    The output is written to given output iterator.
    \return The updated output iterator - value_type of X_monotone_curve_2
  */
  template <class OutputIterator>
  OutputIterator construct_conic (
    const Rational& r,
    const Rational& s,
    const Rational& t,
    const Rational& u,
    const Rational& v,
    const Rational& w,
    OutputIterator o) const
    {
      typedef CGAL::Polynomial< Coefficient >           Polynomial_1;
      typedef CGAL::Polynomial< Polynomial_1 >          Polynomial_2;
      typedef CGAL::Polynomial_traits_d< Polynomial_1 > Polynomial_traits_1;
      typedef CGAL::Polynomial_traits_d< Polynomial_2 > Polynomial_traits_2;
        
      typedef std::list<CGAL::Object> Object_list;
      typedef CGAL::Exponent_vector EV;
      typedef std::pair<EV, Coefficient> Monomial;
      typedef CGAL::Fraction_traits< Rational > FT;

      typename FT::Decompose decompose;

      // getting integer from the rational coefficients.
      Coefficient rn, rd;
      Coefficient sn, sd;
      Coefficient tn, td;
      Coefficient un, ud;
      Coefficient vn, vd;
      Coefficient wn, wd;
      decompose(r, rn, rd);
      decompose(s, sn, sd);
      decompose(t, tn, td);
      decompose(u, un, ud);
      decompose(v, vn, vd);
      decompose(w, wn, wd);
        
      Coefficient d = rd*sd*td*ud*vd*wd;
      // need to multiply each numenator with the multiplication of the
      // denominator (except the current denom).
      rn *= sd*td*ud*vd*wd;
      sn *= rd*td*ud*vd*wd;
      tn *= rd*sd*ud*vd*wd;
      un *= rd*sd*td*vd*wd;
      vn *= rd*sd*td*ud*wd;
      wn *= rd*sd*td*ud*vd;
        
      std::vector<Monomial> monomials;
        
      monomials.push_back(Monomial(EV(2, 0), rn));
      monomials.push_back(Monomial(EV(0, 2), sn));
      monomials.push_back(Monomial(EV(1, 1), tn));
      monomials.push_back(Monomial(EV(1, 0), un));
      monomials.push_back(Monomial(EV(0, 1), vn));
      monomials.push_back(Monomial(EV(0, 0), wn));
        
      typename Polynomial_traits_2::Construct_polynomial 
        construct_polynomial;
      Polynomial_2 result = construct_polynomial(monomials.begin(),
                                                 monomials.end());
        
      Object_list x_objs;
        
      Curve_2 cur = 
        Base::instance().kernel().construct_curve_2_object()(result);
        
      this->make_x_monotone_2_object()
        (cur, std::back_inserter(x_objs));
        
        
      for (Object_list::iterator it = x_objs.begin();
           it != x_objs.end(); ++it)
      {
        const X_monotone_curve_2 *xcurve;
        if ((xcurve = object_cast<X_monotone_curve_2> (&(*it))) != NULL)
        {
          *o++ = *xcurve;

        }
        else
        {
          CGAL_error_msg("This function supports only X_monotone_curve_2");
        }

      }

      return o;
    }
};

// template <class T>
// std::ostream& 
// operator<< (std::ostream& os, 
//             const CGAL::_Anisotropic_site<T>& site)
// {
//   os << "Center: " << site.center() << std::endl;
//   os << "            / " << site.A() << " " << site.B() << " \\" << std::endl;
//   os << "Matrix: M = \\ " << site.B() << " " << site.C() << " /" << std::endl;
//   os << "R: " << site.D() << std::endl;
//   return os;
// }

} //namespace CGAL

#endif // CGAL_ANISOTROPIC_VORONOI_TRAITS_2_H
