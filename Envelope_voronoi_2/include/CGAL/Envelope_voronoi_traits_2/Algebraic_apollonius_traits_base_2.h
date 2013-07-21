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
// Author(s): Ophir Setter          <ophir.setter@cs.tau.ac.il>
//            Eric Berberich        <eric@mpi-inf.mpg.de>
//

/*! \file
 * This is the file contains the base class for (currently) 
 * Algebraic_apollonius_traits_2 and Algebraic_furthest_apollonius_traits_2.
 * These two traits classes have a lot in common - this base class contains
 * the similar things in both classes.
 */

#ifndef CGAL_ALGEBRAIC_APOLLONIUS_TRAITS_BASE_2_H
#define CGAL_ALGEBRAIC_APOLLONIUS_TRAITS_BASE_2_H

#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/compare/tribool.hpp>


#include <CGAL/basic.h>
#include <CGAL/Fraction_traits.h>

/*! \todo Move the function to a more "algebraic" place and remove this
  include. For the time being it is better to have a strange include then
  a code duplication. */
#include <CGAL/Voronoi_diagram_of_lines_3/vol_tools.h>
#include <CGAL/Envelope_voronoi_2/envelope_voronoi_assertions.h>

namespace CGAL {

namespace Envelope_voronoi_traits_2
{
  
/*!
 * \class _Apollonius_disk a rational circle is used as a site.
 * \todo This class only refers to the apollonius diagram, and this base class
 *       is planned to be the based to more diagrams. Move _Apollonius_disk to
 *       a different location and change the name of this file and class so
 *       it will say what it really is (base class for diagrams whose bisectors
 *       are composed of algebraic subcurves that are based on the CKvA.)
 */
template <typename Rational>
class _Apollonius_disk
{
public:
  typedef std::pair< Rational, Rational> Rat_point_2;
  
  friend bool operator< (const _Apollonius_disk& s1,
                         const _Apollonius_disk& s2)
  {
    // we can use <, > operators, but it is better to use compare since
    // it can save one of the comparisons (we usually need both < and >).
    CGAL::Comparison_result res = CGAL::compare(s1.center().first, 
                                                s2.center().first);
    if (res == CGAL::SMALLER)
      return true;
    if (res == CGAL::LARGER)
      return false;

    res = CGAL::compare(s1.center().second, s2.center().second);
    if (res == CGAL::SMALLER)
      return true;
    if (res == CGAL::LARGER)
      return false;

    res = CGAL::compare(s1.r(), s2.r());
    if (res == CGAL::SMALLER)
      return true;
    if (res == CGAL::LARGER)
      return false;

    return false;
  }

  // isn't this the default? Maybe I missed something in the standard.
  friend bool operator== (const _Apollonius_disk& s1,
                          const _Apollonius_disk& s2)
    {
    return (s1.center() == s2.center() &&
            s1.r() == s2.r());
  }
  
  friend std::ostream & operator<<(std::ostream & os,
                                   const _Apollonius_disk& site)
  {
    return os << "(" << site.center().first << ", " <<
      site.center().second << "), " << " " << site.r();
  }
    
  friend std::istream & operator>>(std::istream & is,
                                   _Apollonius_disk& site)
  {
    Rational x, y;
    is >> x >> y >> site._r;
    CGAL_envelope_voronoi_assertion(site._r >= 0);
    site._center = Rat_point_2(x, y);
    return is;
  }
    
protected:
  Rat_point_2 _center;

  // There is a way to support circles whose SQUARED radius is rational and 
  // not their radius. Somehow, I have the fealing it will be a lot slower.
  Rational _r;
    
public:
  _Apollonius_disk()
  {}
    
  _Apollonius_disk(const Rat_point_2& center, const Rational &r)
  : _center(center), _r(r)
  {}
        
  const Rat_point_2& center() const
  {
    return _center;
  }
        
  const Rational& r() const
  {
    return _r;
  }
};
}


/*!
 * \class Algebraic_apollonius_traits_base_2 The traits class
 * Distance_ret_req_ is algebraic structure tag that determines the 
 * requirement from the return type of the distance (in case the input 
 * arguments are rational.
 * \todo maybe pass the distance function just as functor (you loose state).
 */
template <class CurvedKernel_, class Derived_, 
          class Site_2_, class Distance_ret_req_>
class Algebraic_apollonius_traits_base_2
: public CurvedKernel_
{
public:
  typedef CurvedKernel_                                  Curved_kernel_2;
  typedef Derived_                                       Derived;
  typedef Site_2_                                        Site_2;
  typedef Distance_ret_req_                              Distance_ret_req;

  typedef Algebraic_apollonius_traits_base_2
  <Curved_kernel_2, Derived, Site_2, Distance_ret_req>   Self;
  typedef Curved_kernel_2                                Base;

  typedef typename Curved_kernel_2::Curve_kernel_2       Algebraic_kernel_2;
  typedef typename Algebraic_kernel_2::Coefficient       Coefficient;
  typedef typename CGAL::
    Get_arithmetic_kernel<Coefficient>::
    Arithmetic_kernel::Field_with_sqrt                   Field_with_sqrt;
  typedef typename Algebraic_kernel_2::Coordinate_1      Coordinate_1;

  typedef typename Coordinate_1::Rational                Rational;
  typedef typename Algebraic_kernel_2::
    Algebraic_kernel_1::Solve_1                          Solve_1;

  typedef typename Base::Point_2                         Point_2;
  typedef typename Base::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Base::Curve_2                         Curve_2;
  typedef typename Base::Multiplicity                    Multiplicity;
  typedef typename Base::Has_boundary_category           Has_boundary_category;

  typedef std::pair<Rational, Rational>                  Rat_point_2;

  typedef boost::numeric::interval<Field_with_sqrt>      Interval;
  
  typedef CGAL::Polynomial< Coefficient > Polynomial_1;
  typedef CGAL::Polynomial< Polynomial_1 > Polynomial_2;
  typedef CGAL::Polynomial_traits_d< Polynomial_1 > Polynomial_traits_1;
  typedef CGAL::Polynomial_traits_d< Polynomial_2 > Polynomial_traits_2;


  // we determine the number type to use for the distance function 
  // in case of a rational point according to the req. of the distance 
  // function that is specified by the type Distance_ret_req.
  typedef typename boost::mpl::if_< 
    boost::is_base_of<Field_with_sqrt_tag, Distance_ret_req>,
    Field_with_sqrt,
    Rational
    >::type                                              Distance_nt;
  
  

  // Field_with_sqrt is CORE::Expr !!!!
  static inline void convert_to_CORE(const Point_2 &p,
                                     Field_with_sqrt &x,
                                     Field_with_sqrt &y)
  {
    Coordinate_1 point_x = p.x();
  
    typedef typename CGAL::Algebraic_structure_traits< 
      Field_with_sqrt >::Root_of Root_of;
    Solve_1 solve;
    
    std::vector< Coordinate_1 > real_roots;
    solve(point_x.polynomial(), std::back_inserter(real_roots), true);
    int kid = std::distance(real_roots.begin(),
                            std::find(real_roots.begin(),
                                      real_roots.end(),
                                      point_x));
    CGAL_envelope_voronoi_assertion(real_roots[kid] == point_x);

    Root_of root_of;
    x = root_of(kid+1, point_x.polynomial().begin(),
                point_x.polynomial().end());
    CGAL_envelope_voronoi_assertion(point_x.compare(x) == CGAL::EQUAL);

    std::vector< Field_with_sqrt > roots;
    
    typedef CGAL::Polynomial< Field_with_sqrt > Poly_sqrt_1;
    typedef CGAL::Polynomial_traits_d< Poly_sqrt_1 > PT_sqrt_1;


    std::vector<Poly_sqrt_1> polys;
    typename PT_sqrt_1::Construct_polynomial construct_poly;
    polys.push_back(construct_poly(x));
    polys.push_back(construct_poly(Field_with_sqrt(0), Field_with_sqrt(1)));

    typename Polynomial_traits_2::Substitute substitue;
    Poly_sqrt_1 res_poly = substitue(p.curve().polynomial_2(),
                                     polys.begin(), polys.end());
    solve_quadratic< Field_with_sqrt, Field_with_sqrt >
      (res_poly, std::back_inserter(roots));
      
    y = roots[p.arcno()];
  }

protected:

  /*! \todo Move this function to somewhere more "algebraic". */
  template < class Field_, class FieldWithSqrt, class NT_, 
    class OutputIterator >
    static inline int solve_quadratic(Polynomial< NT_ > p, OutputIterator it) 
    {
      
      typedef Field_                                         Field;
      typedef FieldWithSqrt                                  Field_with_sqrt;
      typedef NT_                                            NT;
      typedef OutputIterator                                 Output_iterator;
    
      switch (p.degree()) 
      {
      case 0:
        CGAL_envelope_voronoi_precondition(!p.is_zero());
        return 0;
      case 1: 
      {
        Field p1(p[1]), p0(p[0]);
        *it = (- p0 / p1);
        return 1;
      }
      case 2: 
      {
        NT discri = p[1]*p[1] - NT(4)*p[0]*p[2];
        switch (sign(discri))
        {
        case NEGATIVE:
          // no roots
          return 0;
        case ZERO: 
        {
          // one root
          Field p1(p[1]), p2(p[2]);
          *it = ((- p1) / (Field(2)*p2));
          return 1;
        }
        case ::CGAL::POSITIVE:
          // two roots; compute them _in_ascending_order_
          Field p0(p[0]), p1(p[1]), p2(p[2]), discr(discri);
            
          Field d = Field(2) * p2;
          Field_with_sqrt sqrt_discr = sqrt(Field_with_sqrt(discr));
          if (d > Field(0)) {
            *it   = ((-p1 - sqrt_discr) / d);
            *++it = ((-p1 + sqrt_discr) / d);
          } else {
            *it   = ((-p1 + sqrt_discr) / d);
            *++it = ((-p1 - sqrt_discr) / d);
          }
          return 2;
        }}
      default:
        CGAL_error_msg("solve_quadratic(p) called with p.degree() != 0, 1, 2");
        return -1; // never reached
      }
    }

public:
  typedef Site_2                                 Xy_monotone_surface_3;
  typedef Site_2                                 Surface_3;
    
  class Make_xy_monotone_3
  {
  public:
        
    template <class OutputIterator>
      OutputIterator operator()(const Surface_3& s,
                                bool is_lower,
                                OutputIterator o) 
    {
      *o++ = s;
      return o;
    }
  };
    
  Make_xy_monotone_3 make_xy_monotone_3_object()
  {
    return Make_xy_monotone_3();
  }
    
  class Compare_z_at_xy_3
  {
  protected:
    const Base *m_base;
      
  public:
        
  Compare_z_at_xy_3(const Base* base)
    : m_base(base) {}

    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2)
    {
      Rational xl = m_base->kernel().lower_boundary_x_2_object() (p.xy());
      Rational xu = m_base->kernel().upper_boundary_x_2_object() (p.xy());
        
      Rational yl = m_base->kernel().lower_boundary_y_2_object() (p.xy());
      Rational yu = m_base->kernel().upper_boundary_y_2_object() (p.xy());
        
      Interval xi(xl, xu);
      Interval yi(yl, yu);

      // I have no idea why this try-catch is needed, but for some reason
      // on my gcc 4.2.3 an "boost::interval: uncertain comparison" exception
      // is thrown even though we picked the tribool comparison.
      // My only reaction is WTF?!?!
      try
      {
        Interval d1 = Derived::distance(xi, yi, h1);
        Interval d2 = Derived::distance(xi, yi, h2);
        
        using namespace boost::numeric::interval_lib::compare::tribool;
        std::cout << "comparing using tribool" << std::endl;
        boost::tribool res = (d1 < d2);
        std::cout << "END comparing using tribool" << std::endl;
        
        if (res)
          return SMALLER;
        if (!res)
          return LARGER;
      }
      catch(const boost::numeric::interval_lib::comparison_error &)
      {}

      Field_with_sqrt px, py;
      convert_to_CORE(p, px, py);
      Field_with_sqrt de1 = Derived::distance(px, py, h1);
      Field_with_sqrt de2 = Derived::distance(px, py, h2);
      return CGAL::compare(de1, de2);
    }
      
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2)
    {
      // Take a representative from the interior of the point 
      // compare surfaces at that point 
      Point_2 p = m_base->construct_interior_vertex_2_object() (cv);
      return (*this)(p, h1, h2);
    }
      
      
    Comparison_result operator()(const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2)
        
    {
      Distance_nt de1 = Derived::distance(Distance_nt(0), 
                                          Distance_nt(0), h1);
      Distance_nt de2 = Derived::distance(Distance_nt(0), 
                                          Distance_nt(0), h2);
      
      return CGAL::compare(de1, de2);
    }
      
  };
  
  Compare_z_at_xy_3 compare_z_at_xy_3_object()
  {
    return Compare_z_at_xy_3(this);
  }

  class Compare_z_at_xy_above_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    {
      std::pair<Rational, Rational> p_above = 
        Voronoi_diagram_of_lines_3::
        rational_point_above_or_below_arc<Base>(cv, true);

      Distance_nt de1 = Derived::template distance<Distance_nt>(
        p_above.first, p_above.second, h1);
      Distance_nt  de2 = Derived::template distance<Distance_nt>(
        p_above.first, p_above.second, h2);

      return CGAL::compare(de1, de2);
    }
      
  };

  Compare_z_at_xy_above_3 compare_z_at_xy_above_3_object() const
  {
    return Compare_z_at_xy_above_3();
  }

  class Compare_z_at_xy_below_3
  {
  protected: 
    const Self *m_self;
  public:

  Compare_z_at_xy_below_3(const Self *self)
    : m_self(self) {}
      
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2)
    {
      Comparison_result res = CGAL::opposite(m_self->
                                             compare_z_at_xy_above_3_object() 
                                             (cv, h1, h2));
#ifndef CGAL_NDEBUG 
      std::pair<Rational, Rational> p_below = 
        Voronoi_diagram_of_lines_3::
        rational_point_above_or_below_arc<Base>(cv, false);

      Distance_nt de1 = Derived::template distance<Distance_nt>(
        p_below.first, p_below.second, h1);
      Distance_nt de2 = Derived::template distance<Distance_nt>(
        p_below.first, p_below.second, h2);

      assert(res == CGAL::compare(de1, de2));
#endif
      return res;
    }

  };

  Compare_z_at_xy_below_3 compare_z_at_xy_below_3_object()
  {
    return Compare_z_at_xy_below_3(this);
  }


  class Construct_projected_boundary_2
  {
  public:

    template <class OutputIterator>
      OutputIterator operator()(const Xy_monotone_surface_3& s,
                                OutputIterator o) const
    {
      return o;
    }
  };

  Construct_projected_boundary_2 
    construct_projected_boundary_2_object()
  {
    return Construct_projected_boundary_2();
  }

  template <class OutputIterator>
    OutputIterator apollonius_bisector(const Rational &xi,
                                       const Rational &yi,
                                       const Rational &ri,
                                       const Rational &xj,
                                       const Rational &yj,
                                       const Rational &rj,
                                       Self & traits,
                                       OutputIterator o) const
  {
    typedef CGAL::Exponent_vector EV;
    typedef std::pair<EV, Coefficient> Monomial;
    typedef CGAL::Fraction_traits< Rational > FT;
    typedef std::list<CGAL::Object> Object_list;

    typename FT::Decompose decompose;

    Comparison_result r_res = CGAL::compare(ri, rj);
    if (r_res == EQUAL)
    {
      Rational dx = xj - xi;
      Rational dy = yj - yi;
      Rational u = 2*dx;
      Rational v = 2*dy;
      Rational w = -dx*(xi + xj) -dy*(yi + yj);
      
      std::vector<Monomial> monomials;

      // convert to integers
      Coefficient un, ud;
      Coefficient vn, vd;
      Coefficient wn, wd;
      decompose(u, un, ud);
      decompose(v, vn, vd);
      decompose(w, wn, wd);

      Coefficient d = ud*vd*wd;
      // need to multiply each numenator with the multiplication of the
      // denominator (except the current denom).
      un *= vd*wd;
      vn *= ud*wd;
      wn *= ud*vd;

      monomials.push_back(Monomial(EV(1, 0), un));
      monomials.push_back(Monomial(EV(0, 1), vn));
      monomials.push_back(Monomial(EV(0, 0), wn));

      typename Polynomial_traits_2::Construct_polynomial 
        construct_polynomial;
      Polynomial_2 result = construct_polynomial(monomials.begin(),
                                                 monomials.end());

      Object_list x_objs;
      std::list<X_monotone_curve_2> x_conix;
      
      Curve_2 cur = Base::instance().kernel().
        construct_curve_2_object()(result);
      
      traits.make_x_monotone_2_object()
        (cur, std::back_inserter(x_objs));
      
      for (Object_list::iterator it = x_objs.begin();
           it != x_objs.end(); ++it)
      {
        const X_monotone_curve_2 *xcurve;
        if ((xcurve = object_cast<X_monotone_curve_2> (&(*it))) != NULL)
        {
          *o++ = make_object(std::make_pair(*xcurve, Multiplicity(1)));

        }
      }

      return o;
      
    }

    /* The following maple code produces the bisector of the two circles. 
       We decided to keep to old, hand-generated code.
             
       with(LinearAlgebra); with(codegen);
       p1 := Vector([X1, Y1])
       p2 := Vector([X2, Y2])
       X := Vector([x, y])
       XX1 := DotProduct(X - p1, X - p1, conjugate = false)
       XX2 := DotProduct(X - p2, X - p2, conjugate = false)
       /                                      2            \
       Q := expand\(XX1 + XX2 + (-1) (R2 - R1) (R2 - R1))  - 4 XX1 XX2/
       a := coeff(coeff(Q, x, 2), y, 0)
       b := coeff(coeff(Q, x, 0), y, 2)
       c := coeff(coeff(Q, x, 1), y, 1)
       d := coeff(coeff(Q, x, 1), y, 0)
       e := coeff(coeff(Q, x, 0), y, 1)
       e := coeff(coeff(Q, x, 0), y, 0)
       A := array([a, b, c, d, e])
       C(A, optimized)
    */

    // if one site is completly inside the other site return empty
    // object.
    const Rational xmx = xi - xj;
    const Rational ymy = yi - yj;
    const Rational rmr = ri - rj;
      
    const Rational xx_2 = xmx * xmx;
    const Rational yy_2 = ymy * ymy;
    const Rational RS = rmr * rmr;
          
    CGAL_envelope_voronoi_assertion_msg (xx_2 + yy_2 != RS,                \
                        "temp. not supporting singular points");

    // temp temp
    if (xx_2 + yy_2 == RS)
      return o;

    if (xx_2 + yy_2 < RS)
      return o;

    // The bisector is a hyperbolic arc with:
    //  r = 4 [(xi-xj)^2 - RS]
    //  s = same as r (exchange x and y)
    //  t = 8 (xi - xj) (yi - yj)
    //  u = 4 [-(yi + yj)(yi - yj)(xi - xj) - (xi + xj)(xi - xj)^2 + 
    //      RS(xi + xj)]
    //  v = same as u (exchange x and y)
    //  w := (xi^2 - xj^2)^2 + (yi^2 - yj^2)^2 + (xi^2 - yj^2)^2 +
    //	     (yi^2 - xj^2)^2 - (xi^2 - yi^2)^2 - (xj^2 - yj^2)^2
    // 	      -2RS (xi^2 + yi^2 + xj^2 + yj^2) + RS^2
    // where:
    //  RS = (ri-rj)^2
    const Rational r = 4 * (xx_2 - RS);
    const Rational s = 4 * (yy_2 - RS);

    const Rational t = 8 * xmx * ymy;

    const Rational xpx = xi + xj;
    const Rational ypy = yi + yj;

    const Rational u = 4 * (RS*xpx - ypy*ymy*xmx - xpx*xmx*xmx);
    const Rational v = 4 * (RS*ypy - xpx*xmx*ymy - ypy*ymy*ymy);

    const Rational xi_2 = xi * xi;
    const Rational xj_2 = xj * xj;
    const Rational yi_2 = yi * yi;
    const Rational yj_2 = yj * yj;

    const Rational xi2xj2 = (xi_2 - xj_2);
    const Rational yi2yj2 = (yi_2 - yj_2);
    const Rational xi2yj2 = (xi_2 - yj_2);
    const Rational yi2xj2 = (yi_2 - xj_2);
    const Rational xi2yi2 = (xi_2 - yi_2);
    const Rational xj2yj2 = (xj_2 - yj_2);
    const Rational w = xi2xj2 * xi2xj2 + yi2yj2 * yi2yj2
      + xi2yj2 * xi2yj2 + yi2xj2 * yi2xj2
      - xi2yi2 * xi2yi2 - xj2yj2 * xj2yj2
      - 2 * RS * (xi_2 + yi_2 + xj_2 + yj_2)
      + RS*RS;

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
    std::list<X_monotone_curve_2> x_conix;

    Curve_2 cur = Base::instance().kernel().construct_curve_2_object()(result);

    traits.make_x_monotone_2_object()
      (cur, std::back_inserter(x_objs));

    for (Object_list::iterator it = x_objs.begin();
         it != x_objs.end(); ++it)
    {
      const X_monotone_curve_2 *xcurve;
      if ((xcurve = object_cast<X_monotone_curve_2> (&(*it))) != NULL)
      {
        x_conix.push_back(*xcurve);
      }
      else
      {
        CGAL_error_msg("The bisector at this point "
                       "should be a hyperbola");
      }
    }
          

    CGAL_envelope_voronoi_assertion(x_conix.size() == 2 || x_conix.size() == 4);

    if (x_conix.size() == 2)
    {
      // we have only two x-mono curves. we take the correct branch
      // according to the sites (on the side of the smaller site)
      Comparison_result res;

      // if both arcs have arcno 0 then we have a vertical asymptot
      if (x_conix.front().arcno() == x_conix.back().arcno())
      {
        // we need to compare the x-coords
        Comparison_result res;
        if (r_res == SMALLER) // ri < rj
        {
          res = CGAL::compare(xi, xj);
        }
        else
        {
          res = CGAL::compare(xj, xi);
        }
        CGAL_envelope_voronoi_assertion(res != EQUAL);

        if (res == SMALLER)
          *o++ = make_object(std::make_pair(x_conix.front(),
                                            Multiplicity(1)));
        else
          *o++ = make_object(std::make_pair(x_conix.back(),
                                            Multiplicity(1)));
      }
      else
      {
        if (r_res == SMALLER) // ri < rj
        {
          res = CGAL::compare(yi, yj);
        }
        else
        {
          res = CGAL::compare(yj, yi);
        }
        CGAL_envelope_voronoi_assertion(res != EQUAL);
              
        // we don't have vert asymptot
        if (res == SMALLER && x_conix.front().arcno() == 0)
          *o++ = make_object(std::make_pair(x_conix.front(),
                                            Multiplicity(1)));
        else
          *o++ = make_object(std::make_pair(x_conix.back(),
                                            Multiplicity(1)));
      }
    }
    else // x_conix.size() == 4
    {
      CGAL_envelope_voronoi_assertion(x_conix.size() == 4);

      Comparison_result res;
      if (r_res == SMALLER) // ri < rj
      {
        res = CGAL::compare(xi, xj);
      }
      else
      {
        res = CGAL::compare(xj, xi);
      }
      CGAL_envelope_voronoi_assertion(res != EQUAL);

      // we can assume that the order of the curves from make_x_monotone
      // is lexicographic
      typename std::list<X_monotone_curve_2>::iterator it =
        x_conix.begin();
      if (res == SMALLER)
      {
        // we take the first two x-mono curves because we need the left
        // subcurves.
        *o++ = make_object(std::make_pair(*it, Multiplicity(1)));
        ++it;
        *o++ = make_object(std::make_pair(*it, Multiplicity(1)));
      }
      else
      {
        // we take the right x-mono curves.
        ++it; ++it;
        *o++ = make_object(std::make_pair(*it, Multiplicity(1)));
        ++it;
        *o++ = make_object(std::make_pair(*it, Multiplicity(1)));
      }
    }

    return o;
  }
};

} //namespace CGAL

#endif
