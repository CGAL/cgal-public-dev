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
// $Id:  $
// 
//
// Author(s): Ophir Setter          <ophir.setter@post.tau.ac.il>
//

/*! \file Moebius_diagram_traits_2.h
  This is the file containing a traits class to create moebius voronoi diagrams
  using the lower envelope code.
 */

#ifndef CGAL_MOEBIUS_DIAGRAM_TRAITS_2_H
#define CGAL_MOEBIUS_DIAGRAM_TRAITS_2_H


#include <CGAL/basic.h>
#include <CGAL/Sqrt_extension.h>

#include <CGAL/Envelope_3/Envelope_base.h>
#include <CGAL/Arr_circle_linear_traits_2.h>

#include <CGAL/Envelope_voronoi_2/envelope_voronoi_assertions.h>

namespace CGAL {

/*!
 * \class Moebius_diagram_traits_2 Traits class for computing Mobius diagram
 */
template <class Kernel_, bool Filter_ = true>
class Moebius_diagram_traits_2 
: public Arr_circle_linear_traits_2< Kernel_, Filter_ >
{
public:
  typedef Moebius_diagram_traits_2< 
    Kernel_, Filter_ >                                Self;
  typedef Arr_circle_linear_traits_2< 
    Kernel_, Filter_ >                                Base;
  typedef typename Base::Point_2                      Point_2;
  typedef typename Point_2::NT                        NT;
  typedef typename Point_2::CoordNT                   CoordNT;

  typedef Kernel_                                     Kernel;
  typedef typename Kernel::FT                         FT;
  typedef typename Kernel::Point_2                    Kernel_point_2;
  typedef typename Kernel::Line_2                     Line_2;
  typedef typename Kernel::Circle_2                   Circle_2;

  typedef typename Base::Curve_2                      Curve_2;
  typedef typename Base::X_monotone_curve_2           X_monotone_curve_2;
  
  /*!
   * \class The class represents an anistropic voronoi diagram site.
   *        The distance function from an anistropic site is:
   *        
   *        $d(x) = M * (x - p) - R$
   *
   *        Where M, R are scalars.
   *
   */  
  class _Moebius_site
  {
  public:
  protected:
    Kernel_point_2 _point;
    FT      _M;
    FT      _R;
        
  public:
    // default constructor
    _Moebius_site()
    {}

    _Moebius_site(const Kernel_point_2& point, const FT& M, const FT& R)
      : _point(point), _M(M), _R(R)
    {}
        
    const Kernel_point_2& point () const
    {
      return _point;
    }
        
    const FT& M() const
    {
      return _M;
    }

    const FT& R() const
    {
      return _R;
    }

    /*!
     * Computes the distance from the moebius site.
     * The distance from an anistropic diagram is defined as:
     * 
     * $d(x) = M * (x - p)^2 - R$
     *
     * \param p The point to calculate the distance from.
     */
    CoordNT distance(const Point_2& p) const
    {
      CGAL_envelope_voronoi_precondition(p.x().gamma() == p.y().gamma() || \
                                         p.x().gamma() == 0 ||          \
                                         p.y().gamma() == 0);

      // \todo move the Sqrt_extension that will make this a lot easier.
      // we use the fact that the center of the mobius is on a rational point.
      typedef CGAL::Sqrt_extension<NT, NT>               Sqrt_ext;
      

      Sqrt_ext x(p.x().alpha() - point().x(), p.x().beta(), p.x().gamma());
      Sqrt_ext y(p.y().alpha() - point().y(), p.y().beta(), p.y().gamma());

      Sqrt_ext x_2 = CGAL::square(x);
      Sqrt_ext y_2 = CGAL::square(y);
      
      // just to make sure
      Sqrt_ext squared_dist = x_2 + y_2;
      Sqrt_ext res = M() * squared_dist - R();

      return CoordNT(res.a0(), res.a1(), res.root());
    }
    

    friend std::ostream& 
    operator<< (std::ostream& os, 
                _Moebius_site& site)
      {
        os << "Point: " << site.point() << std::endl;
        os << "M: " << site.M() << "\tR: " << site.R() << std::endl;
        return os;
      }
    
    friend std::istream& 
    operator>> (std::istream& is, 
                _Moebius_site& site)
      {
        Kernel_point_2 p;
        FT M, R;
        
        is >> p >> M >> R;
        site = _Moebius_site(p, M, R);
        
        return is;
      }

  };

     
  typedef _Moebius_site                                 Site_2;
  
  /* construct the bisector between two points */
  class Construct_bisector_2
  {
    typedef typename Base::Make_x_monotone_2            Make_x_monotone_2;
  protected:
    const Make_x_monotone_2 _make_x;
  public:
    Construct_bisector_2()
    {}

    template <class OutputIterator>
      OutputIterator operator()(const Site_2& s1,
                                const Site_2& s2,
                                OutputIterator o) const
    {
      typedef std::list<Object> X_curve_list;

      Kernel ker;

      typename Kernel::Compute_x_2 compute_x = ker.compute_x_2_object();
      FT x1 = compute_x(s1.point());
      FT x2 = compute_x(s2.point());
      typename Kernel::Compute_y_2 compute_y = ker.compute_y_2_object();
      FT y1 = compute_y(s1.point());
      FT y2 = compute_y(s2.point());


      FT A = s1.M() - s2.M();
      if (sign(A) == ZERO)
      {
        if (s1.point() == s2.point())
          return o;
        *o++ = X_monotone_curve_2(CGAL::bisector(s1.point(), s2.point()));
        return o;
      }
      
      // we calculate the center and radius of the bisecting circle.
      FT D = -2 * (x1 - x2);
      FT E = -2 * (y1 - y2);
      FT F = (x1*x1 + y1*y1) - (x2*x2 + y2*y2) - (s1.R() - s2.R());

      FT x = -D / 2*A;
      FT y = -E / 2*A;
      FT r_2 = (-F + (D*D + E*E) / (4*A*A)) / A;
      if (sign(r_2) == NEGATIVE)
        return o;

      Kernel_point_2 p (x, y);
      Circle_2 circle = ker.construct_circle_2_object()(p, r_2);
      Curve_2 curve(circle);
      X_curve_list l;
      (*const_cast<typename Base::Make_x_monotone_2*>(&_make_x))
        (curve, std::back_inserter(l));
      
      for (typename X_curve_list::iterator it = l.begin();
           it != l.end(); ++it)
      {
        X_monotone_curve_2 xcurve;
        if (assign(xcurve, *it))
        {
          *o++ = xcurve;
        }
      }

      return o;
    }
  };

  Construct_bisector_2 construct_bisector_2_object() const
  {
    return Construct_bisector_2();
  }

  /* compares distance above x-monotone curve */
  class Compare_distance_above_2
  {
  public:
    Comparison_result operator()(const Site_2& h1,
                                 const Site_2& h2, 
                                 const X_monotone_curve_2& cv) const
      
    {
      if (cv.orientation() == COLLINEAR)
      {
        Kernel ker;
        Comparison_result res = ker.compare_y_2_object() 
          (h2.point(), h1.point());
        if (res == EQUAL)
          return ker.compare_x_2_object() (h2.point(), h1.point());
        return res;
      }
      else
      {
        Comparison_result res = CGAL::compare(h1.M(), h2.M());
        bool right = cv.is_directed_right();
        bool orient = (cv.orientation() == COUNTERCLOCKWISE);
        if (right == orient)
          return opposite(res);
        else
          return res;
      }
    }
  };

  Compare_distance_above_2 compare_distance_above_2_object() const
  {
    return Compare_distance_above_2();
  }

  class Compare_distance_at_point_2
  {
  public:
    Comparison_result operator()(const Site_2& h1,
                                 const Site_2& h2,
                                 const Point_2& p) const
    {
      return CGAL::compare(h1.distance(p), h2.distance(p));
    }
  };
  
  Compare_distance_at_point_2 compare_distance_at_point_2_object() const
  {
    return Compare_distance_at_point_2();
  }

  class Construct_point_on_x_monotone_2
  {
  public:
    Point_2 operator() (const X_monotone_curve_2& xcurve) const
    {
      // construct isolating intervals of both end points of the xcurve
      // to get a rational in the middle. Then intersect the line in the 
      // middle with the arc to get the point.
      // If one of the vertices is really a rational and not a sqrt extension
      // it is easier.
      // Handle unbounded curves also + vertical lines.

      CGAL_error();
      return Point_2();
    }
  };
  
  Construct_point_on_x_monotone_2 construct_point_on_x_monotone_2_object() const
  {
    return Construct_point_on_x_monotone_2();
  }

  /* because one site can't dominate other site in this case, we can always 
     return equal */
  class Compare_dominance_2
  {
  public:
    Comparison_result operator()(const Site_2& h1,
                                 const Site_2& h2) const
    {
      return EQUAL;
      // return compare(h1.M(), h2.M());
    }
  };
  
  Compare_dominance_2 compare_dominance_2_object() const
  {
    return Compare_dominance_2();
  }
};

} //namespace CGAL

#endif
