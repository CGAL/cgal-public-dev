// Copyright (c) 2008  GeometryFactory Sarl (France).
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
// $URL: svn+ssh://rozapoga@scm.gforge.inria.fr/svn/cgal/trunk/GraphicsView/include/CGAL/Qt/PainterOstream.h $
// $Id: PainterOstream.h 46265 2008-10-14 10:43:45Z lrineau $
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_ARR_PAINTER_OSTREAM_H
#define CGAL_QT_ARR_PAINTER_OSTREAM_H

#include <QPainter>
#include <QPen>
#include <QRectF>
#include <QPainterPath>
#include <CGAL/Qt/Converter.h>
//#include <CGAL/Qt/ArrConverter.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Minkowski_sum_with_circle_2.h>

#define SHOW_ARCS // comment this definition to see full circles in UI

namespace CGAL {

  template <class CK>
  class  Circular_arc_point_2;

//  template <class CK>
//  class  Circular_arc_2;

  template <class CK>
  class  Line_arc_2;

namespace Qt {

template<typename ArrTraits>
struct ArrTraitsHelper
{
  typedef ArrTraits                                     Traits_2;
  typedef typename Traits_2::Kernel                     Kernel;
  typedef Converter<Kernel>                             ConverterType;

  typedef typename CGAL::Point_2<Kernel>        Point_2;
  typedef typename CGAL::Line_2<Kernel>         Line_2;
  typedef typename CGAL::Segment_2<Kernel>  Segment_2;
//  typedef typename CGAL::Segment_2<
//    typename CGAL::Cartesian<typename Traits_2::CoordNT> 
//    > Segment_2;
  typedef typename CGAL::Circle_2<Kernel>       Circle_2;
//  typedef typename CGAL::Circular_arc_2<Kernel>       Circular_arc_2;
  
  static bool is_linear(const typename Traits_2::X_monotone_curve_2& c) { return c.is_linear(); }
  static bool is_circular(const typename Traits_2::X_monotone_curve_2& c) { return c.is_circular(); }

  static Line_2 supporting_line(const typename Traits_2::X_monotone_curve_2& c)
  {
      // std::clog << "non-conic seg" << std::endl;
  
      return c.supporting_line();  
  }

  static Circle_2 supporting_circle(const typename Traits_2::X_monotone_curve_2& c)
  {
      // std::clog << "non-conic arc" << c.source() << c.target() << std::endl;
      return c.supporting_circle();
  }
  
//  static Circular_arc_2 circular_arc(const typename Traits_2::X_monotone_curve_2& c)
//  {
//      return return Circular_arc_2(supporting_circle(c), c.source(), c.target());
//  }
};


typedef Minkowski_sum_with_circle_2<CGAL::Polygon_2<
  CGAL::Cartesian<CGAL::Gmpq > 
  //CGAL::Cartesian<CORE::BigRat > 
    > >::Conic_traits_2 Conic_traits_2;

typedef CGAL::Cartesian<CGAL::CORE_algebraic_number_traits::Algebraic >
    Alg_kernel;

/*
Error C2764:

template<typename NT>
struct ArrTraitsHelper<typename Minkowski_sum_with_circle_2<typename CGAL::Polygon_2<
    typename CGAL::Cartesian<NT > 
    > >::Conic_traits_2 >
{
  typedef typename Minkowski_sum_with_circle_2<typename CGAL::Polygon_2<
    typename CGAL::Cartesian<NT > 
    > >::Conic_traits_2 Traits_2;
*/

template<>
struct ArrTraitsHelper<Conic_traits_2 >
{
  typedef Conic_traits_2    Traits_2;

  //  typedef Traits_2::Alg_kernel  Kernel;
  typedef Alg_kernel  Kernel;

  typedef Converter<Kernel>                             ConverterType;

  typedef CGAL::Point_2<Kernel>         Point_2;
  typedef CGAL::Line_2<Kernel>          Line_2;
  typedef CGAL::Segment_2<Kernel>       Segment_2;
  typedef CGAL::Circle_2<Kernel>        Circle_2;
  //typedef CGAL::Circular_arc_2<Kernel>       Circular_arc_2;
  
  // conic curve equasion:
  // r x^2 + s y^2 + t x y + u x + v y + w = 0
  static bool is_linear(const Traits_2::X_monotone_curve_2& c) 
  { 
    // conic line: r = s = t = 0
    return ((c.t() == 0) && (c.r() == 0) && (c.s() == 0)); 
  }

  static bool is_circular(const Traits_2::X_monotone_curve_2& c) 
  { 
    // conic circle: r = s, t = 0
    return ((c.t() == 0) && (c.r() == c.s())); 
  }

  static Segment_2 supporting_line(const Traits_2::X_monotone_curve_2& c)
  {
      // std::clog << "conic seg" << std::endl;
      return Segment_2(c.source(), c.target());
  }

  static Circle_2 supporting_circle(const Traits_2::X_monotone_curve_2& c)
  {
    // std::clog << "conic arc" << c.source() << c.target() << std::endl;

    // verify r = s
    CGAL_assertion(c.r() == c.s());
    
//    std::clog << c.r() << "*x^2 + " 
//              << c.s() << "*y^2 + " 
//              << c.t() << "*xy + " 
//              << c.u() << "*x + " 
//              << c.v() << "*y + " 
//              << c.w() << " = 0"
//              << std::endl;

    Kernel::FT s(c.s()), u(c.u()), v(c.v()), w(c.w());
    // circle (x - a)^2 + (y - b)^2 = r^2
    // center (a, b), radius r
    // a = - u / 2s, b = - v / 2s
    Kernel::FT minus_a = (u / (2 * s));
    Kernel::FT minus_b = (v / (2 * s));
    Point_2 center( - minus_a , - minus_b);
 //   std::clog << "center(" << center.x() << "," << center.y() << ")" << std::endl;

    Kernel::FT sqr_a = (minus_a * minus_a);
    Kernel::FT sqr_b = (minus_b * minus_b);

    // r^2 = a^2 + b^2 - w/s
    Kernel::FT sqr_radius = sqr_a + sqr_b - (w / s);
    
    //std::clog << "point_on_arc(" << c.source().x() << "," << c.source().y() << ")" << std::endl;
    //Kernel::FT diff_x = c.source().x() + minus_a;
    //Kernel::FT diff_y = c.source().y() + minus_b;
    //Kernel::FT sqr_radius = diff_x * diff_x + diff_y * diff_y;
    CGAL_assertion(sqr_radius >= 0);

//    std::clog << " < = > ? " << std::endl;
    //std::clog << "(x + " << minus_a << ")^2  + " 
    //          << "(y + " << minus_b << ")^2  = " 
    //          << sqr_radius 
    //          << std::endl;

    return Circle_2(center, sqr_radius, c.orientation());
  }
  
//  static Circular_arc_2 circular_arc(const Traits_2::X_monotone_curve_2& c)
//  {
//      return Circular_arc_2(supporting_circle(c), c.source(), c.target());
//  }

};


template <typename ArrTraits>
class ArrPainterOstream {

public:
     
  typedef ArrTraitsHelper<ArrTraits> TraitsHelper;

  typedef typename TraitsHelper::Kernel Kernel;
  typedef typename TraitsHelper::ConverterType ConverterType;

private:
  QPainter* qp;
  ConverterType convert;
  
public:
  ArrPainterOstream(QPainter* p, QRectF rect = QRectF())
    : qp(p), convert(rect)
  {}
/*
  ArrPainterOstream& operator<<(const Point_2<Kernel>& p)
  {
    qp->drawPoint(convert(p));
    return *this;
  }
  */
  
  ArrPainterOstream& operator<<(const Segment_2<Kernel>& s)
  {
    // std::clog << "Segment_2" << std::endl;  
    qp->drawLine(convert(s.source()), convert(s.target()));
    return *this;
  }
  
  
  ArrPainterOstream& operator<<(const Ray_2<Kernel>& r)
  {
    // std::clog << "Ray_2" << std::endl;  
    qp->drawLine(convert(r));
    return *this;
  }

  
  ArrPainterOstream& operator<<(const Line_2<Kernel>& l)
  {
    // std::clog << "Line_2" << std::endl; 
    qp->drawLine(convert(l));
    return *this;
  }

/*
  ArrPainterOstream& operator<<(const Triangle_2<Kernel>& t)
  {
    qp->drawPolygon(convert(t));
    return *this;
  }
*/
  ArrPainterOstream& operator<<(const Iso_rectangle_2<Kernel>& r)
  {
    // std::clog << "Iso_rectangle_2" << std::endl; 
    qp->drawRect(convert(r));
    return *this;
  }


  ArrPainterOstream& operator<<(const Circle_2<Kernel>& c)
  {
    // std::clog << "Circle_2" << std::endl;
    qp->drawEllipse(convert(c.bbox()));
    return *this;
  }
      
  ArrPainterOstream& operator<<(const typename ArrTraits::X_monotone_curve_2& c)
  {
    // std::clog << "X_monotone_curve_2" << std::endl;    
    if(TraitsHelper::is_linear(c))
    {
      // TODO: fix so only supporting_line() is needed
  //    std::clog << "line " << TraitsHelper::supporting_line(c) << std::endl;
  //    std::clog << "line " << c.source() << " --> " << c.target() 
  //              << std::endl;

//      (*this) << TraitsHelper::supporting_line(c);
      (*this) << Segment_2<Kernel>(
        Point_2<Kernel>(to_double(c.source().x()), to_double(c.source().y())),
        Point_2<Kernel>(to_double(c.target().x()), to_double(c.target().y()))
        );

    }
    else
    if(TraitsHelper::is_circular(c))
    {
    //  std::clog << "circle " << TraitsHelper::supporting_circle(c) 
    //            << std::endl;

#ifdef SHOW_ARCS
      drawArc(TraitsHelper::supporting_circle(c), c.source(), c.target());
#else
      (*this) << TraitsHelper::supporting_circle(c);
#endif

     }

    return *this;
  }


/*
  ArrPainterOstream& operator<<(const Circular_arc_point_2<Kernel>& p)
  {
    typedef typename Kernel::Point_2   Point_2;
    (*this) << Point_2(to_double(p.x()), to_double(p.y()));
    return *this;
  }
*/


  ArrPainterOstream& drawArc(const typename TraitsHelper::Circle_2& circ, const typename TraitsHelper::Traits_2::Point_2& source, const typename TraitsHelper::Traits_2::Point_2& target)
  {
    // std::clog << "drawArc" << std::endl;             
    const typename TraitsHelper::Point_2 & center = circ.center();

    double asource = std::atan2( -to_double(source.y() - center.y()),
                 to_double(source.x() - center.x())); 
    double atarget = std::atan2( -to_double(target.y() - center.y()),
                 to_double(target.x() - center.x()));

    std::swap(asource, atarget);
    double aspan = atarget - asource;

    if(aspan < 0.)
      aspan += 2 * CGAL_PI;

    const double coeff = 180*16/CGAL_PI;
    qp->drawArc(convert(circ.bbox()), 
        (int)(asource * coeff), 
             (int)(aspan * coeff));
    return *this;
  }

  ArrPainterOstream& operator<<(const Circular_arc_2<Kernel>& arc)
  {
    // std::clog << "Circular_arc_2" << std::endl;        
    const typename ArrTraits::Circle_2 & circ = arc.supporting_circle();
    const typename ArrTraits::Point_2 & center = circ.center();
    const typename ArrTraits::Circular_arc_point_2 & source = arc.source();
    const typename ArrTraits::Circular_arc_point_2 & target = arc.target();

    double asource = std::atan2( -to_double(source.y() - center.y()),
                 to_double(source.x() - center.x())); 
    double atarget = std::atan2( -to_double(target.y() - center.y()),
                 to_double(target.x() - center.x()));

    std::swap(asource, atarget);
    double aspan = atarget - asource;

    if(aspan < 0.)
      aspan += 2 * CGAL_PI;

    const double coeff = 180*16/CGAL_PI;
    qp->drawArc(convert(circ.bbox()), 
        (int)(asource * coeff), 
             (int)(aspan * coeff));
    return *this;
  }

  ArrPainterOstream& operator<<(const Line_arc_2<Kernel>& arc)
  {
    // std::clog << "Line_arc_2" << std::endl;        
 
    (*this) << Segment_2<Kernel>(Point_2<Kernel>(to_double(arc.source().x()), to_double(arc.source().y())),
                Point_2<Kernel>(to_double(arc.target().x()), to_double(arc.target().y())));
     return *this;
  }

/*
  void draw_parabola_segment(const  Point_2<Kernel>& center, const Line_2<Kernel>& line, 
                 const  Point_2<Kernel>& source, const Point_2<Kernel>& target)
  {
    const Point_2<Kernel> proj_source = line.projection(source);
    const Point_2<Kernel> proj_target = line.projection(target);
    const Point_2<Kernel> intersection = circumcenter(proj_source,
                         proj_target,
                         center);
    // Property: "intersection" is the intersection of the two tangent
    // lines in source and target.
    QPainterPath path;
    path.moveTo(convert(source));
    path.quadTo(convert(intersection), convert(target));
    qp->drawPath(path);
  }
  */ 
};

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ARR_PAINTER_OSTREAM_H
