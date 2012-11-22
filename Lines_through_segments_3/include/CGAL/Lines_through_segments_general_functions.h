// Copyright (c) 2010  Tel-Aviv University (Israel).
// All rights reserved.if 0
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
// $Id: $
// 
//
// Author(s)     : Asaf Porat          <asafpor1@post.tau.ac.il>

#ifndef LINES_THROUGH_SEGMENTS_GENERAL_FUNCTIONS_H
#define LINES_THROUGH_SEGMENTS_GENERAL_FUNCTIONS_H

#include <string>

#include <CGAL/basic.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Lines_through_segments_bounded_seg.h>
#include <CGAL/Lines_through_segments_bounded_segs_vector.h>
#include <CGAL/Lines_through_segments_traits_2_adapt.h>
#include <CGAL/Lines_through_segments_point_adapt.h>
#include <CGAL/Lines_through_segments_isolated_points.h>

namespace CGAL {


template <typename Lines_through_segments_traits_3>
class Lines_through_segments_general_functions
{
public:

   enum {
      CGAL_QUERY_SUCCEED = 0,
      CGAL_QUERY_FAILD = 1,
   };

  typedef Lines_through_segments_traits_3               Traits_3;
private:
  typedef typename Traits_3::Traits_arr_on_plane_2      Traits_arr_on_plane_2;
  typedef typename Traits_3::Traits_arr_on_sphere_2     Traits_arr_on_sphere_2;
  typedef typename Traits_3::Alg_kernel                 Alg_kernel;
  typedef typename Traits_3::Rational_kernel            Rational_kernel;

  typedef typename Alg_kernel::FT                       Algebraic;
  typedef typename Rational_kernel::FT                  Rational;

  typedef Polyhedron_3<Rational_kernel>                 Rational_polyhedron_3;

  typedef typename Alg_kernel::Point_3                  Alg_point_3;
  typedef typename Alg_kernel::Line_3                   Alg_line_3;
  typedef typename Alg_kernel::Segment_3                Alg_segment_3;
  typedef typename Alg_kernel::Plane_3                  Alg_plane_3;
    
  typedef typename Rational_kernel::Point_3             Rational_point_3;
  typedef typename Rational_kernel::Line_3              Rational_line_3;
  typedef typename Rational_kernel::Segment_3           Rational_segment_3;
  typedef typename Rational_kernel::Plane_3             Rational_plane_3;
  typedef typename Rational_kernel::Point_2             Rational_point_2;
      
  typedef typename Traits_arr_on_plane_2::Curve_2       Rational_arc_2;
      
  typedef Lines_through_segments_point_adapt_2<Traits_3, typename Traits_arr_on_plane_2::Point_2,Algebraic> Point_2;

public:
  bool has_the_same_supporting_but_not_intersect(const Rational_segment_3& s1,
                                                 const Rational_segment_3& s2)
  {
    return
      (!s1.has_on(s2.source()) && !s1.has_on(s2.target()) &&
       s1.supporting_line().has_on(s2.source()) && 
       s1.supporting_line().has_on(s2.target()));
  }
      
  /*************************************************************
   * The following function gets a intersection point on the arrangement, 
   * and returns the line that crosses S1 at S1_t and S2 and S2_t.
   **************************************************************/
  template <typename Point_on_arr_2>
  bool add_line_to_output(const Rational_segment_3& s1,
                          const Rational_segment_3& s2,
                          bool s1_s2_intersect,
                          const Point_on_arr_2& intersection_point,
                          Alg_line_3& common_line)
  {

    Algebraic S1_t;
    Algebraic S2_t;
    Lines_through_segments_get_algebraic_number_adapt<
    Traits_3>
      get_algebraic_number_adapt;
        
#if ARR_ON_SUR_DEBUG
    std::cout << "\n\n"
              << "*************************************NEW LINE PLANE"
              << "*****************\n\n"
              << std::endl;
#endif        
    Point_2 int_adapt(intersection_point);
    S1_t = get_algebraic_number_adapt(int_adapt.x());
    S2_t = get_algebraic_number_adapt(int_adapt.y());
         
#if ARR_ON_SUR_DEBUG
    std::cout << "S1_t = " << S1_t << std::endl;
    std::cout << "S2_t = " << S2_t << std::endl;
#endif
    int status = 
      get_line_from_intersection_point(S1_t, S2_t, s1, s2, common_line);
         
    if (status != CGAL_QUERY_SUCCEED)
    {
      /* In case S1 and S2 intersect, the lines that passes through the 
       * intersection point are handled at the arr on sphere case.
       */
      if (!s1_s2_intersect)
      {
        CGAL_error_msg("Unexpected error");
      }
      return false;
    }
    return true;
  }
  /*************************************************************
   * The following function gets a line and a segment at 3 space,
   * and return true if the line intersect with the segment.
   *
   * The segment is represented by a line where the segement is from the first 
   * point on the line and the second point on the line
   *************************************************************/
  template <typename Kernel, typename Line_3>
  bool do_intersect_line_segment(const Line_3& L_1,
                                 const Rational_segment_3& S_2,
                                 const Kernel* kernel)
  {
   
    typedef typename Kernel::Point_3 Point_3;
    typedef typename Kernel::Segment_3 Segment_3;
    typedef typename Kernel::FT FT;
   
    Point_3 p1 = Point_3(S_2.source().x(),
                         S_2.source().y(),
                         S_2.source().z());

    if (S_2.source() == S_2.target())
    {
      if (!L_1.has_on(p1))
      {
        return false;
      }
    }
    else
    {
      Point_3 p2 = Point_3(S_2.target().x(), S_2.target().y(), S_2.target().z());
      Line_3 L_2_alg = Line_3(p1,p2);
      
      CGAL::Object result;
      Point_3 ipoint;
      
      result = kernel->intersect_3_object()(L_1,L_2_alg);
        
      
      if (CGAL::assign(ipoint, result))
      {
        FT x1 = S_2.source().x();
        FT y1 = S_2.source().y();
        FT z1 = S_2.source().z();
      
        FT x2 = S_2.target().x();
        FT y2 = S_2.target().y();
        FT z2 = S_2.target().z();
      
        Point_3 alg_point0(x1,y1,z1);
        Point_3 alg_point1(x2,y2,z2);
        Segment_3 S_2(alg_point0,alg_point1);


#if LINES_DEBUG
        std::cout << "ipoint = " << ipoint << std::endl;
        if ((x2-x1) != 0)
          std::cout<<"L3.t = "<<((ipoint.x()-x1)/(x2-x1))<<std::endl;
        if ((y2-y1) != 0)
          std::cout<<"L3.t = "<<((ipoint.y()-y1)/(y2-y1))<<std::endl;
        if ((z2-z1) != 0)
          std::cout<<"L3.t = "<<((ipoint.z()-z1)/(z2-z1))<<std::endl;
#endif

        if (!S_2.has_on(ipoint))
        {
#if LINES_DEBUG
          std::cout <<"Point on line but not on segment"<<std::endl;
#endif         
          return false;
        }
      }
      else
      {
        Line_3 iline;
        if (CGAL::assign(iline, result))
        {
          return true;
        }
        else
        {
          return false;
        }
      }
    }
   
    return true;
  }

  template <typename Kernel, typename Line_3>
  bool do_intersect_line_segment(const Line_3& L_1,
                                 const Rational_polyhedron_3& p,
                                 const Kernel* kernel)
  {
    return true;
  }
      
  bool do_intersect_segments(const Rational_segment_3& S1,
                             const Rational_segment_3& S2,
                             Rational_point_3& intersection_point_S1S2,
                             const Rational_kernel* rat_kernel)
  {
    CGAL::Object result = 
      rat_kernel->intersect_3_object()(S1.supporting_line(),
                                       S2.supporting_line());
    if (CGAL::assign(intersection_point_S1S2, result))
    {
      if (S1.has_on(intersection_point_S1S2) && 
          S2.has_on(intersection_point_S1S2))
      {
        return true;
      }
    }
    return false;
  }

  /*************************************************************
   * The following function gets the scalar to the line vector from
   * a point on the line.
   **************************************************************/
  template <typename Line_3>
  void get_scalar_from_point_on_line(const Rational_point_3& ipoint,
                                     const Line_3& line,
                                     Rational& result)
  {
    if (line.point(1).x() - line.point(0).x() != 0)
    {
      result = ((ipoint.x()-line.point(0).x())/
                (line.point(1).x()-line.point(0).x()));
    }
    else if (line.point(1).y() - line.point(0).y() != 0)
    {
      result = ((ipoint.y()-line.point(0).y())/
                (line.point(1).y()-line.point(0).y()));
    }
    else if (line.point(1).z() - line.point(0).z() != 0)
    {
      result = ((ipoint.z()-line.point(0).z())/
                (line.point(1).z()-line.point(0).z()));
    }
  }

  /*************************************************************
   * The following function gets the values of two points one on S1 and the
   * other on S2 (S1_t and S2_t) and return the line that passes through these 
   * 2 points.
   * If the points are the same point the  function return CGAL_QUERY_FAILD
   *  otherwise the function return CGAL_QUERY_SUCCEED and the line that passes
   * through the 2 points.
   * 
   * Input:
   *         S1_t      - This values determines a single point on S1.
   *         S2_t      - This values determines a single point on S1.
   *         S1        - Line in 3 space.
   *         S2        - Line in 3 space.
   *
   * Output:
   *         common_line - A line which passes through the 2 lines.
   *
   *************************************************************/
   
  int get_line_from_intersection_point(const Algebraic& S1_t,
                                       const Algebraic& S2_t,
                                       const Rational_segment_3& S1,
                                       const Rational_segment_3& S2,
                                       Alg_line_3& common_line)
  {
    Rational_point_3 S1_P1 = S1.source();
    Rational_point_3 S1_P2 = S1.target();
   
    Rational_point_3 S2_P1 = S2.source();
    Rational_point_3 S2_P2 = S2.target();

    Algebraic x1 = S1_P1.x() + S1_t * (S1_P2.x() - S1_P1.x());
    Algebraic y1 = S1_P1.y() + S1_t * (S1_P2.y() - S1_P1.y());
    Algebraic z1 = S1_P1.z() + S1_t * (S1_P2.z() - S1_P1.z());
   
    Algebraic x2 = S2_P1.x() + S2_t * (S2_P2.x() - S2_P1.x());
    Algebraic y2 = S2_P1.y() + S2_t * (S2_P2.y() - S2_P1.y());
    Algebraic z2 = S2_P1.z() + S2_t * (S2_P2.z() - S2_P1.z());
   
    Alg_point_3 P1(x1,y1,z1);
    Alg_point_3 P2(x2,y2,z2);

    /* The point is the intersection point of S1 and S2. 
       Hence the common line will have to be calculated using S3.*/
    if (P1.x() == P2.x() && P1.y() == P2.y() && P1.z() == P2.z())
    {
#if LINES_DEBUG
      std::cout << "POINT_ON_THE_INTERSCTION" << std::endl;
      std::cout << "P1 = (" << P1.x() << "," << P1.y() << "," << P1.z()
                << std::endl;
      std::cout << "P2 - (" << P2.x() << "," << P2.y() << "," << P2.z()
                << std::endl;
#endif
      return CGAL_QUERY_FAILD;
    }
   
    common_line = Alg_line_3(P1,P2);
  
    return CGAL_QUERY_SUCCEED;
  }

   int get_line_from_intersection_point(const Rational& S1_t,
                                        const Rational& S2_t,
                                        const Rational_segment_3& S1,
                                        const Rational_segment_3& S2,
                                        Rational_line_3& common_line)
   {
      Rational_point_3 S1_P1 = S1.source();
      Rational_point_3 S1_P2 = S1.target();
      
      Rational_point_3 S2_P1 = S2.source();
      Rational_point_3 S2_P2 = S2.target();
      
    Rational x1 = S1_P1.x() + S1_t * (S1_P2.x() - S1_P1.x());
    Rational y1 = S1_P1.y() + S1_t * (S1_P2.y() - S1_P1.y());
    Rational z1 = S1_P1.z() + S1_t * (S1_P2.z() - S1_P1.z());
   
    Rational x2 = S2_P1.x() + S2_t * (S2_P2.x() - S2_P1.x());
    Rational y2 = S2_P1.y() + S2_t * (S2_P2.y() - S2_P1.y());
    Rational z2 = S2_P1.z() + S2_t * (S2_P2.z() - S2_P1.z());
   
    Rational_point_3 P1(x1,y1,z1);
    Rational_point_3 P2(x2,y2,z2);

    /* The point is the intersection point of S1 and S2. 
       Hence the common line will have to be calculated using S3.*/
    if (P1.x() == P2.x() && P1.y() == P2.y() && P1.z() == P2.z())
    {
#if LINES_DEBUG
      std::cout << "POINT_ON_THE_INTERSCTION" << std::endl;
      std::cout << "P1 = (" << P1.x() << "," << P1.y() << "," << P1.z()
                << std::endl;
      std::cout << "P2 - (" << P2.x() << "," << P2.y() << "," << P2.z()
                << std::endl;
#endif
      return CGAL_QUERY_FAILD;
    }
   
    common_line = Rational_line_3(P1,P2);
  
    return CGAL_QUERY_SUCCEED;
  }

  /*************************************************************
   * The following function gets a point on the sphere
   * and return the line that passes through the 4 lines.
   * 
   * Input:
   *         point     - Point on the sphere.
   *         cpoint    - center of the sphere.
   *
   * Output:
   *         common_line - A line which passes through the 4 lines.
   *
   *************************************************************/

  void get_line_from_intersection_point(const Rational_point_3& cpoint,
                                        const Rational& x,
                                        const Rational& y,
                                        const Rational& z,
                                        Rational_line_3& common_line)
  {
    common_line = Rational_line_3(Rational_point_3(cpoint.x(),cpoint.y(),cpoint.z()),
                                  Rational_point_3(cpoint.x()+x,cpoint.y()+y,cpoint.z()+z));
  }
};

} //namespace CGAL

#endif //LINES_THROUGH_SEGMENTS_GENERAL_FUNCTIONS_H
