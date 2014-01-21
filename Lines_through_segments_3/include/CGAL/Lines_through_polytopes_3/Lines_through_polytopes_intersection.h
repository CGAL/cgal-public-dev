// Copyright (c) 2010  Tel-Aviv University (Israel).
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
// $Id: $
// 
//
// Author(s)     : Asaf Porat          <asafpor1@post.tau.ac.il>

#ifndef LINE_THROUGH_POLYTOPES_INTERSECTION_H
#define LINE_THROUGH_POLYTOPES_INTERSECTION_H

#include <CGAL/Polygon_2.h>
#include <CGAL/Lines_through_segments_general_functions.h>
#include <CGAL/Lines_through_segments_arr_gen_func.h>

namespace CGAL {
   
template <typename Lines_through_segments_traits_3>
class Lines_through_polytopes_intersection
{
  typedef typename Lines_through_segments_traits_3::Algebraic_NT Algebraic;
  typedef typename Lines_through_segments_traits_3::Rational_NT  Rational;

  typedef typename Lines_through_segments_traits_3::Traits_arr_on_plane_2
  Traits_arr_on_plane_2;
  typedef typename Lines_through_segments_traits_3::Traits_arr_on_sphere_2
  Traits_arr_on_sphere_2;
  typedef typename Lines_through_segments_traits_3::Alg_kernel Alg_kernel;
  typedef typename Lines_through_segments_traits_3::Rational_kernel Rational_kernel;
  typedef typename Alg_kernel::FT Algebraic;
  typedef typename Rational_kernel::FT  Rational;

         
  typedef typename Alg_kernel::Point_3                    Alg_point_3;
  typedef typename Alg_kernel::Line_3                     Alg_line_3;
  typedef typename Alg_kernel::Segment_3                  Alg_segment_3;
  typedef typename Alg_kernel::Plane_3                    Alg_plane_3;
  typedef typename Alg_kernel::Point_2                    Alg_point_2;
  typedef typename Alg_kernel::Line_2                     Alg_line_2;
         
  typedef typename Rational_kernel::Point_3               Rational_point_3;
  typedef typename Rational_kernel::Line_3                Rational_line_3;
  typedef typename Rational_kernel::Segment_3             Rational_segment_3;
  typedef typename Rational_kernel::Plane_3               Rational_plane_3;
  typedef typename Rational_kernel::Point_2               Rational_point_2;

  typedef CGAL::Polyhedron_3<Rational_kernel>             Rational_polyhedron_3;
  typedef CGAL::Polygon_2<Rational_kernel>                Rational_polygon_2;
  typedef CGAL::Polygon_2<Alg_kernel>                     Alg_polygon_2;
  typedef typename Rational_polyhedron_3::Halfedge_handle Halfedge_handle;
  typedef typename Rational_polyhedron_3::Edge_const_iterator   Edge_iterator;
  typedef typename Rational_polyhedron_3::Vertex_handle   Vertex_handle;
  typedef typename Rational_polyhedron_3::Facet_iterator  Facet_iterator;
         
  typedef typename Rational_polyhedron_3::Halfedge_around_vertex_const_circulator
  Halfedge_around_vertex_const_circulator;
         
  typedef typename Rational_polyhedron_3::Halfedge_around_facet_const_circulator
  Halfedge_around_facet_const_circulator;
         
private:
  const Alg_kernel* m_alg_kernel;
         

public:
  Lines_through_polytopes_intersection(const Alg_kernel* alg_kernel)
  {
    m_alg_kernel = alg_kernel;
  }
         
  CGAL::Bounded_side operator()(Rational_polyhedron_3& poly,const Alg_line_3& line)
  {
    Alg_point_3 point_on_boundary;
    bool point_on_boundary_valid = false;
            
    CGAL::Bounded_side result = CGAL::ON_UNBOUNDED_SIDE;
    for (Facet_iterator facet_it = poly.facets_begin();
         facet_it != poly.facets_end();
         facet_it++)
    {
      Halfedge_around_facet_const_circulator hafc = facet_it->facet_begin();
      Rational_point_3 p1(hafc->vertex()->point());
      hafc++;
      Rational_point_3 p2(hafc->vertex()->point());
      hafc++;
      Rational_point_3 p3(hafc->vertex()->point());
               
      Alg_point_3 p1_alg(p1.x(),p1.y(),p1.z());
      Alg_point_3 p2_alg(p2.x(),p2.y(),p2.z());
      Alg_point_3 p3_alg(p3.x(),p3.y(),p3.z());
               

      std::cout << p1_alg << std::endl;
      std::cout << p2_alg << std::endl;
      std::cout << p3_alg << std::endl;
               
      Alg_plane_3 plane(p1_alg, p2_alg, p3_alg);
               
      /*************DEBUG*******************/
      hafc = facet_it->facet_begin();
      do
      { 
        Alg_point_3 temp_p3(hafc->vertex()->point().x(),
                            hafc->vertex()->point().y(),
                            hafc->vertex()->point().z());
                  
        hafc++;
                     
        Alg_line_3 line_of_edge(temp_p3,
                                Alg_point_3(hafc->vertex()->point().x(),
                                            hafc->vertex()->point().y(),
                                            hafc->vertex()->point().z()));
                     
        Alg_line_3 iline;
        CGAL::Object res_line = m_alg_kernel->intersect_3_object()(plane,line_of_edge);
        if (!CGAL::assign(iline,res_line))
        {
          std::cout << "ERRROR line and Plane do not intersect = " << line_of_edge << std::endl;
          CGAL_error_msg("ERROR");
        }
        else
          std::cout << "GOOD" << std::endl;
      }
      while(hafc != facet_it->facet_begin());
               
               
      /*************END DEBUG***************/
      /* There is a bug at the polygon_2 class.
         If ipoint is on one of the lines of the polygon edges 
         but no the edge,
         the query bounded_side returns on bounded side, while
         the result should be on unbounded side.
         In this case we skip the facet.
      */
                  
      bool skip_facet = false;
               
               
      Alg_point_3 ipoint;
      Alg_line_3 iline;
               
      CGAL::Object int_result = 
        m_alg_kernel->intersect_3_object()(plane, line);
               
      if (CGAL::assign(ipoint, int_result))
      {
        /* Draw a line from ipoint to some point at polus infinity on
           the plane. 
           If the line intersect one line than the point is inside
           the facet, else if the point is on one of the lines than
           the point is on the boundary of the facet.
           else the point is outside the facet.
        */
        std::cout << "-----------------------------------------------------\n"
                  << "NEW FACET"
                  << std::endl;
                  
        std::list<Alg_point_2> v_points;
        hafc = facet_it->facet_begin();
        do
        { 
          Alg_point_3 temp_p3(hafc->vertex()->point().x(),
                              hafc->vertex()->point().y(),
                              hafc->vertex()->point().z());
                     
          Alg_point_2 temp_p(plane.to_2d(temp_p3));
          v_points.push_back(temp_p);
          hafc++;
                     
          Alg_segment_3 temp_seg(temp_p3,
                                 Alg_point_3(hafc->vertex()->point().x(),
                                             hafc->vertex()->point().y(),
                                             hafc->vertex()->point().z()));
                     
          if (m_alg_kernel->do_intersect_3_object()(Alg_line_3(temp_seg.source(),temp_seg.target()), line) &&
              !temp_seg.has_on(ipoint))
          {
            Alg_line_3 line_of_edge(temp_seg.source(),temp_seg.target());
            skip_facet = true;
            std::cout << "SKIP FACET" << std::endl;
            std::cout << "ipoint = " << ipoint << std::endl;
            std::cout << "temp_seg  = " << temp_seg << std::endl;
            std::cout << "line.has_on(ipoint) = " << line.has_on(ipoint) << std::endl;
            std::cout << "line_of_edge.has_on(ipoint) = " << line_of_edge.has_on(ipoint) << std::endl;
            std::cout << "Plane = " << plane << std::endl;
                        
            CGAL::Object res_line = m_alg_kernel->intersect_3_object()(plane,line_of_edge);
            if (!CGAL::assign(iline,res_line))
              std::cout << "ERRROR line and Plane do not intersect = " << line_of_edge << std::endl;
            break;
          }

          if (temp_seg.has_on(ipoint))
          {
            Alg_point_2 temp_p(plane.to_2d(ipoint));
            Alg_line_2 temp_line111(plane.to_2d(temp_seg.source()),
                                    plane.to_2d(temp_seg.target()));

                        
            if (!temp_seg.has_on(ipoint))
              std::cout << "WWWWWWWWEEEEEEEEEEEEEEEEEEERRRRRRRRRRRRRRRRRRRRROOOOOOOOOOOOOORRRRRRRRRR" << std::endl;

            if (!line.has_on(ipoint))
              std::cout << "MMMMMMEEEEEEEEEEEEEEEEEEERRRRRRRRRRRRRRRRRRRRROOOOOOOOOOOOOORRRRRRRRRR" << std::endl;

            if (!temp_line111.has_on(temp_p))
              std::cout << "EEEEEEEEEEEEEEEEEEERRRRRRRRRRRRRRRRRRRRROOOOOOOOOOOOOORRRRRRRRRR" << std::endl;

          }
                     
          std::cout << "temp_seg.has_on(ipoint) = "
                    << temp_seg.has_on(ipoint)
                    <<std::endl;

          std::cout << "Do intersect edge line = "
                    << m_alg_kernel->do_intersect_3_object()(Alg_line_3(temp_seg.source(),temp_seg.target()), line)
                    <<std::endl;
        }
        while(hafc != facet_it->facet_begin());
                  
        if (!skip_facet)
        {
          Alg_polygon_2 pgn(v_points.begin(),v_points.end());
          std::cout << "polygon_2 = " << pgn << std::endl;
                     
          Alg_point_2 temp_p(plane.to_2d(ipoint));
                     
          if (!pgn.is_convex())
            std::cout << change_color(CGAL_BLUE,"ERROR NOT CONVEX") << std::endl;
                     
          if (!pgn.is_simple())
            std::cout << change_color(CGAL_BLUE,"ERROR NOT CONVEX") << std::endl;
          switch (pgn.bounded_side(temp_p))
          {
           case CGAL::ON_BOUNDED_SIDE:
            return CGAL::ON_BOUNDED_SIDE;

           case CGAL::ON_BOUNDARY:
            {
              if (point_on_boundary_valid && ipoint != point_on_boundary)
                return CGAL::ON_BOUNDED_SIDE;
              else
                point_on_boundary_valid = true;
              point_on_boundary = ipoint;
              result = CGAL::ON_BOUNDARY;
            }
            break;
                        
           case  CGAL::ON_UNBOUNDED_SIDE:
           default:
            break;
          }
        }
      }
      else if (CGAL::assign(iline, int_result))
      {
        std::cout << change_color(CGAL_RED,"TODO") << std::endl;
      }
    }
    return result; 
  }
};
   
} //namespace CGAL

#endif /*LINE_THROUGH_POLYTOPES_INTERSECTION_H*/
