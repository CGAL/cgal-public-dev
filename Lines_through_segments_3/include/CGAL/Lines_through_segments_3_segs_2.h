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

#ifndef LINE_THROUGH_SEGMENTS_3_SEGS_2_H
#define LINE_THROUGH_SEGMENTS_3_SEGS_2_H

#include <CGAL/Lines_through_segments_general_functions.h>
#include <CGAL/Lines_through_segments_arr_gen_func.h>
#include <CGAL/Lines_through_segments_3/observer.h>
#include <CGAL/Lines_through_segments_arr_plane_faces.h>
#include <CGAL/Lines_through_segments_traits_2_adapt.h>
#include <CGAL/Lines_through_segments_point_adapt.h>

/*************************************************************
 * This class purpose is to find all lines that passes through 3 segments in 
 * the plane.
 * The class returns an arrangement that represents the all the lines that 
 * passes through the 3 segments parametrized by the first to segments.
 *
 *************************************************************/
namespace CGAL {

template <typename Traits_3_,
          typename Lines_through_segments_arr_plane_faces,
          typename Isolated_points_on_plane,
          typename Point_on_plane_and_line_pair,
          typename Arc_end_points,
          typename With_segments>
class Lines_through_segments_3_segs_2 {
public:
  typedef Traits_3_                                      Traits_3;
private:  
  typedef typename Traits_3::Alg_kernel                  Alg_kernel;
  typedef typename Traits_3::Rational_kernel             Rational_kernel;
  typedef typename Traits_3::Traits_arr_on_plane_2       Traits_arr_on_plane_2;
      
  typedef typename Alg_kernel::FT                        Algebraic;
  typedef typename Rational_kernel::FT                   Rational;

  typedef typename Alg_kernel::Point_3                   Alg_point_3;
  typedef typename Alg_kernel::Line_3                    Alg_line_3;
  typedef typename Alg_kernel::Segment_3                 Alg_segment_3;
  typedef typename Alg_kernel::Plane_3                   Alg_plane_3;
      
  typedef typename Rational_kernel::Point_3              Rational_point_3;
  typedef typename Rational_kernel::Line_3               Rational_line_3;
  typedef typename Rational_kernel::Segment_3            Rational_segment_3;
  typedef typename Rational_kernel::Plane_3              Rational_plane_3;
  typedef typename Rational_kernel::Point_2              Rational_point_2;
      
  /***************************************************/
  /*    Arrangement on plane typedefs.               */
  /***************************************************/
      
  /* Extended each edge its creator line, Extend each vertex with color.*/
  typedef Lines_through_segments_arr_ext_dcel<Traits_arr_on_plane_2, Rational_segment_3> 
    Dcel_on_plane;
  typedef CGAL::Arrangement_2<Traits_arr_on_plane_2, Dcel_on_plane>
    Arrangement_on_plane_2;

  typedef Lines_through_segments_point_adapt_2<Traits_3,
    typename Traits_arr_on_plane_2::Point_2, Algebraic>  Point_2;

  typedef typename Traits_arr_on_plane_2::Curve_2        Rational_arc_2;
  typedef Lines_through_segments_arr_observer<Arrangement_on_plane_2>
    Lines_through_segments_arr_observer_on_plane;

  typedef typename Arrangement_on_plane_2::Halfedge_handle 
                                                         Halfedge_handle;    
  typedef typename Arrangement_on_plane_2::Face_iterator Face_iterator;
  typedef typename Arrangement_on_plane_2::Face_handle   Face_handle;
  typedef Lines_through_segments_traits_on_plane_adapt<Traits_3>
    Traits_2_adapt;
         
private:
  const Rational_segment_3* m_S1;
  const Rational_segment_3* m_S2;
  const Rational_kernel *m_rat_kernel;
  const Alg_kernel* m_alg_kernel;
  Lines_through_segments_traits_on_plane_adapt<Traits_3> m_traits_2_adapt;
  typedef Lines_through_segments_general_functions<Traits_3> LTS_g_func;
   
  LTS_g_func m_g_func;
  Lines_through_segments_arr_gen_func<Traits_3, With_segments> m_arr_g_func;
  bool m_bound_s1_s2;
    
public:
  Lines_through_segments_3_segs_2(){}

  /*************************************************************
   * The following constructor handles the case where the 3 segments are on
   * the same plane.
   * It gets a segment and adds a plane which represents all the lines that
   * passes through S3 and S1 and S2.
   * 
   * The function finds the lines of the end points of S1 with the end points
   * of S2, for each line it finds the intersection with S3.
   * with the end points of S3 we will get <= 4 points/Curves 
   * (we will get points only for 2 end points (S1 and S2)
   * that are at the bottom of S1 and S2, or both on the upper part of S1 
   * and S2.
   *
   * Next, the function finds the convex hull of these 6 objects, 
   * only two of them will create
   * the convex hull since the intermidiate points creates curves inside the
   * outer points. 
   * 
   * The function finds the two most outer points on S3 and insert the face 
   * between the two curves/points,
   * that were generated by these points. 
   *
   * Input:
   *      S3 - segment.
   *
   * Output:
   *
   *************************************************************/
  void operator() (const Rational_segment_3* S1,
                   const Rational_segment_3* S2,
                   const Rational_segment_3& S3,
                   bool bound_s1_s2,
                   const Rational_kernel* rat_kernel,
                   const Alg_kernel* alg_kernel,
                   const Rational_point_3& intersection_point_S1S2,
                   bool S1_S2_intersect,
                   Isolated_points_on_plane& isolated_points_on_plane,
                   Lines_through_segments_arr_plane_faces& arr_plane_faces)
  {
    m_bound_s1_s2 = bound_s1_s2;
         
    m_arr_g_func = 
      Lines_through_segments_arr_gen_func<Traits_3, With_segments>(alg_kernel);
    m_S1 = S1;
    m_S2 = S2;
    m_rat_kernel = rat_kernel;
    m_alg_kernel = alg_kernel;
        
    /* An isolated point is a point on one of the square corners. 
       It represents an arc that passes through the corner, and bounds a 
       face with counter > 0.
    */
    Isolated_points_on_plane local_isolated_points;
    Isolated_points_on_plane local_isolated_points_f; 
    Isolated_points_on_plane local_isolated_points_s;
    typename std::list<Rational_arc_2>  arcs[2];
    std::list<Arc_end_points > end_points[2];
    Rational min_point = 1;
    Rational max_point = 0;
    CGAL::Object result;
        
    /* Get the arcs of the face that represents all the lines that passes 
       through S1, S2 and S3.
       Only the arcs at the arrangement are found at this stage. Later the 
       bounded linear segments on
       the the square [0,0] ,[1,1] are added in order to close the faces.
            
       Note:
       In degenerate cases the arcs can be "shrank" to points on the sides 
       of the square [0,0], [1,1].
    */

    get_arr_arcs(S3,
                 arcs,
                 end_points,
                 intersection_point_S1S2,
                 S1_S2_intersect,
                 min_point,
                 max_point,
                 local_isolated_points_f,
                 local_isolated_points_s);
         
    typename std::list<Rational_arc_2> arcs_to_insert;
    get_arr_faces_2_intersections(S3, arcs, arr_plane_faces);
    return;
  }

private:
  /*************************************************************
   * The function finds the combination of lines at the plane that passes 
   * throuh the end points of S1 with the end points of S2,
   * for each line, it finds the intersection with S3.
   * with the end points of S3 we will get <= 4 points/Curves 
   * (we will get points only for 2 end points (S1 and S2)
   * that are at the bottom of S1 and S2, or both on the upper part of S1 
   * and S2.
   *************************************************************/

  void get_arr_arcs(const Rational_segment_3& S3,
                    typename std::list<Rational_arc_2>  arcs[2],
                    std::list<Arc_end_points > end_points[2],
                    const Rational_point_3& intersection_point_S1S2,
                    bool S1_S2_intersect,
                    Rational& min_point,
                    Rational& max_point,
                    Isolated_points_on_plane& local_isolated_points_f,
                    Isolated_points_on_plane& local_isolated_points_s)
  {
    Rational_line_3 S1_S2[4];
    Rational_point_3 points[2];
    Rational_point_3 ipoint;
    CGAL::Object result;

    S1_S2[0] = Rational_line_3(m_S1->source(),m_S2->source());
    S1_S2[1] = Rational_line_3(m_S1->source(),m_S2->target());
    S1_S2[2] = Rational_line_3(m_S1->target(),m_S2->source());
    S1_S2[3] = Rational_line_3(m_S1->target(),m_S2->target());         

#if ARR_ON_SUR_DEBUG
    std::cout<<"3 lines on the same plane case" << std::endl;
#endif

    /* Find all the lines that passes through the start point of S3 and S1
     * and S2.
     */
    m_arr_g_func.get_all_lines_through_point_and_2_lines
      (*m_S1,*m_S2,S3,
       m_bound_s1_s2, S3.source(),
       *m_rat_kernel,intersection_point_S1S2,
       local_isolated_points_f, std::back_inserter(arcs[0]),
       true,&end_points[0],
       S1_S2_intersect);
         
    if (arcs[0].size() > 0)
    {
      points[0] = S3.source();
      min_point = 0;
    }

    /* Find all the lines that passes through the end point of S3 and S1 
       and S2. */
    m_arr_g_func.get_all_lines_through_point_and_2_lines
      (*m_S1,*m_S2,S3,
       m_bound_s1_s2, S3.target(),
       *m_rat_kernel,intersection_point_S1S2,
       local_isolated_points_s,
       std::back_inserter(arcs[1]), true, &end_points[1],
       S1_S2_intersect);
         
    if (arcs[1].size() > 0)
    {
      points[1] = S3.target();
      max_point = 1;
    }
         
    /* If all the lines that passes through the 3 lines doesn't pass 
       through the source or target of S3. */
    if (min_point != 0 || max_point != 1)
    {
      Rational scalar;
      for (int ii = 0; ii < 4; ++ii)
      {
        result = m_rat_kernel->intersect_3_object()(S1_S2[ii],
                                                    S3.supporting_line());
        if (CGAL::assign(ipoint, result) && S3.has_on(ipoint))
        {
          m_g_func.get_scalar_from_point_on_line(ipoint,
                                                 S3.supporting_line(),scalar);
          if (scalar < min_point)
          {
            min_point = scalar;
            points[0] = ipoint;
          }
             
          if (scalar > max_point)
          {
            max_point = scalar;
            points[1] = ipoint;
          }
        }
      }
       
      if (min_point != 0 && min_point < 1)
      {
        m_arr_g_func.get_all_lines_through_point_and_2_lines
          (*m_S1,*m_S2,S3, 
           m_bound_s1_s2, points[0],
           *m_rat_kernel,
           intersection_point_S1S2,
           local_isolated_points_f,
           std::back_inserter(arcs[0]), true, &end_points[0],
           S1_S2_intersect);
      }
            
      CGAL_assertion(arcs[0].size() <= 2);
      if (max_point != 1 && max_point > 0)
      {
        m_arr_g_func.get_all_lines_through_point_and_2_lines
          (*m_S1,*m_S2,S3,          
           m_bound_s1_s2, points[1],
           *m_rat_kernel,
           intersection_point_S1S2,
           local_isolated_points_s,
           std::back_inserter(arcs[1]), true, &end_points[1],
           S1_S2_intersect);
      }
      CGAL_assertion(arcs[1].size() <= 2);
    }
  }

private:
      
  /*************************************************************
   * Creates an arragenment from the arcs at arcs_to_insert list.
   *************************************************************/
    
private:
  /*************************************************************
   * The following function handle the case were S1 and S2 intersect, and
   * S3 intersect one or both of S1 and S2.
   **************************************************************/
  void get_arr_faces_2_intersections(const Rational_segment_3& S3,
                                     typename std::list<Rational_arc_2> arcs[2],
                                     Lines_through_segments_arr_plane_faces&
                                     arr_plane_faces)
  {
    typename std::list<Rational_arc_2> arcs_to_insert;
    Rational_arc_2 arc;
         
    /* Push the 4 lines that bound the square [0,0], [1,1]. */
         
    m_traits_2_adapt.create_segment_on_plane_arr(arc,
                                                 Rational_point_2(Rational(0),
                                                                  Rational(0)),
                                                 Rational_point_2(Rational(0),
                                                                  Rational(1)),
                                                 &S3);
    arcs_to_insert.push_back(arc);
         

    m_traits_2_adapt.create_segment_on_plane_arr(arc,
                                                 Rational_point_2(Rational(0),
                                                                  Rational(1)),
                                                 Rational_point_2(Rational(1),
                                                                  Rational(1)),
                                                 &S3);
    arcs_to_insert.push_back(arc);
        
        
    m_traits_2_adapt.create_segment_on_plane_arr(arc,
                                                 Rational_point_2(Rational(1),
                                                                  Rational(1)),
                                                 Rational_point_2(Rational(1),
                                                                  Rational(0)),
                                                 &S3);
    arcs_to_insert.push_back(arc);
        
        
    m_traits_2_adapt.create_segment_on_plane_arr(arc,
                                                 Rational_point_2(Rational(1),
                                                                  Rational(0)),
                                                 Rational_point_2(Rational(0),
                                                                  Rational(0)),
                                                 &S3);
    arcs_to_insert.push_back(arc);

    /* Push the arcs - for each verify that it is not one of the square 
       [0,0],[1,1] sides,
       These arcs will be inserted later.*/
    typename std::list<Rational_arc_2>::iterator arcs_it;
         
    for (arcs_it = arcs[0].begin();arcs_it != arcs[0].end();++arcs_it)
    {
      if (!is_segment_a_side(*arcs_it))
        arcs_to_insert.push_back(*arcs_it);
    }

    for (arcs_it = arcs[1].begin();arcs_it != arcs[1].end();++arcs_it)
    {
      if (!is_segment_a_side(*arcs_it))
      {
        arcs_to_insert.push_back(*arcs_it);
      }
            
    }  
    remove_faces(S3,arr_plane_faces,arcs_to_insert);
  }

  /* Get over all of the faces and for each face determine weather its a 
   *  part of the arrangement.
   */

  void remove_faces(const Rational_segment_3& S3,
                    Lines_through_segments_arr_plane_faces& arr_plane_faces,
                    typename std::list<Rational_arc_2> arcs_to_insert)
  {
    std::list<Halfedge_handle> edges_to_erase;
        
    Arrangement_on_plane_2* temp_arr_on_plane = 
      new Arrangement_on_plane_2();
    Lines_through_segments_arr_observer_on_plane* temp_obs_on_plane =
      new Lines_through_segments_arr_observer_on_plane(*temp_arr_on_plane);

    temp_obs_on_plane->set_is_plane(true);
    temp_obs_on_plane->set_last_inserted_segment(&S3);

    insert (*temp_arr_on_plane, arcs_to_insert.begin(), 
            arcs_to_insert.end());
         
    typename Arrangement_on_plane_2::Edge_iterator   eit;
    for (eit = temp_arr_on_plane->edges_begin();
         eit != temp_arr_on_plane->edges_end();
         ++eit)
    {
      Point_2 source(eit->source()->point());
      Point_2 target(eit->target()->point());
            
      if (is_on_square_boudns(eit) &&
          !is_line_at_arr_valid(source,
                                target,S3))
      {
        edges_to_erase.push_back(eit);
      }
    }
         
    typename std::list<Halfedge_handle>::iterator it;
    for(it = edges_to_erase.begin();
        it != edges_to_erase.end();
        ++it)
    {                
      temp_arr_on_plane->remove_edge(*it);
    }
#if ARR_ON_SUR_DEBUG        
    for (eit = temp_arr_on_plane->edges_begin();
         eit != temp_arr_on_plane->edges_end();
         ++eit)
    {
      std::cout << change_color(CGAL_BLUE,"Curve = ",eit->curve()) 
                << std::endl;
    }
#endif
    Face_iterator fit;
    for (fit = temp_arr_on_plane->faces_begin();
         fit != temp_arr_on_plane->faces_end();
         ++fit)
    {
       /* Mark the face as not plane face. */
      if (!fit->is_unbounded() && !is_face_at_arr_valid(fit,S3))
      {
#if ARR_ON_SUR_DEBUG        
        std::cout << change_color(CGAL_CYAN,"Set face count 0") 
                  << std::endl;
#endif
        fit->clear();
      }
    }
    arr_plane_faces.add_element(temp_arr_on_plane,temp_obs_on_plane);
  }

  /* Return true if the segment is one of the square sides */ 
  bool is_segment_a_side(const Rational_arc_2& arc)
  {
    if (m_traits_2_adapt.orientation(arc,m_alg_kernel) == CGAL::COLLINEAR)
    {
      Point_2 source(arc.source());
      Point_2 target(arc.target());
      /* 0,0 -> 0,1 */
      return ((source.x() == Rational(0) && source.y() == Rational(0) &&
               target.x() == Rational(0) && target.y() == Rational(1)) ||
              /* 0,1 -> 0,0 */
              (source.x() == Rational(0) && source.y() == Rational(1) &&
               target.x() == Rational(0) && target.y() == Rational(0)) ||
              /* 0,1 -> 1,1 */
              (source.x() == Rational(0) && source.y() == Rational(1) &&
               target.x() == Rational(1) && target.y() == Rational(1)) ||
              /* 1,1 -> 0,1 */
              (source.x() == Rational(1) && source.y() == Rational(1) &&
               target.x() == Rational(0) && target.y() == Rational(1)) ||
              /* 1,1 -> 1,0 */
              (source.x() == Rational(1) && source.y() == Rational(1) &&
               target.x() == Rational(1) && target.y() == Rational(0)) ||
              /* 1,0 -> 1,1 */
              (source.x() == Rational(1) && source.y() == Rational(0) &&
               target.x() == Rational(1) && target.y() == Rational(1)) ||
              /* 1,0 -> 0,0 */
              (source.x() == Rational(1) && source.y() == Rational(0) &&
               target.x() == Rational(0) && target.y() == Rational(0)) ||
              /* 0,0 -> 1,0 */
              (source.x() == Rational(0) && source.y() == Rational(0) &&
               target.x() == Rational(1) && target.y() == Rational(0)));
    }
    return false;
  }
  /*************************************************************
   * The following function gets a curve on the arrangement and return true 
   * if the curve is a segment on the bounded box [0,0],[1,1]
   **************************************************************/
  bool is_on_square_boudns(Halfedge_handle e)
  {
    Point_2 source(e->curve().source());
    Point_2 target(e->curve().target());
    return ((source.x() == Rational(0) && target.x() == Rational(0)) ||
            (source.x() == Rational(1) && target.x() == Rational(1)) ||
            (source.y() == Rational(0) && target.y() == Rational(0)) ||
            (source.y() == Rational(1) && target.y() == Rational(1)));
  }

  /*************************************************************
   * The following function gets 2 points on the arrangement and return true 
   * if the segment that conects these 2 points represents 3D lines 
   * that crosses the 3 lines.
   **************************************************************/
  template <typename Point_a_2, typename Point_b_2>
  bool is_line_at_arr_valid(const Point_a_2& first_p,
                            const Point_b_2& second_p,
                            const Rational_segment_3& S3)
  {
    Lines_through_segments_get_algebraic_number_adapt<Traits_3>
      get_algebraic_number_adapt;
         
    Algebraic first_p_x = get_algebraic_number_adapt(first_p.x());
    Algebraic second_p_x = get_algebraic_number_adapt(second_p.x());
    Algebraic first_p_y = get_algebraic_number_adapt(first_p.y());
    Algebraic second_p_y = get_algebraic_number_adapt(second_p.y());
         
    Alg_line_3 common_line;
    CGAL_assertion_code(int status =)
       m_g_func.get_line_from_intersection_point((first_p_x + second_p_x)/2,
                                                 (first_p_y + second_p_y)/2,
                                                 *m_S1, *m_S2, common_line);
       
    CGAL_assertion(status == LTS_g_func::CGAL_QUERY_SUCCEED);
        
    if (!m_g_func.do_intersect_line_segment(common_line, S3, m_alg_kernel))
    {
      return false;
    }
    return true;
  }
      
  /*************************************************************
   * The following function gets a face at the arrangement and return true
   * if its interior represent lines that passes through m_S1, m_S2 and S3.
   * 
   **************************************************************/
  bool is_face_at_arr_valid(const Face_handle face,
                            const Rational_segment_3& S3)
  {
    Algebraic max_y = -1;
    
    if (face->is_unbounded())
      return false;
        
    typename Arrangement_on_plane_2::Ccb_halfedge_circulator curr;
    typename Arrangement_on_plane_2::Ccb_halfedge_circulator circ =
      face->outer_ccb();
    typename Arrangement_on_plane_2::Halfedge_const_handle max_e;
         
    curr = circ;
    /* Choose non vertical edge and save it at max_e*/
    do {
      if (curr->source()->point().x() != curr->target()->point().x())
      {
        max_e = curr;
        curr = circ;
      }
      else
        ++curr;
    }
    while (curr != circ);

    /* Draw a vertical segment  from the middle of the edge.
       Iterate on all of the edges and find the intersection point of the 
       edges.
       Get the only intersection point.
    */
         
    Point_2 mid_p;
    m_traits_2_adapt.get_mid_point(max_e->curve(),
                                   mid_p);

    curr = circ = face->outer_ccb();
    max_y = -1;
    typename Arrangement_on_plane_2::Halfedge_const_handle temp_edge;
    do {
      temp_edge = curr;
      if (temp_edge != max_e && temp_edge->twin() != max_e)
      {

        if (curr->curve().left().x() <= mid_p.x() &&
            curr->curve().right().x() >= mid_p.x())
        {
          /* Get the y value of curve at x = mid_p.x() */
          Algebraic temp_y = m_traits_2_adapt.get_y_val(curr->curve(),
                                                        mid_p.x());

          if (temp_y > max_y)
          {
            max_y = temp_y;
          }
        }
      }
      curr++;
    } while (curr != circ);
         
    CGAL_assertion(max_y != Algebraic(-1));
    /* Check if the middle point between the two points represent a valid
       line that passes through the three segments.
    */       

    Alg_line_3 common_line;
    CGAL_assertion_code(int status =)
       m_g_func.get_line_from_intersection_point(mid_p.x(),
                                                 (max_y + mid_p.y())/2,
                                                 *m_S1, *m_S2, common_line);
         
    CGAL_assertion(status == LTS_g_func::CGAL_QUERY_SUCCEED);
    
    if (!m_g_func.do_intersect_line_segment(common_line, S3, m_alg_kernel))
    {
      return false;
    }
    return true;
  }
};

} //namespace CGAL

#endif /* LINE_THROUGH_SEGMENTS_3_SEGS_2_H */
