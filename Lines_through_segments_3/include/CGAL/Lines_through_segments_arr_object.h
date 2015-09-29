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

#ifndef LINE_THROUGH_SEGMENTS_ARR_OBJECT_H
#define LINE_THROUGH_SEGMENTS_ARR_OBJECT_H

/*! \file
 * This class represent 2 dimensional Rational Object of the following types:
 * 1. Hyperbola:
 *    The coefficients of the variables at the hyperbola are rational numbers.
 * 
 *       y = (x * a1 + a0) / 
 *           (x * b1 + b0)
 *
 * 2. Degenerate Hyperbola - represented by 2 line segments (horizontal and
 *    vertical).
 *
 * 3. Segment: y = x * a1 + a0.
 *
 */

#include <CGAL/Lines_through_segments_bounded_segs_vector.h>
#include <CGAL/Lines_through_segments_general_functions.h>
#include <CGAL/Lines_through_segments_arr_gen_func.h>
#include <CGAL/Lines_through_segments_traits_2_adapt.h>

namespace CGAL {

template <typename Traits_3_, 
          typename With_segments>
class Lines_through_segments_arr_object {
public:
   typedef Traits_3_                                     Traits_3;
private:
  typedef Lines_through_segments_arr_object             Self;
  
  typedef typename Traits_3::Alg_kernel       Alg_kernel;
  typedef typename Traits_3::Rational_kernel  Rational_kernel;
  typedef typename Alg_kernel::FT             Algebraic;
  typedef typename Rational_kernel::FT        Rational;
  typedef typename Rational_kernel::Point_3   Rational_point_3;
  typedef typename Rational_kernel::Line_3    Rational_line_3;
  typedef typename Rational_kernel::Segment_3 Rational_segment_3;
  typedef typename Rational_kernel::Plane_3   Rational_plane_3;
  typedef typename Rational_kernel::Point_2   Rational_point_2;
  typedef Lines_through_segments_bounded_segs_vector<Rational>
                                                       Bounded_segs_vector;

  typedef Lines_through_segments_rbound_unbound_union<Rational>
                                                        Rbound;
  
  typedef Lines_through_segments_bounded_seg<Rational>  Bounded_seg;
  typedef Lines_through_segments_rbound_unbound_union<Rational>  LTS_rbound;  

  typedef typename Traits_3::Traits_arr_on_plane_2
  Traits_arr_on_plane_2;
#define  CGAL_SEGMENT_INVALID -314
  static const unsigned int CGAL_MAX_ITERATIONS = 6;
  static const unsigned int CGAL_VALID_ITERATION = (CGAL_MAX_ITERATIONS+1);
      
private:
  Rational m_a1;
  Rational m_a0;
  Rational m_b1;
  Rational m_b0;

  /* Since segments are used, only parts of the hyperbola are used. */
  Bounded_segs_vector m_bounded_segs;
      
  const Rational_segment_3* m_creator_segment;
  const Rational_segment_3* m_S1;
  const Rational_segment_3* m_S2;
  
  enum Obj_on_arr_type {
    CGAL_NOT_DETERMINED_YET = 0,
    CGAL_HYPERBOLA = 1,
    CGAL_DEGENERATE_HYPERBOLA = 2,
    CGAL_SEGMENT = 3
  };
      
  int m_obj_on_arr_type;
  Rational m_vertical_segment_x_coordinate;
  Rational m_horizontal_segment_y_coordinate;
  Bounded_segs_vector m_vertical_segment_bounded_segs;
  Bounded_segs_vector m_horizontal_segment_bounded_segs;
  const Rational_kernel *m_rat_kernel;

  /* Used for lines through polytope when s2 is only bounded later. */
  bool m_bound_s1_s2;
  
  Lines_through_segments_traits_on_plane_adapt<Traits_3> m_traits_2_adapt;
  typedef Lines_through_segments_general_functions<Traits_3> LTS_g_func;
  typedef Lines_through_segments_arr_gen_func<Traits_3, With_segments>
    LTS_arr_g_func;

  LTS_arr_g_func m_arr_g_func;
  LTS_g_func m_g_func;
      
public:    
  Lines_through_segments_arr_object()
  {
    m_a1 = 0;m_a0 = 0;m_b1 = 0;m_b0 = 0;
    m_vertical_segment_x_coordinate = CGAL_SEGMENT_INVALID;
    m_horizontal_segment_y_coordinate = CGAL_SEGMENT_INVALID;
    m_obj_on_arr_type = CGAL_NOT_DETERMINED_YET;
    m_rat_kernel = NULL;
    m_bound_s1_s2 = true;
  }

  Lines_through_segments_arr_object(const Rational_segment_3* S1,
                                    const Rational_segment_3* S2,
                                    const Rational_segment_3* S3,
                                    const Rational_kernel* rat_kernel,
                                    bool bound_s1_s2)
  {
    m_vertical_segment_x_coordinate = CGAL_SEGMENT_INVALID;
    m_horizontal_segment_y_coordinate = CGAL_SEGMENT_INVALID;
    m_obj_on_arr_type = CGAL_NOT_DETERMINED_YET;
    m_rat_kernel = rat_kernel;
    m_bound_s1_s2 = bound_s1_s2;
    this->create_arr_object(S1,S2,S3);
  }
      
  Lines_through_segments_arr_object(const Self& _hyperbola)
  {
    m_a1 = _hyperbola.m_a1;
    m_a0 = _hyperbola.m_a0;
    m_b1 = _hyperbola.m_b1;
    m_b0 = _hyperbola.m_b0;
    m_S1 = _hyperbola.m_S1;
    m_S2 = _hyperbola.m_S2;

    /* Since segments are used, only parts of the hyperbola are used. */
    m_bounded_segs = Bounded_segs_vector(_hyperbola.m_bounded_segs);
    m_creator_segment = _hyperbola.m_creator_segment;
         
    m_vertical_segment_x_coordinate =
      _hyperbola.m_vertical_segment_x_coordinate;
    m_horizontal_segment_y_coordinate =
      _hyperbola.m_horizontal_segment_y_coordinate;
    m_obj_on_arr_type = CGAL_NOT_DETERMINED_YET;

    m_rat_kernel = _hyperbola.m_rat_kernel;
  }
            
  /*************************************************************
   * Function description:
   * -------------------- 
   * The following functions checks if the hyperbola is degenerate.
   * A hyperbola is degenrate if one of the following cases holds.
   *
   * 1. S1 and L3 intersect.
   * 2. S2 and L3 intersect.
   *
   *
   *
   *************************************************************/
  bool is_degernate_hyperbola()
  {
    return (m_rat_kernel->do_intersect_3_object()
            (m_S1->supporting_line(), m_creator_segment->supporting_line()) ||
            m_rat_kernel->do_intersect_3_object()
            (m_S2->supporting_line(), m_creator_segment->supporting_line()));
  }

  void set_degenerate_hyperbola()
  {     
    Rational_point_3 i_point;
    /* Degenerate case 2 */
    /* S1 and S3 intersecting. */
    if (m_rat_kernel->do_intersect_3_object()
        (m_S1->supporting_line(), m_creator_segment->supporting_line()))
    {
      if (m_bound_s1_s2)
      {
        if (m_g_func.do_intersect_segments(*m_S1,*m_creator_segment,
                                           i_point,m_rat_kernel))
        {
                 
          m_vertical_segment_bounded_segs.add_bounded_seg
            (Bounded_seg(Rbound(Rational(1)), Rbound(Rational(0)) ,true, true));
        }
      }
      else
      {
        if (m_g_func.do_intersect_line_segment(m_S1->supporting_line(),
                                               *m_creator_segment,
                                               m_rat_kernel))
        {
          m_vertical_segment_bounded_segs.
            add_bounded_seg(Bounded_seg(
                               LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY,
                               LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY,
                               false, false));
        }
      }
           
           
      calc_degenerate_hyperbola(m_S1,m_S2,
                                m_horizontal_segment_bounded_segs,
                                m_vertical_segment_x_coordinate,
                                m_horizontal_segment_y_coordinate,
                                m_bound_s1_s2);
            
#if OBJ_ON_ARR_DEBUG
      std::cout << "m_rat_kernel->do_intersect_3_object()(S1,CL)" 
                << std::endl;
#endif
    }
         
    /* Degenerate case 1 */
    else if (m_rat_kernel->do_intersect_3_object()
             (m_S2->supporting_line(), m_creator_segment->supporting_line()))
    {
      Rational_line_3 L2(m_S2->source(),m_S2->target());

      if (((m_bound_s1_s2 && 
           m_g_func.do_intersect_segments(*m_S2,*m_creator_segment,
                                          i_point,m_rat_kernel))) ||
          (!m_bound_s1_s2 &&
          m_g_func.do_intersect_line_segment(L2, *m_creator_segment,
                                             m_rat_kernel)))
      {
        if (m_bound_s1_s2)
        {
          m_horizontal_segment_bounded_segs.add_bounded_seg
            (Bounded_seg(Rbound(Rational(1)), Rbound(Rational(0)), true, true));
        }
        else
        {
          m_horizontal_segment_bounded_segs.add_bounded_seg
            (Bounded_seg(LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY,
                         LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY, 
                         false, false));
        }
      }
        
      calc_degenerate_hyperbola(m_S2,m_S1,
                                m_vertical_segment_bounded_segs,
                                m_horizontal_segment_y_coordinate,
                                m_vertical_segment_x_coordinate,
                                m_bound_s1_s2);
#if OBJ_ON_ARR_DEBUG
      std::cout << "m_rat_kernel->do_intersect_3_object()(S2,CL)" 
                << std::endl;
#endif
    }
#if OBJ_ON_ARR_DEBUG
    std::cout << "m_vertical_segment_x_coordinate = " 
              << m_vertical_segment_x_coordinate
              << m_vertical_segment_bounded_segs << std::endl;
    std::cout << "m_horizontal_segment_y_coordinate = " 
              << m_horizontal_segment_y_coordinate 
              << m_horizontal_segment_bounded_segs << std::endl;
#endif
  }
  
  /*************************************************************
   * The following constructor creates a hyperbola equation of  segment
   * passing through 3 Segments in 3 space.
   *
   * Input:
   *       S1         - Segment in 3 space.
   *       S2         - Segment in 3 space.
   *       L3         - Segment in 3 space.
   *
   * Output:
   *       
   * 
   * Function Description:
   * -------------------- 
   * Given 3 blue segment S1, S2, L3.
   * S1:(S1.x1 + S1.t * S1.x2, S1.y1 + S1.t * S1.y2, S1.z1 + S1.t * S1.z2)
   * S2:(S2.x1 + S2.t * S2.x2, S2.y1 + S2.t * S2.y2, S2.z1 + S2.t * S2.z2)
   * L3:(L3.x1 + L3.t * L3.x2, L3.y1 + L3.t * L3.y2, L3.z1 + L3.t * L3.z2)
   *
   * The function finds a general equation of a red segment intersecting all
   * 3 blue segments. 
   * There are infinity such red segments, all of them are defined by 
   * hyperbola.
   *
   * A segment equation of general segment which intersects with both S1 and
   * S2
   * will look like: (1-t)*I1 + t*I2,
   * Where I1, I2 are the intersection points of this segment with 
   * S1 and S2 respectivly.
   * If this segment also intersect L3, the point I3 = (1-t)*I1 + t*I2.
   * In order to get the hyperbola the following 3 equation are solved:
   *    (1-t)*(S1.x1 + S1.t * S1.x2) + t*(S2.x1 + S2.t * S2.x2) = 
   *            (L3.x1 + L3.t * L3.x2)
   *    (1-t)*(S1.y1 + S1.t * S1.y2) + t*(S2.y1 + S2.t * S2.y2) = 
   *            (L3.y1 + L3.t * L3.y2)
   *    (1-t)*(S1.z1 + S1.t * S1.z2) + t*(S2.z1 + S2.t * S2.z2) = 
   *            (L3.z1 + L3.t * L3.z2)
   *
   * From the third eqution we will get L3.t in terms of t.
   *
   *   Error check 1 - Verify that (L3.z2 != 0)
   *    If (L3.z2 == 0) We swap the equations s.t L3.z2 won't be 0.
   *                  (First we try to swap with the second equation and 
   *                   then with the first.)
   *
   *        If S1.z2 == 0 && S2.z2 == 0 && L3.z2 == 0 the segments are at 
   *        the same plane or on parallel planes.
   *
   *  L3.t = (((1-t)*(S1.z1 + S1.t * S1.z2) + 
   *            t*(S2.z1 + S2.t * S2.z2) - L3.z1)/L3.z2)
   *
   * By placing L3.t at the first and second equations we get.
   * 
   *  (1-t)*(S1.x1 + S1.t * S1.x2) + t*(S2.x1 + S2.t * S2.x2) =
   *  (L3.x1 + (((1-t)*(S1.z1 + S1.t * S1.z2) + 
   *      t * (S2.z1 + S2.t * S2.z2) - L3.z1)/L3.z2) * L3.x2)
   *    
   *  (1-t)*(S1.y1 + S1.t * S1.y2) + t*(S2.y1 + S2.t * S2.y2) =
   *  (L3.y1 + (((1-t)*(S1.z1 + S1.t * S1.z2) + t*(S2.z1 + S2.t * S2.z2) - 
   *      L3.z1)/L3.z2) * L3.y2)
   *
   *
   * Get S1.t from the first equation:
   * ----------------------------------
   *
   *  (1-t)*(S1.x1 + S1.t * S1.x2) + t*(S2.x1 + S2.t * S2.x2)=
   *  (L3.x1 + (((1-t)*(S1.z1 + S1.t * S1.z2) + t*(S2.z1 + S2.t * S2.z2) - 
   *      L3.z1)/L3.z2) * L3.x2)
   *
   * Multiply by L3.z2:
   *  L3.z2 * (1-t) * (S1.x1 + S1.t * S1.x2) + t*(S2.x1 + S2.t * S2.x2)=
   *(L3.z2 * L3.x1 + ((1-t)*(S1.z1 + S1.t * S1.z2) + t*(S2.z1 + S2.t * S2.z2)
   *      - L3.z1) * L3.x2)
   *
   * Open the first parentheses multiply by 1-t:
   *L3.z2 * (S1.x1 + S1.t * S1.x2 - t * S1.x1 - t * S1.t * S1.x2 + t * (S2.x1
   *     + S2.t * S2.x2)) = 
   *  L3.x1 * L3.z2 + ((1-t)*(S1.z1 + S1.t * S1.z2) + t * 
   *      (S2.z1 + S2.t * S2.z2) - L3.z1) * L3.x2
   *
   * Open the second parentheses multiply by 1-t:
   *  L3.z2 * (S1.x1 + S1.t * S1.x2 - t * S1.x1 - t * S1.t * S1.x2 + t * 
   *     (S2.x1 + S2.t * S2.x2)) = 
   *  L3.x1 * L3.z2 + L3.x2 * (S1.z1 + S1.t * S1.z2 - t * S1.z1 - 
   *      t * S1.t * S1.z2 + t * S2.z1 + t * S2.t * S2.z2 - L3.z1)
   *
   *  S1.t * (L3.z2 * S1.x2 - L3.z2 * t * S1.x2) +
   *  (L3.z2 * S1.x1 - L3.z2 * t * S1.x1 + L3.z2 * t * S2.x1 + 
   *     L3.z2 * t * S2.t * S2.x2) = 
   *  S1.t * (L3.x2 * S1.z2  - L3.x2 * t * S1.z2) + L3.x1 * L3.z2 + 
   *     L3.x2 * S1.z1  - L3.x2 * t * S1.z1 +
   *  L3.x2 * t * S2.z1 + L3.x2 * t * S2.t * S2.z2 - L3.x2 * L3.z1
   *
   * Isolate S1.t - get it out of the parentheses:
   *  S1.t * (L3.z2 * S1.x2 - t * L3.z2 * S1.x2) + 
   *       L3.z2 * (S1.x1 - t * S1.x1 + t*(S2.x1 + S2.t * S2.x2)) =
   *  L3.x1 * L3.z2 + S1.t * (L3.x2 * S1.z2 - t * L3.x2 * S1.z2) +
   *  L3.x2 * (S1.z1 - t * S1.z1  + t * S2.z1 + t * S2.t * S2.z2 - L3.z1)
   *
   * Move S1.t to the left side of the equation:
   *  S1.t * ((L3.z2 * S1.x2 - t * L3.z2 * S1.x2) - 
   *        (L3.x2 * S1.z2 - t * L3.x2 * S1.z2)) =
   *  (L3.x1 * L3.z2 + L3.x2 * 
   *    (S1.z1 - t * S1.z1  + t * S2.z1 + t * S2.t * S2.z2 - L3.z1)) - 
   *    L3.z2 * (S1.x1 - t * S1.x1 + t*(S2.x1 + S2.t * S2.x2))
   *
   *  S1.t * (L3.z2 * S1.x2 - t * L3.z2 * S1.x2 - L3.x2 * S1.z2  + 
   *      t * L3.x2 * S1.z2) =
   * L3.x1 * L3.z2 + L3.x2 * S1.z1  - L3.x2 * t * S1.z1 + L3.x2 * t * S2.z1 +
   *  L3.x2 * t * S2.t * S2.z2 - L3.x2 * L3.z1 - L3.z2 * S1.x1 + 
   *      L3.z2 * t * S1.x1 - L3.z2 * t * S2.x1 - L3.z2 * t * S2.t * S2.x2)
   *
   *
   *  S1.t =
   *  ((L3.x1 * L3.z2 + L3.x2 * 
   *         (S1.z1 - t * S1.z1  + t * S2.z1 + t * S2.t * S2.z2 - L3.z1) - 
   *  (L3.z2 * (S1.x1 - t * S1.x1 + t*(S2.x1 + S2.t * S2.x2)))) / 
   *  ((L3.z2 * S1.x2 - t* L3.z2 * S1.x2) - (L3.x2 * S1.z2 - t* L3.x2 * 
   *  S1.z2)))
   *
   *  S1.t  =
   *  (L3.x1 * L3.z2 + L3.x2 * S1.z1  - t * L3.x2 * S1.z1 + t * 
   *       L3.x2 * S2.z1 + t * L3.x2 * S2.t * S2.z2 - L3.x2 * L3.z1 -
   *  L3.z2 * S1.x1 + t * L3.z2 * S1.x1 - t * L3.z2 * S2.x1 - 
   *         t * L3.z2 * S2.t * S2.x2) /
   * (L3.z2 * S1.x2 - t * L3.z2 * S1.x2 - L3.x2 * S1.z2  + t * L3.x2 * S1.z2)
   *
   *  Let C1 = L3.z2 * S1.x2
   *  Let C2 = L3.x2 * S1.z2
   *  Let C3 = C1 - C2
   *  (C3 - t * C3)
   * 
   *  S1.t =
   *     ((L3.x1 * L3.z2 + L3.x2 * (S1.z1 - t * S1.z1  + t * S2.z1 + 
   *     t * S2.t * S2.z2 - L3.z1) - 
   *     (L3.z2 * (S1.x1 - t * S1.x1 + t*(S2.x1 + S2.t * S2.x2)))) /
   *     (C3 * (1 - t)))
   *
   *  S1.t  =
   *     (L3.x1 * L3.z2 + L3.x2 * S1.z1  - 
   *     t * L3.x2 * S1.z1 + t * L3.x2 * S2.z1 +
   *     t * L3.x2 * S2.t * S2.z2 - L3.x2 * L3.z1 - 
   *     L3.z2 * S1.x1 + t * L3.z2 * S1.x1 - t * L3.z2 * S2.x1 -
   *     t * L3.z2 * S2.t * S2.x2) / (C3 - t * C3)
   *
   * Sort the variables into 2 groups scalars and variables with 
   * one degree (of t).
   *
   *  S1.t  =
   *     (L3.x1 * L3.z2 + L3.x2 * S1.z1  - L3.x2 * L3.z1 - L3.z2 * S1.x1 + 
   *     t * L3.z2 * S1.x1 -
   *     t * L3.z2 * S2.x1 - t * L3.z2 * S2.t * S2.x2 - t * L3.x2 * S1.z1 + 
   *     t * L3.x2 * S2.z1 + t * L3.x2 * S2.t * S2.z2) / (C3 - t * C3)
   *
   *  Let C4 = L3.x1 * L3.z2 + L3.x2 * S1.z1 - L3.x2 * L3.z1 - L3.z2 * S1.x1
   *
   *  S1.t  =
   *     (C4 + t * L3.z2 * S1.x1 - t * L3.z2 * S2.x1 - 
   *     t * L3.z2 * S2.t * S2.x2 - t * L3.x2 * S1.z1 +
   *     t * L3.x2 * S2.z1 + t * L3.x2 * S2.t * S2.z2) / (C3 - t * C3)
   *
   *  S1.t  =
   *   (C4 + t * S2.t * (L3.x2 * S2.z2 - L3.z2 * S2.x2) + t * L3.z2 * S1.x1 -
   *     t * L3.z2 * S2.x1 -
   *     t * L3.x2 * S1.z1 + t * L3.x2 * S2.z1) / (C3 - t * C3)
   *
   *  Let C5 = L3.x2 * S2.z2 - L3.z2 * S2.x2
   *
   *  S1.t  =
   *     (C4 + t * S2.t * C5 + t * 
   *     (L3.z2 * S1.x1 - L3.z2 * S2.x1 - L3.x2 * S1.z1 + L3.x2 * S2.z1) / 
   *     (C3 - t * C3)
   * 
   *  S1.t =
   *     (C4 + S2.t * t * C5 - t * 
   *     (- L3.x2 * S1.z1 + L3.x2 * S2.z1 + L3.z2 * S1.x1 - L3.z2 * S2.x1) /
   *     (C3 * (1 - t))
   *     
   * Let C6 = - L3.x2 * S1.z1 + L3.x2 * S2.z1 + L3.z2 * S1.x1 - L3.z2 * S2.x1
   *
   *  S1.t = (C4 + S2.t * t * C5 + t * C6) / (C3 * (1 - t)) 
   *
   *
   *
   * Get S2.t from the second equation:
   * ----------------------------------
   *
   *  (1-t)*(S1.y1 + S1.t * S1.y2) + t*(S2.y1 + S2.t * S2.y2) =
   *  (L3.y1 + (((1-t)*(S1.z1 + S1.t * S1.z2) + 
   *  t*(S2.z1 + S2.t * S2.z2) - L3.z1)/L3.z2) * L3.y2)
   *
   * Multiply by L3.z2.
   *
   *  L3.z2 * (1-t)*(S1.y1 + S1.t * S1.y2) + 
   *   L3.z2 * t * (S2.y1 + S2.t * S2.y2) =
   *  (L3.y1 * L3.z2 + (((1-t)*(S1.z1 + S1.t * S1.z2) + 
   *    t*(S2.z1 + S2.t * S2.z2) - L3.z1)) * L3.y2)
   *     
   * Open the first parentheses multiply by 1-t: 
   *
   *  L3.z2 * (S1.y1 + S1.t * S1.y2 - t * S1.y1 - 
   *    t * S1.t * S1.y2 + t * (S2.y1 + S2.t * S2.y2)) = 
   *  L3.y1 * L3.z2 + (((1-t)*(S1.z1 + S1.t * S1.z2) + 
   *    t*(S2.z1 + S2.t * S2.z2) - L3.z1)) * L3.y2
   *
   *  L3.z2 * S1.y1 + L3.z2 * S1.t * S1.y2 - L3.z2 * t * S1.y1 - 
   *    L3.z2 * t * S1.t * S1.y2 +
   *  L3.z2 * t * S2.y1 + L3.z2 * t * S2.t * S2.y2 =
   *  L3.y1 * L3.z2 + S1.z1 * L3.y2  + S1.t * S1.z2 * L3.y2 - 
   *  t * S1.z1 * L3.y2 - t * S1.t * S1.z2 * L3.y2  + 
   *  t * S2.z1 * L3.y2 + t * S2.t * S2.z2 * L3.y2  - L3.z1 * L3.y2
   *     
   * Open the second parentheses multiply by 1-t:
   *
   *  (L3.z2 * S1.y1 + L3.z2 * S1.t * S1.y2 - L3.z2 * t * S1.y1 - 
   *  L3.z2 * t * S1.t * S1.y2 +
   *  L3.z2 * t * S2.y1 + L3.z2 * t * S2.t * S2.y2) = 
   *  L3.y1 * L3.z2 + S1.z1 * L3.y2 + S1.t * S1.z2 * L3.y2 - 
   *  t * S1.z1 * L3.y2 - t * S1.t * S1.z2 * L3.y2 + 
   *  t * S2.z1 * L3.y2 + t * S2.t * S2.z2 * L3.y2 - L3.z1 * L3.y2
   *  
   * Isolate S2.t at the left side of the equation:
   *
   *  (L3.z2 * S1.y1 + L3.z2 * S1.t * S1.y2 - L3.z2 * t * S1.y1 - 
   *  L3.z2 * t * S1.t * S1.y2 + L3.z2 * t * S2.y1 + 
   *  L3.z2 * t * S2.t * S2.y2) = 
   *  L3.y1 * L3.z2 + S1.z1 * L3.y2 + S1.t * S1.z2 * L3.y2 - 
   *  t * S1.z1 * L3.y2 - t * S1.t * S1.z2 * L3.y2 +
   *  t * S2.z1 * L3.y2 + t * S2.t * S2.z2 * L3.y2 - L3.z1 * L3.y2
   *
   *  (L3.z2 * S1.y1 + L3.z2 * S1.t * S1.y2 - L3.z2 * t * S1.y1 - 
   *  L3.z2 * t * S1.t * S1.y2 + L3.z2 * t * S2.y1) = 
   *  L3.y1 * L3.z2 + S1.z1 * L3.y2 + S1.t * S1.z2 * L3.y2 - 
   *  t * S1.z1 * L3.y2 -
   *  t * S1.t * S1.z2 * L3.y2 + t * S2.z1 * L3.y2 - L3.z1 * L3.y2
   *
   *  S2.t * t =
   *     (L3.y1 * L3.z2 + S1.z1 * L3.y2  + S1.t * S1.z2 * L3.y2 - 
   *     t * S1.z1 * L3.y2 - t * S1.t * S1.z2 * L3.y2  +
   *     t * S2.z1 * L3.y2   - L3.z1 * L3.y2 - L3.z2 * S1.y1 - 
   *     L3.z2 * S1.t * S1.y2 +
   *     L3.z2 * t * S1.y1 + L3.z2 * t * S1.t * S1.y2 - L3.z2 * t * S2.y1) /
   *     (L3.z2 * S2.y2 - S2.z2 * L3.y2)
   *
   *
   *  S2.t * t =
   *   (L3.y1 * L3.z2 + S1.z1 * L3.y2 - L3.z2 * S1.y1 - L3.z1 * L3.y2
   *   - L3.z2 * S1.t * S1.y2 + L3.z2 * t * S1.y1 + L3.z2 * t * S1.t * S1.y2 
   *     - L3.z2 * t * S2.y1 +
   *   S1.t * S1.z2 * L3.y2 - t * S1.z1 * L3.y2 - t * S1.t * S1.z2 * L3.y2  +
   *     t * S2.z1 * L3.y2) /
   *     (L3.z2 * S2.y2 - S2.z2 * L3.y2)
   *
   *  Let C101 = (L3.z2 * S2.y2  - S2.z2 * L3.y2)
   *
   * Sort the variables to 2 groups scalars and variables with one degree:
   *
   * Let C102 = L3.y1 * L3.z2 + S1.z1 * L3.y2 - L3.z2 * S1.y1 - L3.z1 * L3.y2
   *
   *  S2.t * t =
   *     (C102 - L3.z2 * S1.t * S1.y2 + L3.z2 * t * S1.y1 + 
   *     L3.z2 * t * S1.t * S1.y2 -
   *     L3.z2 * t * S2.y1 + S1.t * S1.z2 * L3.y2 - t * S1.z1 * L3.y2 - 
   *     t * S1.t * S1.z2 * L3.y2  + t * S2.z1 * L3.y2) / C101
   *
   *  S2.t =
   *     (C102 - (L3.z2 * S1.t * S1.y2 - L3.z2 * t * S1.y1 - 
   *     L3.z2 * t * S1.t * S1.y2 + L3.z2 * t * S2.y1) +
   *     ( S1.t * S1.z2 * L3.y2 - t * S1.z1 * L3.y2 - 
   *     t * S1.t * S1.z2 * L3.y2 + t * S2.z1 * L3.y2)) / (t * C101))
   *     
   *  S2.t = 
   *     (C102 - (L3.z2 * S1.t * S1.y2 - L3.z2 * t * S1.y1 - 
   *     L3.z2 * t * S1.t * S1.y2 + L3.z2 * t * S2.y1) +
   *     ( S1.t * S1.z2 * L3.y2 - t * S1.z1 * L3.y2 - 
   *     t * S1.t * S1.z2 * L3.y2 + t * S2.z1 * L3.y2)) / (t * C101)
   *     
   *  S2.t = 
   *     (C102 + S1.t * (S1.z2 * L3.y2 - L3.z2 * S1.y2 + L3.z2 * t * S1.y2 - 
   *     t * S1.z2 * L3.y2) - ( - L3.z2 * t * S1.y1 + L3.z2 * t * S2.y1) +
   *     (t * S1.z1 * L3.y2 + t * S2.z1 * L3.y2)) / 
   *     (t * C101)
   *
   *  S2.t * t =
   *     (C102 + S1.t * (t * L3.z2 * S1.y2 - L3.z2 * S1.y2 + 
   *     S1.z2 * L3.y2 - t * S1.z2 * L3.y2) + t * (L3.z2 * S1.y1 - 
   *     L3.z2 * S2.y1 - S1.z1 * L3.y2 + S2.z1 * L3.y2)) / C101
   *
   *  
   *  Let C103 = L3.z2 * S1.y2
   *  Let C104 = S1.z2 * L3.y2
   *    
   *  S2.t * t =
   *     (C102 + S1.t * (t * C103 - C103 + C104 - t * C104) +
   *     t * (L3.z2 * S1.y1 - L3.z2 * S2.y1 - S1.z1 * L3.y2 + S2.z1 * L3.y2))
   *     / C101
   *
   *  Let C105 = (C104 - C103 + C103 * t - t * C104)
   *  Let C105 = (C104 - C103) - t * (C104 - C103)
   *  Let C105 = (C104 - C103) * (1-t)
   *
   *  S2.t * t =
   *     (C102 + (C104 - C103) * (1-t) * S1.t + t * 
   *     (L3.z2 * S1.y1 - S1.z1 * L3.y2 + S2.z1 * L3.y2 - L3.z2 * S2.y1))
   *     / C101
   *
   *  S2.t = 
   *     (C102 + (C104 - C103) * (1-t) * S1.t - 
   *     ( - L3.z2 * t * S1.y1 + L3.z2 * t * S2.y1) +
   *     ( - t * S1.z1 * L3.y2  + t * S2.z1 * L3.y2)) / 
   *     (t * C101)
   *
   *Let C106 =  L3.z2 * S1.y1 - S1.z1 * L3.y2 + S2.z1 * L3.y2 - L3.z2 * S2.y1
   *
   *  S2.t * t =
   *     (C102 + (C104 - C103) * (1-t) * S1.t + t * C106) / C101
   *      
   *  S2.t = (C102 + (C104 - C103) * (1-t) * S1.t + t * C106) /(t * C101)
   *
   *  Place the value we got from the second equation 
   *    for S2.t at the first equation:
   * -----------------------------------------------------------------------
   *      
   *  S1.t =(C4 + S2.t * t * C5 + t * C6) / (C3 * (1 - t))
   S1.t =(C4 + C5 *((C102 + (C104 - C103) * (1-t) * S1.t + t * C106) / C101)
   *     + t * C6) / (C3 * (1 - t))
   *  
   *  Multiply by C101:
   *
   *  S1.t * C101 = (C4 * C101 + C5 *((C102 + (C104 - C103) * (1-t) * S1.t + 
   *     t * C106)) + t * C6 * C101) / (C3 * (1 - t))
   *  
   * Multiply by C3 :
   *
   *  S1.t * C101 * C3 = (C4 * C101 + C5 * C102 + (C5 * C104 - 
   *   C5 * C103) * (1-t) * S1.t + C5 * t * C106) + t * C6 * C101) / (1 - t)
   *
   * Isolate S1 at the left side of the equation.
   *
   *  S1.t * C101 * C3 - (C5 * C104 - C5 * C103) * S1.t = 
   *     (C4 * C101 + C5 * C102 + C5 * t * C106 + t * C6 * C101) / (1 - t)
   *
   *  S1.t * (C101 * C3 - C5 * C104 + C5 * C103) = 
   *     (C4 * C101 + C5 * C102 + C5 * t * C106 + t * C6 * C101) / (1 - t)
   *
   *  Let C202 = C101 * C3 - C5 * C104 + C5 * C103
   *
   *  S1.t * C202 = 
   *       (C4 * C101 + C5 * C102 +  t * (C5 * C106 + C6 * C101)) / (1 - t)
   *   
   *  Let C201 = C5 * C106 + C6 * C101
   *  Let C203 = C4 * C101 + C5 * C102 
   *
   *  S1.t = (C203 +  t * C201) / ((1 - t) * C202)
   *
   * Place the value we got from the first equation for S1.t at the 
   * second equation
   * ------------------------------------------------------------------------
   *  
   *  S1.t * (1 - t) = (C203 +  t * C201) /  C202
   *
   *  S2.t = (C102 + (C104 - C103) * (1-t) * S1.t + t * C106) /(t * C101)
   *
   S2.t = (C102 + (C104 - C103) * ((C203 +  t * C201) /  C202) + t * C106) /
   *      (t * C101)
   * 
   *  S2.t = (C102 * C202 + (C104 - C103) * ((C203 +  t * C201)) + 
   *     C202 * t * C106) / (t * C101 * C202)
   *
   *  S2.t = (C102 * C202 +  (C104 - C103) * C203 + 
   *       (C104 - C103) * t * C201 + C202 * t * C106) / (t * C101 * C202)
   *
   *  S2.t = 
   *   (C102 * C202 +  C104 * C203 - C103 * C203  +  
   * C104 * t * C201 - C103 * t * C201 + C202 * t * C106) / (t * C101 * C202)
   *
   *  Let C301 = C102 * C202 +  C104 * C203 - C103 * C203
   *  
   *  S2.t = (C301  +  t * (C104 * C201 - C103 * C201 + C202 * C106)) /
   *     (t * C101 * C202)
   *
   *  Let C302 =  C104 * C201 - C103 * C201 + C202 * C106
   *
   *  Let C303 = C101 * C202
   *
   *  S2.t = (C301  +  t * C302) / (t * C303)
   *
   * Get the hyperbola equation of L3 at the space defined by S1.t and S2.t:
   * Isolate t at the first equation:
   *
   *  S1.t = (C203 +  t * C201) / ((1 - t) * C202)
   *
   *  ((1 - t) * C202) * S1.t =  (C203 +  t * C201)
   *
   *  C202 * S1.t - t * S1.t * C202 = C203 +  t * C201
   *
   *  C202 * S1.t - C203 = t * C201 + t * S1.t * C202
   * 
   *  C202 * S1.t - C203 = t * (C201 + S1.t * C202)
   *
   *
   *  t = ((C202 * S1.t - C203) / (C201 + S1.t * C202))
   *
   * Place t at the second equation.
   *
   *  S2.t = (C301  +  t * C302) / (t * C303)
   *
   *  S2.t = (C301 / t  + C302) / C303
   *
   *  S2.t = 
   *      (C301 * (C201 + S1.t * C202) / (C202 * S1.t - C203)  + C302) / C303
   *
   *  S2.t = (C301 * (C201 + S1.t * C202) + C302 * (C202 * S1.t - C203)) /
   *      (C303 * (C202 * S1.t - C203))
   *
   *  S2.t = (C301 * C201 + C301 * S1.t * C202 + C302 * C202 * S1.t - 
   *     C302 * C203) / ( C303 * C202 * S1.t - C303 * C203)
   *
   *  S2.t = (C301 * C201 + C301 * S1.t * C202 + C302 * C202 * S1.t - 
   *     C302 * C203) / ( C303 * C202 * S1.t - C303 * C203)
   *
   *  Let C402 = C303 * C202
   *  Let C403 = C303 * C203
   *
   *  S2.t = (C301 * C201 - C203 * C302 + S1.t * (C202 * C302 + C301 * C202))
   *  / (C402 * S1.t - C403)
   *
   *  Let C400 = C301 * C201 - C203 * C302
   *  Let C401 = C202 * C302 + C301 * C202
   *
   *  S2.t = (C400 + S1.t * C401) / 
   *         (C402 * S1.t - C403)
   *
   *   y = (x * a1 + a0) / 
   *       (x * b1 + b0)
   *
   *
   *************************************************************/
  void create_arr_object(const Rational_segment_3* S1,
                         const Rational_segment_3* S2,
                         const Rational_segment_3* S3)
  {
    m_creator_segment = S3;
    m_S1 = S1;
    m_S2 = S2;
    
    if (is_degernate_hyperbola()){
      m_obj_on_arr_type = CGAL_DEGENERATE_HYPERBOLA;
      set_degenerate_hyperbola();
    } else {
      
      
      
      /* In order to get presentation of a segment L in a form of 
       * (x1,y1,z1) + t * (x2,y2,z2),
       * We will get two arbitrary points on L and set (x1,y1,z1) 
       * coordinates to be the first point and (x2,y2,z2) to be (x1,y1,z1)
       *  minus the coordinates of the second point. 
       */
      
      Rational_point_3 S1_P1 = m_S1->source();
      Rational_point_3 S1_P2 = m_S1->target();
            
      Rational_point_3 S2_P1 = m_S2->source();
      Rational_point_3 S2_P2 = m_S2->target();
            
      Rational_point_3 S3_P1 = m_creator_segment->source();
      Rational_point_3 S3_P2 = m_creator_segment->target();
            
      Rational S1_px = S1_P1.x();
      Rational S1_py = S1_P1.y();
      Rational S1_pz = S1_P1.z();      
      Rational S1_dx = S1_P2.x() - S1_px;
      Rational S1_dy = S1_P2.y() - S1_py;
      Rational S1_dz = S1_P2.z() - S1_pz;

      Rational S2_px = S2_P1.x();
      Rational S2_py = S2_P1.y();
      Rational S2_pz = S2_P1.z();      
      Rational S2_dx = S2_P2.x() - S2_px;
      Rational S2_dy = S2_P2.y() - S2_py;
      Rational S2_dz = S2_P2.z() - S2_pz;

      Rational S3_px = S3_P1.x();
      Rational S3_py = S3_P1.y();
      Rational S3_pz = S3_P1.z();      
      Rational S3_dx = S3_P2.x() - S3_px;
      Rational S3_dy = S3_P2.y() - S3_py;
      Rational S3_dz = S3_P2.z() - S3_pz;      

      // This handles the generic case. 
      // Three skew lines that have linear independent directions. 
      // The code is based on maple's C package 
      if(CGAL::orientation(
             S1->source()-S1->target(),
             S2->source()-S2->target(),
             S3->source()-S3->target())
          != CGAL::COPLANAR ){
     
        Rational t1 = S1_dz * S3_dx;
        Rational t3 = S2_px * S3_dy;
        Rational t5 = S1_dy * S3_dx;
        Rational t7 = S1_dx * S2_pz;
        Rational t9 = S1_dx * S3_pz;
        Rational t11 = S1_dx * S3_dz;
        Rational t14 = S1_dy * S3_px;
        Rational t16 = S1_dy * S2_px;
        Rational t19 = S3_px * S1_dz;
        Rational t22 = t1 * S2_py - t3 * S1_dz - t5 * S2_pz + t7 * S3_dy - t9 * S3_dy - t11 * S2_py + t5 * S3_pz - t14 * S3_dz + t16 * S3_dz - t1 * S3_py + t19 * S3_dy + t11 * S3_py;
        Rational t29 = S3_dx * S2_pz;
        Rational t35 = S3_dx * S3_py;
        Rational t40 = S2_px * S3_pz;
        Rational t42 = S3_px * S2_pz;
        Rational t46 = S1_py * S3_px;
        Rational t50 = S2_px * S3_dz;
        Rational t52 = S1_py * S2_px;
        Rational t54 = S1_px * S3_dz;
        Rational t57 = S1_px * S3_dy * S2_pz - t3 * S1_pz + S3_dx * S1_pz * S2_py - t29 * S1_py + S3_dx * S1_py * S3_pz - S3_dx * S2_py * S3_pz + t35 * S2_pz - t35 * S1_pz - S1_px * S3_pz * S3_dy + t40 * S3_dy - t42 * S3_dy + S3_px * S1_pz * S3_dy - t46 * S3_dz + S3_px * S2_py * S3_dz - t50 * S3_py + t52 * S3_dz - t54 * S2_py + t54 * S3_py;
        Rational t59 = S1_dx * S2_dz;
        Rational t60 = t59 * S3_dy;
        Rational t61 = S3_dx * S2_dz;
        Rational t62 = t61 * S1_dy;
        Rational t64 = S1_dx * S2_dy * S3_dz;
        Rational t65 = S1_dz * S2_dx;
        Rational t66 = t65 * S3_dy;
        Rational t67 = S2_dx * S3_dz;
        Rational t68 = t67 * S1_dy;
        Rational t69 = S3_dx * S2_dy;
        Rational t70 = t69 * S1_dz;
        Rational t71 = -t60 + t62 + t64 + t66 - t68 - t70;
        Rational t73 = t61 * S1_py;
        Rational t75 = S3_px * S2_dz;
        Rational t77 = S1_px * S2_dz;
        Rational t78 = t77 * S3_dy;
        Rational t79 = t69 * S1_pz;
        Rational t80 = t67 * S1_py;
        Rational t82 = S2_dx * S3_pz;
        Rational t84 = S1_px * S2_dy;
        Rational t85 = t84 * S3_dz;
        Rational t87 = S3_px * S2_dy;
        Rational t89 = S2_dx * S1_pz;
        Rational t90 = t89 * S3_dy;
        Rational t91 = t73 - t61 * S3_py + t75 * S3_dy - t78 - t79 - t80 + t67 * S3_py - t82 * S3_dy + t85 + t69 * S3_pz - t87 * S3_dz + t90;
        Rational t93 = t84 * S3_pz;
        Rational t94 = t87 * S1_pz;
        Rational t95 = t82 * S1_py;
        Rational t96 = t89 * S3_py;
        Rational t98 = S2_dx * S1_py * S2_pz;
        Rational t99 = S2_dx * S2_py;
        Rational t100 = t99 * S1_pz;
        Rational t101 = t99 * S3_pz;
        Rational t103 = S2_dx * S3_py * S2_pz;
        Rational t105 = S1_px * S2_pz * S2_dy;
        Rational t107 = S2_px * S1_pz * S2_dy;
        Rational t108 = t40 * S2_dy;
        Rational t109 = t42 * S2_dy;
        Rational t110 = t75 * S2_py;
        Rational t111 = t46 * S2_dz;
        Rational t112 = t52 * S2_dz;
        Rational t114 = S2_px * S3_py * S2_dz;
        Rational t115 = t77 * S2_py;
        Rational t116 = t77 * S3_py;
        Rational t117 = t93 - t94 - t95 + t96 + t98 - t100 + t101 - t103 - t105 + t107 - t108 + t109 - t110 + t111 - t112 + t114 + t115 - t116;
        Rational t119 = t59 * S2_py;
        Rational t120 = t59 * S3_py;
        Rational t121 = S1_dy * S2_dx;
        Rational t122 = t121 * S2_pz;
        Rational t123 = t121 * S3_pz;
        Rational t124 = t14 * S2_dz;
        Rational t125 = t9 * S2_dy;
        Rational t126 = t16 * S2_dz;
        Rational t127 = t7 * S2_dy;
        Rational t128 = t65 * S2_py;
        Rational t129 = t65 * S3_py;
        Rational t130 = t19 * S2_dy;
        Rational t132 = S2_px * S1_dz * S2_dy;
        Rational t133 = -t119 + t120 - t122 + t123 - t124 - t125 + t126 + t127 + t128 - t129 + t130 - t132;
        Rational t135 = t50 * S2_dy;
        Rational t136 = t29 * S2_dy;
        Rational t137 = -t135 + t136 + t73 - t78 - t79 - t80 + t85 + t90 + t93 - t94 - t95 + t96 + t98 - t100 + t101;
        Rational t139 = S2_dx * S3_dy * S2_pz;
        Rational t140 = t99 * S3_dz;
        Rational t141 = t61 * S2_py;
        Rational t142 = t3 * S2_dz;
        Rational t143 = -t103 - t105 + t107 - t108 + t109 - t110 + t111 - t112 + t114 + t115 - t116 - t139 + t140 - t141 + t142;
        Rational t146 = -t119 + t120 - t66 - t122 + t123 + t70 - t124 - t125 + t126 + t127 + t60 - t64 + t68 - t62 + t128 - t129 + t130 - t132;
        Rational t148 = t80 - t85 - t90 + t79 - t140 + t139 + t141 - t73 - t136 + t135 - t142 + t78;
        
        m_a1 = t22;
        m_a0 = t57; 
        m_b1 = t71; 
        m_b0 = t91; 
        
        CGAL_precondition(!CGAL::is_zero(t133));
        CGAL_precondition(!CGAL::is_zero(t146));
        CGAL_precondition(!CGAL::is_zero(t71));
        
        Rational b0 = t117 / t133;
        Rational b1 = (t137 + t143) / t146;
        Rational L  = t148 / t71; 
        Rational A  = -  m_b0 / m_b1; 
        
        Rational bmin = std::min(b0,b1);
        Rational bmax = std::max(b0,b1);
               
        
        static const Rational& zero = Rational(0);
        static const Rational& one  = Rational(1);

        if(bmin < L && L < bmax){
          bmin = std::min(bmin,one );
          bmax = std::max(bmax,zero);          
          // return (]-oo,bmin] \cup [bmax,+oo[) \cap [0,1]
          if( bmin >= zero ) this->add_bounded_seg(zero,bmin,A);       
          if( bmax <= one  ) this->add_bounded_seg(bmax,one ,A);      
        }else{
          bmin = std::max(bmin,zero);
          bmax = std::min(bmax,one);
          if(bmin <= bmax) this->add_bounded_seg(bmin,bmax,A);    
        }
        return;
      }


      /* If t is constant the calculation is dependent on S3_t. */
      if ((S1_dx == 0 &&
           S2_dx == 0 &&
           S3_dx == 0))
      {
        /* (1-t) * S1_px + t * S2_px = S3_px */
        CGAL_assertion(S2_px != S1_px);
        Rational t = (S3_px-S1_px)/(S2_px-S1_px);
#if OBJ_ON_ARR_DEBUG
        std::cout << "t1 = " << t << std::endl;
#endif
        Rational C_y((S3_py - (1 - t) * S1_py - t * S2_py));
        Rational C_z((S3_pz - (1 - t) * S1_pz - t * S2_pz));
        Rational S1_y_temp(S1_dy * (1-t));
        Rational S2_y_temp(S2_dy * t);
        Rational S1_z_temp(S1_dz * (1-t));
        Rational S2_z_temp(S2_dz * t);

        calc_obj_on_arr_t_is_constant(S1_y_temp,
                                      S2_y_temp,
                                      S3_dy,
                                      C_y,
                                      S1_z_temp,
                                      S2_z_temp,
                                      S3_dz,
                                      C_z);
      }
      else if (S1_dy == 0 &&
               S2_dy == 0 &&
               S3_dy == 0)
      {
        /* (1-t) * S1_py + t * S2_py = S3_py */
        CGAL_assertion(S2_py != S1_py);
        Rational t = (S3_py - S1_py)/(S2_py - S1_py);
#if OBJ_ON_ARR_DEBUG
        std::cout << "t2 = " << t << std::endl;
#endif
        Rational C_x((S3_px - (1 - t) * S1_px - t * S2_px));
        Rational C_z((S3_pz - (1 - t) * S1_pz - t * S2_pz));
        Rational S1_x_temp(S1_dx * (1-t));
        Rational S2_x_temp(S2_dx * t);
        Rational S1_z_temp(S1_dz * (1-t));
        Rational S2_z_temp(S2_dz * t);
               
        calc_obj_on_arr_t_is_constant(S1_x_temp,
                                      S2_x_temp,
                                      S3_dx,
                                      C_x,
                                      S1_z_temp,
                                      S2_z_temp,
                                      S3_dz,
                                      C_z);
      }
      else if (S1_dz == 0 &&
               S2_dz == 0 &&
               S3_dz == 0)
      {
        /* (1-t) * S1_pz + t * S2_pz = S3_pz */
        CGAL_assertion(S2_pz != S1_pz);
         
        Rational t = (S3_pz - S1_pz)/(S2_pz - S1_pz);
#if OBJ_ON_ARR_DEBUG
        std::cout << "t3 = " << t << std::endl;
#endif
        Rational C_x(-(S3_px - (1 - t) * S1_px - t * S2_px));
        Rational C_y((S3_py - (1 - t) * S1_py - t * S2_py));
        Rational S1_x_temp(S1_dx * (1-t));
        Rational S2_x_temp(S2_dx * t);
        Rational S1_y_temp(S1_dy * (1-t));
        Rational S2_y_temp(S2_dy * t);

        calc_obj_on_arr_t_is_constant(S1_x_temp,
                                      S2_x_temp,
                                      S3_dx,
                                      C_x,
                                      S1_y_temp,
                                      S2_y_temp,
                                      S3_dy,
                                      C_y);
      }
      else
      {
        Rational C1;
        Rational C2;
        Rational C3;
        Rational C4;
        Rational C5;
        Rational C6;
               
        Rational C101;
        Rational C102;
        Rational C103;
        Rational C104;
        Rational C106;
               
        Rational C201;
        Rational C202;
            
        /* Swap the equations until(S3_dz != 0 && C202 != 0 && C101 != 0) 
           run the following 6 iterations:
           (0,1,2)
           (0,2,1)
           (2,0,1)
           (2,1,0)
           (1,2,0)
           (1,0,2)
        */
        unsigned int iteration = 0;
        while (iteration < CGAL_MAX_ITERATIONS)
        {                                  
#if OBJ_ON_ARR_DEBUG
          std::cout << "S1_px = " << S1_px << std::endl;
          std::cout << "S1_py = " << S1_py << std::endl;
          std::cout << "S1_pz = " << S1_pz << std::endl;
          std::cout << "S1_dx = " << S1_dx << std::endl;
          std::cout << "S1_dy = " << S1_dy << std::endl;
          std::cout << "S1_dz = " << S1_dz << std::endl<< std::endl;
                  
          std::cout << "S2_px = " << S2_px << std::endl;
          std::cout << "S2_py = " << S2_py << std::endl;
          std::cout << "S2_pz = " << S2_pz << std::endl;
          std::cout << "S2_dx = " << S2_dx << std::endl;
          std::cout << "S2_dy = " << S2_dy << std::endl;
          std::cout << "S2_dz = " << S2_dz << std::endl<< std::endl;
                  
          std::cout << "S3_px = " << S3_px << std::endl;
          std::cout << "S3_py = " << S3_py << std::endl;
          std::cout << "S3_pz = " << S3_pz << std::endl;
          std::cout << "S3_dx = " << S3_dx << std::endl;
          std::cout << "S3_dy = " << S3_dy << std::endl;
          std::cout << "S3_dz = " << S3_dz << std::endl;

#endif

#if OBJ_ON_ARR_DEBUG
          std::cout << "iteration = " << iteration << std::endl;
#endif
          C1 = S3_dz * S1_dx;
          C2 = S3_dx * S1_dz;
          C3 = C1 - C2;
          C4 = S3_px * S3_dz + S3_dx * S1_pz - 
            S3_dx * S3_pz - S3_dz * S1_px;
          C5 = S3_dx * S2_dz - S3_dz * S2_dx;
          C6 = - S3_dx * S1_pz + S3_dx * S2_pz + S3_dz * S1_px - 
            S3_dz * S2_px;
               
          C101 = S3_dz * S2_dy - S2_dz * S3_dy;
          C102 = S3_py * S3_dz + S1_pz * S3_dy - S3_dz * S1_py - 
            S3_pz * S3_dy;
          C103 = S3_dz * S1_dy;
          C104 = S1_dz * S3_dy;
          C106 = S3_dz * S1_py - S1_pz * S3_dy + S2_pz * S3_dy - 
            S3_dz * S2_py;
               
          C201 = C5 * C106 + C6 * C101;
          C202 = C101 * C3 - C5 * C104 + C5 * C103;

#if OBJ_ON_ARR_DEBUG
          std::cout << change_color(CGAL_YELLOW,"  C101 = ",C101) 
                    << std::endl;
#endif
          if (S3_dz != 0 && C101 != 0)
          {
            iteration = CGAL_VALID_ITERATION;
          }
          else
          {
            /* Swap the second and the third equations. */
            if (iteration % 2 == 0)
            {
              std::swap(S3_dy,S3_dz);
              std::swap(S3_py,S3_pz);
              std::swap(S2_dy,S2_dz);
              std::swap(S2_py,S2_pz);
              std::swap(S1_dy,S1_dz);
              std::swap(S1_py,S1_pz);
            }
            else
            {
              /* Swap the second and the first equations. */
              std::swap(S3_dy,S3_dx);
              std::swap(S3_py,S3_px);
              std::swap(S2_dy,S2_dx);
              std::swap(S2_py,S2_px);
              std::swap(S1_dy,S1_dx);
              std::swap(S1_py,S1_px);
            }
            ++iteration;
          }
        }
        CGAL_assertion(iteration != CGAL_MAX_ITERATIONS);
               
#if OBJ_ON_ARR_DEBUG
        std::cout << "S1_px = " << S1_px << std::endl;
        std::cout << "S1_py = " << S1_py << std::endl;
        std::cout << "S1_pz = " << S1_pz << std::endl;
        std::cout << "S1_dx = " << S1_dx << std::endl;
        std::cout << "S1_dy = " << S1_dy << std::endl;
        std::cout << "S1_dz = " << S1_dz << std::endl<< std::endl;

        std::cout << "S2_px = " << S2_px << std::endl;
        std::cout << "S2_py = " << S2_py << std::endl;
        std::cout << "S2_pz = " << S2_pz << std::endl;
        std::cout << "S2_dx = " << S2_dx << std::endl;
        std::cout << "S2_dy = " << S2_dy << std::endl;
        std::cout << "S2_dz = " << S2_dz << std::endl<< std::endl;

        std::cout << "S3_px = " << S3_px << std::endl;
        std::cout << "S3_py = " << S3_py << std::endl;
        std::cout << "S3_pz = " << S3_pz << std::endl;
        std::cout << "S3_dx = " << S3_dx << std::endl;
        std::cout << "S3_dy = " << S3_dy << std::endl;
        std::cout << "S3_dz = " << S3_dz << std::endl;

#endif

        if (C202 == 0)
        {
          Rational C203 = C4 * C101 + C5 * C102;
          Rational t = - C203 / C201;
          /* 
           * (1-t)*(S1.x1 + S1.t * S1.x2) + t * (S2.x1 + S2.t * S2.x2) = 
           *    (S3.x1 + S3.t * S3.x2)
           * S1.t * S1.x2 * (1-t) + S2.t * S2.x2 * t = S3.t * S3.x2 + 
           *    (S3.x1 - (1 - t) S1.x1 - t * S2.x1)
           * C_x = (S3.x1 - (1 - t) * S1.x1 - t * S2.x1)
           *
           * (1-t)*(S1.y1 + S1.t * S1.y2) + t*(S2.y1 + S2.t * S2.y2) =
           *     (S3.y1 + S3.t * S3.y2)
           * S1.t * S1.y2 * (1-t) + S2.t * S2.y2 * t = S3.t * S3.y2 +
           *     (S3.y1 - (1-t) S1.y1 - t * S2.y1)
           * C_y = (S3.y1 - (1 - t) * S1.y1 - t * S2.y1)
           */


#if OBJ_ON_ARR_DEBUG
          std::cout << "t4 = " << t << std::endl;
#endif
          Rational C_x;
          Rational C_y;
          Rational S1_x_temp;
          Rational S2_x_temp;
          Rational S1_y_temp;
          Rational S2_y_temp;
          int status;
          unsigned int iteration = 0;
                    
          do {
            if (iteration == 1)
            {
              /* Swap the second and the third equations. */
              std::swap(S3_dy,S3_dz);
              std::swap(S3_py,S3_pz);
              std::swap(S2_dy,S2_dz);
              std::swap(S2_py,S2_pz);
              std::swap(S1_dy,S1_dz);
              std::swap(S1_py,S1_pz);
            }
            else if (iteration == 2)
            {
              /* Swap the third and the first equations. */
              std::swap(S3_dz,S3_dx);
              std::swap(S3_pz,S3_px);
              std::swap(S2_dz,S2_dx);
              std::swap(S2_pz,S2_px);
              std::swap(S1_dz,S1_dx);
              std::swap(S1_pz,S1_px);
            }
                       
            C_x = (S3_px - (1 - t) * S1_px - t * S2_px);
            C_y = (S3_py - (1 - t) * S1_py - t * S2_py);
            S1_x_temp = (S1_dx * (1-t));
            S2_x_temp = (S2_dx * t);
            S1_y_temp = (S1_dy * (1-t));
            S2_y_temp = (S2_dy * t);
                       
            status = calc_obj_on_arr_t_is_constant(S1_x_temp,
                                                   S2_x_temp,
                                                   S3_dx,
                                                   C_x,
                                                   S1_y_temp,
                                                   S2_y_temp,
                                                   S3_dy,
                                                   C_y);
            ++iteration;
          }
          while (status != LTS_g_func::CGAL_QUERY_SUCCEED && iteration != 3);
                    
          if (status != LTS_g_func::CGAL_QUERY_SUCCEED)
          {
            std::cout << "Unexpected error " << std::endl;
            CGAL_error_msg("Unepxected error t is constant B3 = 0");
            CGAL_error();
          }

          return;
        }

            
        Rational C203 = C4 * C101 + C5 * C102;
        Rational C301 = C102 * C202 + C104 * C203 - C103 * C203;
        Rational C302 = C104 * C201 - C103 * C201 + C202 * C106;
        Rational C303 = C101 * C202;
            
#if OBJ_ON_ARR_DEBUG
        std::cout << "S1->t = (" << C203 << " + t * " << C201 
                  << ") / ((1-t) * " << C202 << ")" << std::endl;
        std::cout << "S2.t = (" << C301 << " + t * " << C302 
                  << ") / (t * " << C303 << ")" << std::endl;
#endif
        m_a0 = C301 * C201 - C203 * C302;
        m_a1 = C202 * C302 + C301 * C202;
        m_b1 = C303 * C202;
        m_b0 = - C303 * C203;

#if OBJ_ON_ARR_DEBUG
        std::cout << this << std::endl;
#endif
        find_bounded_segs(C203,C201,C202,C301,C302,C303,
                          S1_pz,S1_dz,
                          S2_pz,S2_dz,
                          S3_pz,S3_dz);
      }
    }
  }

  void add_bounded_seg(const Rational& bmin, const Rational& bmax, const Rational& A){
    CGAL_precondition(bmin <= bmax);
    if(bmin <= A && A <= bmax){
      // Asymptote splits the interval 
      this->m_bounded_segs.add_bounded_seg(Bounded_seg(Rbound(bmax),Rbound(A),true,false));
      this->m_bounded_segs.add_bounded_seg(Bounded_seg(Rbound(A),Rbound(bmin),false,true));
    }else{
      this->m_bounded_segs.add_bounded_seg(Bounded_seg(Rbound(bmax),Rbound(bmin),true,true));
    }    
  }
  
  ~Lines_through_segments_arr_object()
  {  
  }
  
  /*************************************************************
   *The following function finds the point on S3, that the segment that goes 
   *  through S1 and S2 goes through it.
   *
   * Input:
   *       S1_t       - Represents a point on the segment S1.
   *       S2_t       - Represents a point on the segment S2.
   *
   * Output:
   *       S3_t       - Represents a point on the segment S3.
   * 
   * Function Description:
   * -------------------- 
   * Given 3 segments S1, S2, S3.
   * S1:(S1.x1 + S1.t * S1.x2, S1.y1 + S1.t * S1.y2, S1.z1 + S1.t * S1.z2)
   * S2:(S2.x1 + S2.t * S2.x2, S2.y1 + S2.t * S2.y2, S2.z1 + S2.t * S2.z2)
   * S3:(S3.x1 + S3.t * S3.x2, S3.y1 + S3.t * S3.y2, S3.z1 + S3.t * S3.z2)
   *
   *  A segment equation of general segment which intersects with both S1 and
   *  S2
   *    will look like: (1-t)*I1 + t*I2,
   * Where I1, I2 are the intersection points of this segment with S1 and S2
   *   respectivly.
   * If this segment also intersect S3, the point I3 = (1-t)*I1 + t*I2.
   * In order to get the hyperbola the following 3 equation are solved:
   *    (1-t)*(S1.x1 + S1.t * S1.x2) + t*(S2.x1 + S2.t * S2.x2) = 
   *       (S3.x1 + S3.t * S3.x2)
   *    (1-t)*(S1.y1 + S1.t * S1.y2) + t*(S2.y1 + S2.t * S2.y2) = 
   *       (S3.y1 + S3.t * S3.y2)
   *    (1-t)*(S1.z1 + S1.t * S1.z2) + t*(S2.z1 + S2.t * S2.z2) = 
   *       (S3.z1 + S3.t * S3.z2)
   *
   * From previous calculations we get that:
   *
   * t = ((C202 * S1.t - C203) / (C201 + S1.t * C202))
   *
   *  
   *
   * From the third eqution we will get S3.t in terms of t.
   *
   *    If (S3.z2 == 0) We swap the equations s.t S3.z2 won't be 0.
   *         (First we try to swap with the second equation and then with the
   *          first.)
   *
   *        If S1.z2 == 0 && S2.z2 == 0 && S3.z2 == 0 the segments are at the
   *           same plane or on parallel planes.
   *
   *  S3.t = (((1-t)*(S1.z1 + S1.t * S1.z2) + t*(S2.z1 + S2.t * S2.z2) - S3.z
   *  1)/S3.z2)
   *
   *  S3.t = (((1-t)*(S1.z1 + S1.t * S1.z2) + t*(S2.z1 + S2.t * S2.z2) - 
   *      S3.z1)/ S3.z2)
   *************************************************************/
  Algebraic getS3_t(Algebraic &S1_t,
                    Algebraic &S2_t)
  {
    /* In order to get presentation of a segment L in a form of (x1,y1,z1)
     *  + t * (x2,y2,z2),
     * We will get two arbitrary points on L and set (x1,y1,z1) coordinates 
     *   to be the first point and
     * (x2,y2,z2) to be (x1,y1,z1) minus the coordinates of the second 
     *  point. 
     */
   
    Rational_point_3 S1_P1 = m_S1->source();
    Rational_point_3 S1_P2 = m_S1->target();

    Rational_point_3 S2_P1 = m_S2->source();
    Rational_point_3 S2_P2 = m_S2->target();

    Rational_point_3 S3_P1 = m_creator_segment->source();
    Rational_point_3 S3_P2 = m_creator_segment->target();
   
    Rational S1_px = S1_P1.x();
    Rational S1_py = S1_P1.y();
    Rational S1_pz = S1_P1.z();
    Rational S1_dx = S1_P2.x() - S1_px;
    Rational S1_dy = S1_P2.y() - S1_py;
    Rational S1_dz = S1_P2.z() - S1_pz;
   
    Rational S2_px = S2_P1.x();
    Rational S2_py = S2_P1.y();
    Rational S2_pz = S2_P1.z();
    Rational S2_dx = S2_P2.x() - S2_px;
    Rational S2_dy = S2_P2.y() - S2_py;
    Rational S2_dz = S2_P2.z() - S2_pz;

    Rational S3_px = S3_P1.x();
    Rational S3_py = S3_P1.y();
    Rational S3_pz = S3_P1.z();
    Rational S3_dx = S3_P2.x() - S3_px;
    Rational S3_dy = S3_P2.y() - S3_py;
    Rational S3_dz = S3_P2.z() - S3_pz;

    if (S3_dz == 0)
    {
      /* Swap the second and the third equations. */
      if (S3_dy != 0)
      {
        std::swap(S3_dy,S3_dz);
        std::swap(S3_py,S3_pz);
        std::swap(S2_dy,S2_dz);
        std::swap(S2_py,S2_pz);
        std::swap(S1_dy,S1_dz);
        std::swap(S1_py,S1_pz);
      }
      /* Swap the first and the third equations. */
      else if (S3_dx != 0)
      {
        std::swap(S3_dx,S3_dz);
        std::swap(S3_px,S3_pz);
        std::swap(S2_dx,S2_dz);
        std::swap(S2_px,S2_pz);
        std::swap(S1_dx,S1_dz);
        std::swap(S1_px,S1_pz);
      }
      else
      {
        /* S3_dz == 0 && S3_dy == 0 && S3_dx == 0 -> The segment is
           a point. */
        CGAL_error();
      }
    }
         
    Rational C1 = S3_dz * S1_dx;
    Rational C2 = S3_dx * S1_dz;
    Rational C3 = C1 - C2;
    Rational C4 = S3_px * S3_dz + S3_dx * S1_pz - S3_dx * S3_pz -
      S3_dz * S1_px;
    Rational C5 = S3_dx * S2_dz - S3_dz * S2_dx;
    Rational C6 = - S3_dx * S1_pz + S3_dx * S2_pz + S3_dz * S1_px -
      S3_dz * S2_px;

    Rational C101 = S3_dz * S2_dy - S2_dz * S3_dy;
    Rational C102 = S3_py * S3_dz + S1_pz * S3_dy - S3_dz * S1_py -
      S3_pz * S3_dy;
    Rational C103 = S3_dz * S1_dy;
    Rational C104 = S1_dz * S3_dy;
    Rational C106 = S3_dz * S1_py - S1_pz * S3_dy + S2_pz * S3_dy -
      S3_dz * S2_py;

    Rational C202 = C101 * C3 - C5 * C104 + C5 * C103;
    Rational C201 = C5 * C106 + C6 * C101;
    Rational C203 = C4 * C101 + C5 * C102;
         
    std::cout << "C201 = " << C201 << std::endl;
    std::cout << "C202 = " << C202 << std::endl;
    std::cout << "C203 = " << C203 << std::endl;
         
    if ((C201 + S1_t * C202) == 0)
    {
      /* t equals infinity => S3_t = -+ infinity - no solutions */
      return Algebraic(-1);
    }
         
    Algebraic t = ((C202 * S1_t - C203) / (C201 + S1_t * C202));
        
    Algebraic S3_t = (((1-t)*(S1_pz + S1_t * S1_dz) + 
                       t*(S2_pz + S2_t * S2_dz) - S3_pz)/S3_dz);
    return S3_t;
  }
        
  Rational_segment_3* get_creator_segment() {return &m_creator_segment;}
      
  std::string to_string()
  {
    std::ostringstream o;
    if (m_obj_on_arr_type == CGAL_DEGENERATE_HYPERBOLA)
    {
      o << "Degenerate hyperbola x = " << m_vertical_segment_x_coordinate
        << " y = " << m_horizontal_segment_y_coordinate << std::endl;

      o << "Vertical segments bounded segments " <<
        m_vertical_segment_bounded_segs << std::endl;

      o << "Horizontal segments bounded segments " <<
        m_horizontal_segment_bounded_segs << std::endl;
            
    }
    else
    {
      o << "y = (" << m_a0 <<" + x * " << m_a1 << ") / (" << m_b0
        << " + x * " << m_b1 <<")" << std::endl;
      o << m_bounded_segs << std::endl;
    }
        
    return o.str();
  }

private:
  /*************************************************************
   * Function description:
   * -------------------- 
   * The following functions gets an inequality and returns the solutions
   *  to this inequality.
   *
   * Input:
   *     Equation of the form (a1 * x + a2)/(b1 * x + b2) >= c.
   *
   * Output:
   *   Solution of the type min_x <= x <= max_x or (x <= min_x || max_x <= x)
   *
   *      The variables min_x_exists and max_x_exists define the type of 
   *        the bound.
   *     If the the value of min_x_exists == true and x is bigger then min_x,
   *         then x >= min_x otherwise x > min_x.  
   *
   *       b1 != 0
   *       --------
   *
   *           c * b1 == a1
   *          -------------
   *             Case 1: If ((b2 + b1 * x) > 0 && b1 > 0) ||
   *                        ((b2 + b1 * x) < 0 && b1 < 0)
   *      
   *                     If (a2 >= c * b2)
   *                         x >= -(b2/b1)
   *                   
   *                    else
   *                          No solutions
   *
   *             Case 2: If ((b2 + b1 * x) < 0 && b1 < 0) ||
   *                        ((b2 + b1 * x) > 0 && b1 > 0)
   *      
   *                     If (a2 >= c * b2)
   *                          No solutions
   *                   
   *                    else
   *                         x <= -(b2/b1)
   *
   *
   *           c * b1 != a1
   *          -------------
   *
   *
   *         Case 1: If ((b2 + b1 * x) > 0 && b1 > 0) ||
   *                    ((b2 + b1 * x) < 0 && b1 < 0)
   *
   *                   (a1 * x + a2) >=  (b2 + b1 * x) * c
   *                   
   *                   a1 * x + a2 >= b2 * c + b1 * x * c
   *  
   *                   a1 * x - b1 * x * c >= b2 * c - a2
   *  
   *                   x * (a1 - b1 * c) >= b2 * c - a2
   *  
   *                   Let CMLS20 = (a1 - b1 * c)
   *                   Let CMLS21 = b2 * c - a2
   *
   *                   x * CMLS20 >= CMLS21
   *
   *                   Case 1.1: CMLS20 > 0 
   *
   *                         x >= (CMLS21/CMLS20) 
   *
   *                         Case 1.1.1: b1 > 0 
   *
   *                               x >= max((CMLS21/CMLS20),(- b2/b1))
   *
   *                               if ((-b2/b1) > (CMLS21/CMLS20))
   * 
   *                                  x > (-b2/b1) 
   *
   *                         Case 1.1.2: b1 < 0 
   *
   *                               - b2/b1 > x >= (CMLS21/CMLS20) 
   *
   *                   Case 1.2: CMLS20 < 0 
   * 
   *                         x <= (CMLS21/CMLS20)
   *
   *                         Case 1.2.1: b1 > 0 
   *
   *                               - b2/b1 < x <= (CMLS21/CMLS20) 
   *
   *                         Case: 1.2.2: b1 < 0 
   *
   *                               x <= min((-b2/b1),(CMLS21/CMLS20))
   *
   *                               if ((-b2/b1) < (CMLS21/CMLS20))
   * 
   *                                  x < (-b2/b1) 
   *
   *         Case 2: If ((b2 + b1 * x) < 0 && b1 > 0) ||
   *                      ((b2 + b1 * x) > 0 && b1 < 0)
   *
   *                   x * CMLS20 <= CMLS21
   *
   *                   Case 2.1: CMLS20 > 0 
   *
   *                         x <= (CMLS21/CMLS20) 
   *
   *                         Case 2.1.1: b1 > 0 
   *
   *                               x <= min((CMLS21/CMLS20),(- b2/b1))
   *
   *                               if ((-b2/b1) < (CMLS21/CMLS20))
   * 
   *                                  x < (-b2/b1) 
   *
   *
   *                         Case 2.1.2: b1 < 0 
   *
   *                               - b2/b1 < x <= (CMLS21/CMLS20) 
   *
   *                   Case 2.2: CMLS20 < 0 
   * 
   *                         x >= (CMLS21/CMLS20)
   *
   *                         Case 2.2.1: b1 > 0 
   *
   *                               - b2/b1 > x >= (CMLS21/CMLS20) 
   *
   *                         Case 2.2.2: b1 < 0 
   *
   *                               x >= max((-b2/b1),(CMLS21/CMLS20))
   *
   *                               if ((-b2/b1) > (CMLS21/CMLS20))
   * 
   *                                  x > (-b2/b1) 
   *
   *     b1 == 0
   *    ---------
   *    (a1 * x + a2)/(b2) >= c.
   *
   *     a1*x/b2 >= c - a2/b2
   *
   *     Case 1: b2/a1 >= 0 
   *
   *                x >= (c * b2 -a2)/a1.
   *
   *     Case 2: b2/a1 < 0 
   *
   *                x <= (c * b2 -a2)/a1.
   *
   **************************************************************/
  void get_bounded_segments_from_inequality(const Rational& a1,
                                            const Rational& a2,
                                            const Rational& b1,
                                            const Rational& b2,
                                            const Rational& c,
                                            Rbound& min_x,
                                            bool& min_x_exists,
                                            Rbound& max_x,
                                            bool& max_x_exists,
                                            bool& splitted_solution)
  {
    Rational CMLS20 = (a1 - b1 * c);
    Rational CMLS21 = b2 * c - a2;
#if OBJ_ON_ARR_DEBUG    
    /* (a1 * x + a2)/(b1 * x + b2) >= c. */
    std::cout << "a1 = " << a1 << std::endl;
    std::cout << "a2 = " << a2 << std::endl;

    std::cout << "b1 = " << b1 << std::endl;
    std::cout << "b2 = " << b2 << std::endl;
    std::cout << "c = " << c << std::endl;
#endif

    if (b1 == 0)
    {
      /* Not hyperbola error */
      CGAL_assertion(a1 != 0);

      /* Case 1: b2/a1 >= 0  => x >= (c * b2 -a2)/a1. */
      if (b2/a1 >= 0)
      {
                
        if (m_bound_s1_s2)
        {
          max_x = Rational(1);
          min_x = max(Rational(0),(c * b2 -a2)/a1);
        }
        else
        {
          max_x = LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY;
          min_x = (c * b2 -a2)/a1;
        }
      }
      /* Case 2: b2/a1 <= 0  => x <= (c * b2 -a2)/a1. */
      else
      {
        if (m_bound_s1_s2)
        {
          min_x = Rational(0);
          max_x = min(Rational(1),(c * b2 -a2)/a1);
        }
        else
        {
          min_x = LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY;
          max_x = (c * b2 -a2)/a1;
        }
      }
    }
    else if (CMLS20 == 0)
    {
      if (b1 > 0)
      {
        if (a2 >= c * b2)
        {
          /* x > -(b2/b1) */
          if (m_bound_s1_s2)
          {
            max_x = Rational(1);
            min_x = max(Rational(0),-(b2/b1));
          }
          else
          {
            max_x = LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY;
            min_x = -(b2/b1);
          }
        }
        else
        {
          /* x < -(b2/b1) */
          if (m_bound_s1_s2)
          {
            min_x = Rational(0);
            max_x = min(Rational(1),-(b2/b1));
          }
          else
          {
            min_x = LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY;
            max_x = -(b2/b1);
          }
        }
      }
      else /* (b1 < 0) */
      {
        if (a2 >= c * b2)
        {
          /* x < -(b2/b1) */
          if (m_bound_s1_s2)
          {
            min_x = Rational(0);
            max_x = min(Rational(1),-(b2/b1));
          }
          else
          {
            min_x = LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY;
            max_x = -(b2/b1);
          }
        }
        else
        {
          /* x >= -(b2/b1) */
          if (m_bound_s1_s2)
          {
            max_x = Rational(1);
            min_x = max(Rational(0),-(b2/b1));
          }
          else
          {
            max_x = LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY;
            min_x = -(b2/b1);
          }
        }
      }
    }
    else if (((CMLS20 > 0) && (b1 > 0)) ||
             ((CMLS20 < 0) && (b1 < 0)))
    {
      /* There is a gap in this case which creates splitted bounded 
         segment. One between 0 and max_x, and the other between max_x and
         1.
      */
      splitted_solution = true;
                  
      /* Case 1.1.1 & Case 2.2.2 :x >= max((-b2/b1),(CMLS21/CMLS20)) */
      if (m_bound_s1_s2)
      {
        min_x = max(Rational(0),(CMLS21/CMLS20));
        min_x = max(min_x.bound(),(- b2/b1));
      }
      else
      {
        min_x = max((- b2/b1),(CMLS21/CMLS20));
      }
            
            
#if OBJ_ON_ARR_DEBUG         
      std::cout << " Case 1.1.1 & Case 2.2.2 (CMLS21/CMLS20) = "
                << (CMLS21/CMLS20) 
                << "  (- b2/b1) = "
                <<(- b2/b1)<<std::endl;
#endif
      /* Case 1.2.2 & Case 2.1.1 : x <= min((CMLS21/CMLS20),(- b2/b1)) */
      if (m_bound_s1_s2)
      {
        max_x = min(Rational(1),(-b2/b1));
        max_x = min(max_x.bound(),(CMLS21/CMLS20));
      }
      else
      {
        max_x = min((CMLS21/CMLS20),(-b2/b1));
      }
            
            
#if OBJ_ON_ARR_DEBUG         
      std::cout << " Case 1.2.2 & Case 2.1.1: (CMLS21/CMLS20) = "
                << (CMLS21/CMLS20) 
                << "  (- b2/b1) = "
                << (- b2/b1)<<std::endl;
      std::cout << "max_x = " << max_x << std::endl;
      std::cout << "min_x = " << min_x << std::endl;
#endif
                  
    }
    else
    {
      splitted_solution = false;
      if ((CMLS21/CMLS20) > (- b2/b1))
      {
        /* Case 2.1.2 & Case 1.2.1 :  - b2/b1 <= x <= (CMLS21/CMLS20) */
        if (m_bound_s1_s2)
        {
          max_x = min(Rational(1),(CMLS21/CMLS20));
          min_x = max(Rational(0),-b2/b1);
        }
        else
        {
          max_x = CMLS21/CMLS20;
          min_x = -b2/b1;
        }
               
#if OBJ_ON_ARR_DEBUG         
        std::cout << " Case 2.1.2 & Case 1.2.1: (CMLS21/CMLS20) = " 
                  << (CMLS21/CMLS20) << "  (- b2/b1) = "
                  << (- b2/b1) << std::endl;
        std::cout << "max_x = " << max_x << std::endl;
        std::cout << "min_x = " << min_x << std::endl;
#endif
      }
      else
      {
        /* Case 1.1.2 & Case 2.2.1 : - b2/b1 >= x >= (CMLS21/CMLS20)*/
        if (m_bound_s1_s2)
        {
          min_x = max(Rational(0),(CMLS21/CMLS20));
          max_x = min(Rational(1),(-b2/b1));
        }
        else
        {
          min_x = CMLS21/CMLS20;
          max_x = -b2/b1;
        }
               
               
#if OBJ_ON_ARR_DEBUG         
        std::cout << " Case 1.1.2 & Case 2.2.1: (CMLS21/CMLS20) = "
                  << (CMLS21/CMLS20) <<"  (- b2/b1) = " 
                  <<(- b2/b1)<<std::endl;
        std::cout << "max_x = " << max_x << std::endl;
        std::cout << "min_x = " << min_x << std::endl;
#endif
      }
    }

    if (b1 == 0)
    {
      min_x_exists = true;
      max_x_exists = true;
    }
    else
    {
      if (                                
          min_x == (- b2/b1) &&
          !(((b2 == 0 && a2 == 0)) || ((b2 != 0) && (a1/b1) == (a2/b2))))
      {
        min_x_exists = false;
      }
      else
      {
        min_x_exists = true;
      }
            
      if (max_x == (- b2/b1) &&
          !(((b2 == 0 && a2 == 0)) || ((b2 != 0) && (a1/b1) == (a2/b2))))
      {
        max_x_exists = false;
      }
      else
      {
        max_x_exists = true;
      }
    }
         
    /* In case the solution is splitted and max_x = min_x the bound is
       [0,1] */
    if (splitted_solution && min_x == max_x &&
        (max_x_exists || min_x_exists))
    {
      splitted_solution = false;
      if (m_bound_s1_s2)
      {
        min_x = 0;
        max_x = 1;
        max_x_exists = true;
        min_x_exists = true;

      }
      else
      {
        min_x = LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY;
        max_x = LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY;
        max_x_exists = false;
        min_x_exists = false;

      }
    }
  }
      
  /*************************************************************
   * Function description:
   * -------------------- 
   * Given 3 blue segments S1, S2, S3.
   * S1:(S1.x1 + S1.t * S1.x2, S1.y1 + S1.t * S1.y2, S1.z1 + S1.t * S1.z2)
   * S2:(S2.x1 + S2.t * S2.x2, S2.y1 + S2.t * S2.y2, S2.z1 + S2.t * S2.z2)
   * S3:(S3.x1 + S3.t * S3.x2, S3.y1 + S3.t * S3.y2, S3.z1 + S3.t * S3.z2)
   *
   * The function finds the bounded segment of the hyperbola according to
   *  the segments.
   * The bounded segment will always be at the positive unit square since
   *   0<=S1_t<=1 and 0<=S2_t<=1
   *
   * A sgement equation of general segment which intersects with both S1 and
   *   S2 will look like: (1-t)*I1 + t*I2,
   * Where I1, I2 are the intersection points of this segment with S1 and
   *   S2 respectivly.
   * If this segment also intersect S3, the point I3 = (1-t)*I1 + t*I2.
   * In order to get the hyperbola the following 3 equation are solved:
   *    (1-t)*(S1.x1 + S1.t * S1.x2) + t*(S2.x1 + S2.t * S2.x2) =
   *      (S3.x1 + S3.t * S3.x2)
   *    (1-t)*(S1.y1 + S1.t * S1.y2) + t*(S2.y1 + S2.t * S2.y2) =
   *      (S3.y1 + S3.t * S3.y2)
   *    (1-t)*(S1.z1 + S1.t * S1.z2) + t*(S2.z1 + S2.t * S2.z2) =
   *      (S3.z1 + S3.t * S3.z2)
   *
   * From the third eqution we will get S3.t in terms of t.
   *
   *  S3.t = 
   *   (((1-t)*(S1.z1 + S1.t * S1.z2) + t*(S2.z1 + S2.t * S2.z2) - S3.z1)
   *    /S3.z2)
   *
   * Given that S3 is a segment, we get 0 <= S3.t <= 1.
   * We will bound the values 
   *
   * From the above calculations we get.
   * S1.t = (C203 +  t * C201) / ((1 - t) * C202)
   * S2.t = (C301 +  t * C302) / (t * C303)
   *
   * Input:
   *      C203,C201,C202,C301,C302,C303 - coefficients for the 
   *        representations of S1.t and S2.t
   *      S1.z1, S1.z2, S2.z1, S2.z2, S3.z1, S3.z2 - coefficients for
   *         the representation of S3.t
   *
   * Output:
   *     max_x,max_y,min_x,min_y      
   *
   * Function detailed description:
   * -------------------- 
   * S3.t =
   *    (((1-t)*(S1.z1 + ((C203 +  t * C201) / ((1 - t) * C202)) * S1.z2) +
   *    t * (S2.z1 + ((C301  +  t * C302) / (t * C303)) * S2.z2) - S3.z1)/
   *    S3.z2)
   * 
   * Multiply by S3.z2.
   *
   * S3.t * S3.z2 =
   *    (((1-t)*(S1.z1 + ((C203 +  t * C201) / ((1 - t) * C202)) * S1.z2) +
   *    t * (S2.z1 + ((C301  +  t * C302) / (t * C303)) * S2.z2) - S3.z1))
   *
   * Open the first parentheses multiply by 1-t:
   *
   * S3.t * S3.z2 =
   *    ((((1-t) * S1.z1 + S1.z2 * ((C203 +  t * C201) / (C202)) ) +
   *    t * (S2.z1 + ((C301  +  t * C302) / (t * C303)) * S2.z2) - S3.z1))
   *  
   * Open the second parentheses multiply by t:
   *
   * S3.t * S3.z2 =
   *    ((1-t) * S1.z1 + (S1.z2 * C203 +  t * S1.z2 * C201) / C202 ) +
   *    (t * S2.z1 + (S2.z2 * C301  +  t * S2.z2 * C302) / C303 - S3.z1)
   *
   *Devide the variables to two groups, ones that depends on t, and the other
   *    of free scalars.
   * Multiply by C303 and C202
   *
   * S3.t * S3.z2 * C202 * C303 =
   *    ((1 - t) * S1.z1 * C202 * C303 + 
   *    (S1.z2 * C203 * C303 +  t * S1.z2 * C201 * C303) ) +
   *    (t * S2.z1 * C202 * C303 + (S2.z2 * C301 * C202 +  
   *    t * S2.z2 * C302 * C202) - S3.z1 * C202 * C303)
   *
   * S3.t * S3.z2 * C202 * C303 =
   *    (S1.z1 * C202 * C303  - t * S1.z1 * C202 * C303 + 
   *    S1.z2 * C203 * C303 +  t * S1.z2 * C201 * C303 +
   *    t * S2.z1 * C202 * C303 + S2.z2 * C301 * C202 + 
   *    t * S2.z2 * C302 * C202 - S3.z1 * C202 * C303)
   *
   * S3.t * S3.z2 * C202 * C303 =
   *    (t * (S1.z2 * C201 * C303 + S2.z1 * C202 * C303 - S1.z1 * C202 * C303
   *    + S2.z2 * C302 * C202) + S1.z1 * C202 * C303 + 
   *    S1.z2 * C203 * C303 + S2.z2 * C301 * C202 - S3.z1 * C202 * C303)
   *
   * Let CMLS1 = (S1.z2 * C201 * C303 + S2.z1 * C202 * C303 - 
   *             S1.z1 * C202 * C303 + S2.z2 * C302 * C202)
   * Let CMLS2 = S1.z1 * C202 * C303 + S1.z2 * C203 * C303 + 
   *             S2.z2 * C301 * C202 - S3.z1 * C202 * C303
   * Let CMLS3 = S3.z2 * C202 * C303
   *
   * S3.t * CMLS3 = (t * CMLS1 + CMLS2)
   *
   * S3.t = (t * CMLS1 + CMLS2)/CMLS3
   *
   *  Limit S1.t by S3.t >= 0 
   * --------------------------
   *
   * => (t * CMLS1 + CMLS2)/CMLS3 >= 0
   *
   * (t * CMLS1/CMLS3) >= -(CMLS2/CMLS3)
   * 
   * if (CMLS3/CMLS1 > 0) (Multiply by (CMLS3/CMLS1))
   *     t >= -(CMLS2/CMLS1)
   * else
   *     t <= -(CMLS2/CMLS1)
   *
   * Isolate t at the following equation:
   *
   *  S1.t = (C203 +  t * C201) / ((1 - t) * C202)
   *
   *  ((1 - t) * C202) * S1.t =  (C203 +  t * C201)
   *
   *  C202 * S1.t - t * S1.t * C202 = C203 +  t * C201
   *
   *  C202 * S1.t - C203 = t * C201 + t * S1.t * C202
   * 
   *  C202 * S1.t - C203 = t * (C201 + S1.t * C202)
   *
   *  t = ((C202 * S1.t - C203) / (C201 + S1.t * C202))
   *
   * Case 1: if (CMLS3/CMLS1 > 0) (Multiply by (CMLS3/CMLS1))
   *         t >= -(CMLS2/CMLS1)
   *
   *         We get that ((C202 * S1.t - C203) / (C201 + S1.t * C202)) >=
   *          -(CMLS2/CMLS1)
   *
   * Case 2: if (CMLS3/CMLS1 < 0) (Multiply by (CMLS3/CMLS1))
   *         t <= -(CMLS2/CMLS1)
   *
   *         We get that ((C202 * S1.t - C203) / (C201 + S1.t * C202)) <= 
   *           -(CMLS2/CMLS1)
   * 
   *  Limit S1.t by S3.t <= 1 
   * --------------------------
   *
   * S3.t = (t * CMLS1 + CMLS2)/CMLS3
   *
   * S3.t <= 1 
   * => (t * CMLS1 + CMLS2)/CMLS3 <= 1
   *
   * (t * CMLS1/CMLS3) <= 1 - (CMLS2/CMLS3)
   * 
   * if (CMLS3/CMLS1 > 0) (Multiply by (CMLS3/CMLS1))
   *     t <= CMLS3/CMLS1 -(CMLS2/CMLS1)
   * else
   *     t >= CMLS3/CMLS1 -(CMLS2/CMLS1)
   *
   * Isolate t at the following equation:
   *
   *  S1.t = (C203 +  t * C201) / ((1 - t) * C202)
   *
   *  ((1 - t) * C202) * S1.t =  (C203 +  t * C201)
   *
   *  C202 * S1.t - t * S1.t * C202 = C203 +  t * C201
   *
   *  C202 * S1.t - C203 = t * C201 + t * S1.t * C202
   * 
   *  C202 * S1.t - C203 = t * (C201 + S1.t * C202)
   *
   *  t = ((C202 * S1.t - C203) / (C201 + S1.t * C202))
   *
   * Case 1: if (CMLS3/CMLS1 > 0) (Multiply by (CMLS3/CMLS1))
   *         t <= CMLS3/CMLS1 -(CMLS2/CMLS1)
   *
   *         Solve the following equation:
   *                   ((C202 * S1.t - C203) / (C201 + S1.t * C202)) <= 
   *                     CMLS3/CMLS1 -(CMLS2/CMLS1)
   *
   * Case 2: if (CMLS3/CMLS1 < 0) (Multiply by (CMLS3/CMLS1))
   *         t >= CMLS3/CMLS1 -(CMLS2/CMLS1)
   *
   *         We get that ((C202 * S1.t - C203) / (C201 + S1.t * C202)) >=
   *              CMLS3/CMLS1 -(CMLS2/CMLS1)
   *
   *********************************************************************/
  void find_bounded_segs(Rational& C203,
                         Rational& C201,
                         Rational& C202,
                         Rational& C301,
                         Rational& C302,
                         Rational& C303,
                         const Rational& S1_pz,
                         const Rational& S1_dz,
                         const Rational& S2_pz,
                         const Rational& S2_dz,
                         const Rational& S3_pz,
                         const Rational& S3_dz)
  {
    Rational CMLS1 = S1_dz * C201 * C303 + S2_pz * C202 * C303 - 
      S1_pz * C202 * C303 + S2_dz * C302 * C202;
    Rational CMLS2 = S1_pz * C202 * C303 + S1_dz * C203 * C303 + 
      S2_dz * C301 * C202 - S3_pz * C202 * C303;
    Rational CMLS3 = S3_dz * C202 * C303;

#if OBJ_ON_ARR_DEBUG
    std::cout << " S3.t = (t * " << CMLS1 << " + " << CMLS2 << ")/" 
              << CMLS3 << std::endl;
#endif
         
    {
      if (CMLS1 == 0)
      {
        /* Devision by 0 S3_t equals infinity */
        CGAL_assertion(CMLS3 != 0);

        /* S3.t * CMLS3 = (t * CMLS1 + CMLS2) */
        /* S3.t * CMLS3 = CMLS2 */
        /* Limit S1.t by 1 >= S3.t >= 0 */
        if (CMLS2/CMLS3 >= 0 && CMLS2/CMLS3 <= 1)
          // (CMLS3 > 0 && CMLS2/CMLS3 >= 0 && CMLS2/CMLS3 <= 1) ||
          // (CMLS3 < 0 && CMLS2/CMLS3 <= 0 && CMLS2/CMLS3 <= 1))
        {
          /* Use the all segment [0,1] as the bounded segment */
          if (m_bound_s1_s2)
          {
            this->m_bounded_segs.add_bounded_seg
              (Bounded_seg(Rbound(Rational(1)), Rbound(Rational(0)), true,
                           true));
          }
          else
          {
            this->m_bounded_segs.add_bounded_seg
              (Bounded_seg(LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY,
                           LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY, 
                           false, false));
          }
        }
        else/* No solutions at this case - no bounded segments.*/
        {
#if OBJ_ON_ARR_DEBUG
          std::cout << "CMLS2 = " << CMLS2 << std::endl;
          std::cout << "CMLS3 = " << CMLS3 << std::endl;
          std::cout << "No solutions " << std::endl;
#endif
        }
      }
      else
      {
        Bounded_segs_vector first_bounded_segs;
        Bounded_segs_vector second_bounded_segs;
        bool first_splitted = false;
        bool second_splitted = false;

        Rbound min_x;
        Rbound max_x;
        bool min_x_exists;
        bool max_x_exists;
            
        /* S3.t >= 0 */ 
        /* Case 1: if (CMLS3/CMLS1 > 0)*/
        if (CMLS3/CMLS1 > 0)
        {
#if OBJ_ON_ARR_DEBUG         
          std::cout << "Case 1: if (CMLS3/CMLS1 > 0) " << std::endl;
#endif
          get_bounded_segments_from_inequality(C202,
                                               -C203,
                                               C202,
                                               C201,
                                               (-CMLS2/CMLS1),
                                               min_x,
                                               min_x_exists,
                                               max_x,
                                               max_x_exists,
                                               first_splitted);
        }
        /* Case 2: if (CMLS3/CMLS1 < 0) */
        else
        {
#if OBJ_ON_ARR_DEBUG         
          std::cout << "Case 2: if (CMLS3/CMLS1 < 0) " << std::endl;
#endif
          get_bounded_segments_from_inequality(-C202,
                                               C203,
                                               C202,
                                               C201,
                                               (CMLS2/CMLS1),
                                               min_x,
                                               min_x_exists,
                                               max_x,
                                               max_x_exists,
                                               first_splitted);
        }
         
         
        if (first_splitted)
        {
          if (m_bound_s1_s2)
          {
            first_bounded_segs.add_bounded_seg
              (Bounded_seg(max_x, Rbound(Rational(0)), max_x_exists, true));
            first_bounded_segs.add_bounded_seg
              (Bounded_seg(Rbound(Rational(1)), min_x, true, min_x_exists));
          }
          else
          {
            first_bounded_segs.add_bounded_seg
              (Bounded_seg(max_x, LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY,
                           max_x_exists, false));
            first_bounded_segs.add_bounded_seg
              (Bounded_seg(LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY,min_x,
                           false,
                           min_x_exists));
          }
        }
        else
        {
          if (m_bound_s1_s2)
          {
            first_bounded_segs.add_bounded_seg(Bounded_seg(max_x, min_x,
                                                           max_x_exists,
                                                           min_x_exists));
          }
          else
          {
            first_bounded_segs.add_bounded_seg(Bounded_seg(max_x, min_x,
                                                           max_x_exists,
                                                           min_x_exists));
          }
        }

        /* S3.t <= 1 */ 
        /* Case 1: if (CMLS3/CMLS1 > 0) */
        if (CMLS3/CMLS1 > 0)
        {
            
#if OBJ_ON_ARR_DEBUG         
          std::cout << "Case 1: if (CMLS3/CMLS1 > 0) " << std::endl;
#endif
          get_bounded_segments_from_inequality(-C202,
                                               C203,
                                               C202,
                                               C201,
                                               - (CMLS3/CMLS1) + 
                                               (CMLS2/CMLS1),
                                               min_x,
                                               min_x_exists,
                                               max_x,
                                               max_x_exists,
                                               second_splitted);
            
        }
        /* Case 2: if (CMLS3/CMLS1 < 0) */
        else
        {
#if OBJ_ON_ARR_DEBUG         
          std::cout << "Case 2: if (CMLS3/CMLS1 < 0)" << std::endl;
#endif
          get_bounded_segments_from_inequality(C202,
                                               -C203,
                                               C202,
                                               C201,
                                               (CMLS3/CMLS1) - 
                                               (CMLS2/CMLS1),
                                               min_x,
                                               min_x_exists,
                                               max_x,
                                               max_x_exists,
                                               second_splitted);
        }

        if (second_splitted)
        {
          if (m_bound_s1_s2)
          {
            second_bounded_segs.add_bounded_seg(Bounded_seg(max_x,
                                                            Rbound(Rational(0)),
                                                            max_x_exists,
                                                            true));
            second_bounded_segs.add_bounded_seg(Bounded_seg(Rbound(Rational(1)),
                                                            min_x, true,
                                                            min_x_exists));
          }
          else
          {
            second_bounded_segs.add_bounded_seg
              (Bounded_seg(max_x, LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY,
                           max_x_exists, false));
            second_bounded_segs.add_bounded_seg
              (Bounded_seg(LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY,min_x, false,
                           min_x_exists));
          }
        }
        else
        {
          second_bounded_segs.add_bounded_seg(Bounded_seg(max_x, min_x,
                                                          max_x_exists,
                                                          min_x_exists));
        }
        this->m_bounded_segs.merge_bounded_segs_vectors(first_bounded_segs,
                                                        second_bounded_segs);
      }
    }
  }
      

  /*************************************************************
   * Function description:
   * -------------------- 
   * The following functions gets a bounded segment and returns all the
   *   arcs of the hyperbola
   * at this segment. The segment is bounded at y = 0, y = 1, x = 0, x = 1.
   *
   * When y = 1 -> b1 * x + b0 = a0 + x * a1
   * b1 * x - x * a1 = a0 - b0
   * x = (a0 - b0) / (b1 - a1)
   *
   * When y = 0 -> a0 + a1 * x = 0 -> x = -(a0/a1)
   *
   *************************************************************/

  void get_bounded_segments_at_positive_unit_square
  (Rational& x_intersection_y_equals_0,
   Rational& x_intersection_y_equals_1,
   Bounded_seg res_bounded_segs[2],
   bool& res_splitted)
  {
    res_splitted = false;

    /* Case 1: The hyperbola intersects with the segments y = 1 and y = 0
       at the positive unit square. */
    if ((0 <= x_intersection_y_equals_1 &&
         x_intersection_y_equals_1 <= 1) &&
        (0 <= x_intersection_y_equals_0 &&
         x_intersection_y_equals_0 <= 1))
    {
      /* Case 1.1: The hyperbola intersect first with the segment y = 0.*/
      if (x_intersection_y_equals_0 < x_intersection_y_equals_1)
      {
                
        /* Check if the hyperbola is "splitted" at the unit square,
           if so we will take two arcs, one from
           (min_x,x_intersection_y_equals_0] and the 
           second [x_intersection_y_equals_1,max_x).
           Otherwise take one arc [x_intersection_y_equals_0,
           x_intersection_y_equals_1]. 
        */
        if ( 
            ( m_b1 != 0) &&
            ( -m_b0/m_b1 >= x_intersection_y_equals_0) &&
            ( -m_b0/m_b1 <= x_intersection_y_equals_1))
        {
          /* Intersect with the segment x = 0. */
#if OBJ_ON_ARR_DEBUG         
          std::cout<<"Case 1.1.1"<<std::endl;
#endif
          res_bounded_segs[0].set_min(Rational(0),true);
          res_bounded_segs[0].set_max(x_intersection_y_equals_0,true);
                  
          res_bounded_segs[1].set_min(x_intersection_y_equals_1,true);
          res_bounded_segs[1].set_max(Rational(1),true);
          res_splitted = true;
                  
        }
        else
        {
#if OBJ_ON_ARR_DEBUG         
          std::cout<<"Case 1.1.2"<<std::endl;
#endif
          res_bounded_segs[0].set_min(x_intersection_y_equals_0,true);
          res_bounded_segs[0].set_max(x_intersection_y_equals_1,true);
        }
      }
      else /* Case 1.2: The hyperbola intersect first with the 
              segment y = 1.*/
      {
        /* Check if the hyperbola is "splitted" at the unit square,
           if so we will take two arcs, one from
           (min_x,x_intersection_y_equals_1] and the second 
           [x_intersection_y_equals_0,max_x).
           Otherwise take one arc [x_intersection_y_equals_0,
           x_intersection_y_equals_1]. 
        */
        if ( 
            ( m_b1 != 0) &&
            ( -m_b0/m_b1 >= x_intersection_y_equals_1) &&
            ( -m_b0/m_b1 <= x_intersection_y_equals_0))
        {
          /* Intersect with the segment x = 0. */
#if OBJ_ON_ARR_DEBUG         
          std::cout<<"Case 1.2.1"<<std::endl;
#endif
          res_bounded_segs[0].set_min(Rational(0),true);
          res_bounded_segs[0].set_max(x_intersection_y_equals_1,true);
                  
          res_bounded_segs[1].set_min(x_intersection_y_equals_0,true);
          res_bounded_segs[1].set_max(Rational(1),true);
          res_splitted = true;
        }
        else
        {
#if OBJ_ON_ARR_DEBUG         
          std::cout<<"Case 1.2.2"<<std::endl;
#endif
          res_bounded_segs[0].set_min(x_intersection_y_equals_1,true);
          res_bounded_segs[0].set_max(x_intersection_y_equals_0,true);
        }

      }
    }
    /* Case 2: The hyperbola intersects only with the segment y = 1 at the
       positive unit square. */
    else if (0 <= x_intersection_y_equals_1 &&
             x_intersection_y_equals_1 <= 1)
    {
       /* if b0 == 0 => y = infinity */
      bool intersect_with_x_equal_0 = (( m_b0 != 0) && 
                                       ( m_a0/ m_b0 >= 0) &&
                                       ( m_a0/ m_b0 <= 1));
      /* if b0 + m_b1== 0 => y = infinity */
      bool intersect_with_x_equal_1 = (( (m_b0 + m_b1) != 0) && 
                                       ( (m_a0 + m_a1)/(m_b0 + m_b1) >= 0) &&
                                       ( (m_a0 + m_a1)/(m_b0 + m_b1) <= 1));

      /* Case 2.1: The Hyperbola intersect with the segment y = 1 at the
         point x = 0*/
      if (x_intersection_y_equals_1 == 0)
      {
#if OBJ_ON_ARR_DEBUG         
        std::cout<<"Case 2.1"<<std::endl;
#endif
        if (intersect_with_x_equal_1)
        {
          res_bounded_segs[0].set_min(Rational(0),true);
          res_bounded_segs[0].set_max(Rational(1),true);
        }
        else
        {
          res_bounded_segs[0].set_min(Rational(0),true);
          res_bounded_segs[0].set_max(0,true);
        }
      }
            
      /* Case 2.2: The Hyperbola intersect with the segment y = 1 at the
         point x = 1*/
      else if (x_intersection_y_equals_1 == 1)
      {
#if OBJ_ON_ARR_DEBUG         
        std::cout<<"Case 2.2"<<std::endl;
#endif
        if (intersect_with_x_equal_0)
        {
          res_bounded_segs[0].set_min(Rational(0),true);
          res_bounded_segs[0].set_max(Rational(1),true);
        }
        else
        {
          res_bounded_segs[0].set_min(1,true);
          res_bounded_segs[0].set_max(Rational(1),true);
        }
      }
      else if (intersect_with_x_equal_0)
      {
        /* Case 2.3: The hyperbola intersect with the segment x = 0
           at the positive unit sqaure.
           y = (a0 + x * a1)/(m_b1 * x + b0) 
           y = (a0/b0)
        */
#if OBJ_ON_ARR_DEBUG         
        std::cout<<"Case 2.3"<<std::endl;
#endif
         
        res_bounded_segs[0].set_min(Rational(0),true);
        res_bounded_segs[0].set_max(x_intersection_y_equals_1,true);
      }
      else /* Case 2.4: The hyperbola intersect with the segment 
              x = 1 at the positive unit sqaure. */
      {
#if OBJ_ON_ARR_DEBUG         
        std::cout<<"Case 2.4"<<std::endl;
#endif
        res_bounded_segs[0].set_min(x_intersection_y_equals_1,true);
        res_bounded_segs[0].set_max(Rational(1),true);
      }

    }
    /* Case 3: The hyperbola intersects only with the segment y = 0 at
       the positive unit square. */
    else if (0 <= x_intersection_y_equals_0 &&
             x_intersection_y_equals_0 <= 1)
    {
       /* if b0 == 0 => y = infinity */
      bool intersect_with_x_equal_0 = (( m_b0 != 0) && 
                                       ( m_a0/ m_b0 >= 0) &&
                                       ( m_a0/ m_b0 <= 1));

      /* if b0 + b1== 0 => y = infinity */
      bool intersect_with_x_equal_1 = (( (m_b0 + m_b1) != 0) && 
                                       ( (m_a0 + m_a1)/(m_b0 + m_b1) >= 0) &&
                                       ( (m_a0 + m_a1)/(m_b0 + m_b1) <= 1));

      /* Case 3.1: The Hyperbola intersect with the segment y = 1 at
         the point x = 0*/
      if (x_intersection_y_equals_0 == 0)
      {
#if OBJ_ON_ARR_DEBUG         
        std::cout<<"Case 3.1"<<std::endl;
#endif
        if (intersect_with_x_equal_1)
        {
          res_bounded_segs[0].set_min(Rational(0),true);
          res_bounded_segs[0].set_max(Rational(1),true);
        }
        else
        {
          res_bounded_segs[0].set_min(Rational(0),true);
          res_bounded_segs[0].set_max(0,true);
        }
      }
            
      /* Case 2.2: The Hyperbola intersect with the segment y = 1 at
         the point x = 1*/
      else if (x_intersection_y_equals_0 == 1)
      {
#if OBJ_ON_ARR_DEBUG         
        std::cout<<"Case 3.2"<<std::endl;
#endif
        if (intersect_with_x_equal_0)
        {
          res_bounded_segs[0].set_min(Rational(0),true);
          res_bounded_segs[0].set_max(Rational(1),true);
        }
        else
        {
          res_bounded_segs[0].set_min(1,true);
          res_bounded_segs[0].set_max(Rational(1),true);
        }
      }
      else if (intersect_with_x_equal_0)
      {
        /* Case 3.3: The hyperbola intersect with the segment x = 0
           at the positive unit sqaure.*/
#if OBJ_ON_ARR_DEBUG         
        std::cout << "Case 3.3 x_intersection_y_equals_0 = " 
                  << x_intersection_y_equals_0 <<std::endl;
#endif
               
        res_bounded_segs[0].set_min(Rational(0),true);
        res_bounded_segs[0].set_max(x_intersection_y_equals_0,true);
      }
      else /* Case 3.2: The hyperbola intersect with the segment x = 1 at
              the positive unit sqaure. */
      {
#if OBJ_ON_ARR_DEBUG         
        std::cout << "Case 3.4 x_intersection_y_equals_0 = "
                  << x_intersection_y_equals_0<<std::endl;
#endif

        res_bounded_segs[0].set_min(x_intersection_y_equals_0,true);
        res_bounded_segs[0].set_max(Rational(1),true);
      }
    }
    else /* Case 4: The hyperbola doesn't intersect with the segments
            y = 1 and y = 0 at the positive unit square. */
    {
      /* If the hyperbola is below or above the positive unit square 
         return.
         If the hyperbola doesn't intersect with the segment x = 0 in
         the positive unit square => its below or above the unit square.
      */
      if(( m_b0 != 0) &&
         ( m_a0/ m_b0 >= 0) &&
         ( m_a0/ m_b0 <= 1))
      {
#if OBJ_ON_ARR_DEBUG         

        std::cout<<"Case 4.1"<<std::endl;
#endif
        res_bounded_segs[0].set_min(Rational(0),true);
        res_bounded_segs[0].set_max(Rational(1),true);
      }
      else
      {
#if OBJ_ON_ARR_DEBUG         
        std::cout<<"Case 4.2"<<std::endl;
#endif
        /* Set the segment unvalid */
        res_bounded_segs[0].set_min(1,true);
        res_bounded_segs[0].set_max(-1,true);
      }
    }
  }
      

  /*************************************************************
   * Function description:
   * -------------------- 
   * The following functions gets a bounded segment and returns all the arcs
   *  of the hyperbola
   * at this segment. The segment is bounded at y = 0, y = 1, x = 0, x = 1.
   *
   * When y = 1 -> b1 * x + b0 = a0 + x * a1
   * b1 * x - x * a1 = a0 - b0
   * x = (a0 - b0) / (b1 - a1)
   *
   * When y = 0 -> a0 + a1 * x = 0 -> x = -(a0/a1)
   * 
   *
   *
   *************************************************************/
public:
  template <typename Rational_arc_2,typename Isolated_points>
  void get_all_arcs_in_positive_unit_square(std::list<Rational_arc_2>& arcs,
                                            Isolated_points& isolated_points)
  {
    typedef typename Isolated_points::IP_point_and_line_pair 
      Point_and_line_pair;
    typedef typename Point_and_line_pair::PL_Point Point;
    Rational_arc_2 t_arc;
        
#if OBJ_ON_ARR_DEBUG
    //        std::cout << change_color(CGAL_BLUE,*this) << std::endl;
#endif
    if (this->m_obj_on_arr_type == CGAL_DEGENERATE_HYPERBOLA)
    {

           
      /*
        Degenerate hyperbola can't be inserted as regular hyperbola 
        since the vertical segment must be bounded.
      */
      if (m_vertical_segment_x_coordinate <= 1 && 
          m_vertical_segment_x_coordinate >= 0)
      {
        this->get_all_arcs_ver_seg(arcs, isolated_points, 
                                   *m_creator_segment);
      }
            
      if (m_horizontal_segment_y_coordinate <= 1 && 
          m_horizontal_segment_y_coordinate >= 0)
      {
        this->get_all_arcs_horizon_seg(arcs,isolated_points,
                                       *m_creator_segment);
      }
    }
    else if (m_obj_on_arr_type == CGAL_SEGMENT)
    {
      /*   y = (x * a1 + a0) / b0 */
      Rational y_coord[2];
      Rational x_coord[2];
            
      typename Bounded_segs_vector::iterator it;
      for (it = this->m_bounded_segs.begin();
           it != this->m_bounded_segs.end();
           ++it)
      {
        x_coord[0] = (*it).get_min().bound();
        x_coord[1] = (*it).get_max().bound();
               
        y_coord[0] = (m_a1 * x_coord[0] + m_a0)/m_b0;
        y_coord[1] = (m_a1 * x_coord[1] + m_a0)/m_b0;

        /* Check If the segment intersect with the unit square,
         * The segment intersect with y = 0 or y = 1
         * The segment is already bounded by the segments x=1 and x=0.
         */
        if (! ((y_coord[0] < 0 && y_coord[1] < 0) ||
               (y_coord[0] > 1 && y_coord[1] > 1)))
        {
          /* Bound the segemnt to the positive unit square. */
          if (y_coord[0] < y_coord[1])
          {
            if (y_coord[0] < 0)
            {
              y_coord[0] = Rational(0);
              x_coord[0] = - (m_a0/m_a1);
                        
            }
                     
            if (y_coord[1] > 1)
            {
              y_coord[1] = Rational(1);
              x_coord[1] = (m_b0 - m_a0)/m_a1;
            }
          }
          else
          {
            if (y_coord[1] < 0)
            {
              y_coord[1] = Rational(0);
              x_coord[1] = - (m_a0/m_a1);
                        
            }
                     
            if (y_coord[0] > 1)
            {
              y_coord[0] = Rational(1);
              x_coord[0] = (m_b0 - m_a0)/m_a1;
            }
          }
          if (x_coord[0] == x_coord[1] &&
              y_coord[0] == y_coord[1])
          {
            typedef typename Point_and_line_pair::PL_Point Point_2;
            typename Traits_arr_on_plane_2::Point_2 temp_p;
            m_traits_2_adapt.construct_point(temp_p, x_coord[0], y_coord[0]);
            isolated_points.add_element
              (Point_and_line_pair(Point_2(temp_p,
                                           x_coord[0], y_coord[0]),
                                   m_creator_segment));
          }
          else
          {
            m_traits_2_adapt.create_segment_on_plane_arr
              (t_arc,
               Rational_point_2(x_coord[0],y_coord[0]),
               Rational_point_2(x_coord[1],y_coord[1]),
               m_creator_segment);
                     
#if OBJ_ON_ARR_DEBUG
            std::cout << "Arc = " << t_arc << std::endl;
#endif
            arcs.push_back (t_arc);
          }
        }
      }
            
    }
    else /*m_obj_on_arr_type == CGAL_HYPERBOLA*/
    {
      Bounded_seg us_bounded_segs[2];
      bool us_splitted;
       
      Rational x_intersection_y_equals_1 = -1;
      if ((m_b1 - m_a1) != 0)
        x_intersection_y_equals_1 =(m_a0 - m_b0)/(m_b1 - m_a1);
         
      Rational x_intersection_y_equals_0 = -1;
      if (m_a1 != 0)
        x_intersection_y_equals_0 = -(m_a0/m_a1);   
       
#if OBJ_ON_ARR_DEBUG
      std::cout << "x_intersection_y_equals_1 = "
                << x_intersection_y_equals_1 << std::endl;
      std::cout << " x_intersection_y_equals_0 = "
                << x_intersection_y_equals_0 << std::endl;
       
#endif
       
      get_bounded_segments_at_positive_unit_square(x_intersection_y_equals_0,
                                                   x_intersection_y_equals_1,
                                                   us_bounded_segs,
                                                   us_splitted);
         
      Bounded_segs_vector temp_bounded_segs2;
      Bounded_segs_vector temp_bounded_segs(this->m_bounded_segs);
            
      temp_bounded_segs2.add_bounded_seg(Bounded_seg(us_bounded_segs[0]));
      if (us_splitted)
      {
        temp_bounded_segs2.add_bounded_seg(Bounded_seg(us_bounded_segs[1]));
      }
         
      this->m_bounded_segs.merge_bounded_segs_vectors(temp_bounded_segs,
                                                      temp_bounded_segs2);
       
      this->get_all_arcs_of_hyp(arcs, isolated_points, *m_creator_segment);
    }
  }


  template <typename Rational_arc_2,
            typename Isolated_points, 
            typename Ext_obj>
  void get_all_arcs(std::list<Rational_arc_2>& arcs,
                    Isolated_points& isolated_points,
                    const Ext_obj& ext_obj)
  {
    if (this->m_obj_on_arr_type == CGAL_DEGENERATE_HYPERBOLA)
    {
      /*
        Degenerate hyperbola can't be inserted as regular hyperbola 
        since the vertical segment must be bounded.
      */
      this->get_all_arcs_ver_seg(arcs,isolated_points,ext_obj);
      this->get_all_arcs_horizon_seg(arcs,isolated_points,ext_obj);
    }
    else if (m_obj_on_arr_type == CGAL_SEGMENT)
    {
      /*   y = (x * a1 + a0) / b0 */
      Rational y_coord[2];
      Rational x_coord[2];
      Rational_arc_2 t_arc;
            
      typename Bounded_segs_vector::iterator it;
      for (it = this->m_bounded_segs.begin();
           it != this->m_bounded_segs.end();
           ++it)
      {
        x_coord[0] = (*it).get_min().bound();
        x_coord[1] = (*it).get_max().bound();
               
        y_coord[0] = (m_a1 * x_coord[0] + m_a0)/m_b0;
        y_coord[1] = (m_a1 * x_coord[1] + m_a0)/m_b0;

               
        if (x_coord[0] == x_coord[1] &&
            y_coord[0] == y_coord[1])
        {
          typedef typename Isolated_points::IP_point_and_line_pair 
            Point_and_line_pair;
                  
          typedef typename Point_and_line_pair::PL_Point Point_2;
          typename Traits_arr_on_plane_2::Point_2 temp_p;
          m_traits_2_adapt.construct_point(temp_p, x_coord[0], y_coord[0]);

          isolated_points.add_element
            (Point_and_line_pair(Point_2(temp_p,
                                         x_coord[0], y_coord[0]), &ext_obj));
        }
        else
        {
          m_traits_2_adapt.create_segment_on_plane_arr
            (t_arc,
             Rational_point_2(x_coord[0],y_coord[0]),
             Rational_point_2(x_coord[1],y_coord[1]),
             m_creator_segment);
                  
#if OBJ_ON_ARR_DEBUG
          std::cout << "Arc = " << t_arc << std::endl;
#endif
          arcs.push_back (t_arc);
        }
      }
    }
    else /*m_obj_on_arr_type == CGAL_HYPERBOLA*/
    {
           
      if (m_b1 != Rational(0))
      {
        Rbound vertical_asym_x = -m_b0 / m_b1;
        Bounded_segs_vector temp_bounded_segs;
              
        /* Split curves with vertical asymptotes to two curves. */
        typename Bounded_segs_vector::iterator it;
        for (it = this->m_bounded_segs.begin();
             it != this->m_bounded_segs.end();
             ++it)
        {
                 
          if (vertical_asym_x > it->get_min() &&
              vertical_asym_x < it->get_max())
          {
            temp_bounded_segs.add_bounded_seg(Bounded_seg
                                              (vertical_asym_x,
                                               it->get_min(), false,
                                               it->segment_contain_min()));

            temp_bounded_segs.add_bounded_seg(Bounded_seg
                                              (it->get_max(),
                                               vertical_asym_x,
                                               it->segment_contain_max(),
                                               false));
          }
          else
          {
            temp_bounded_segs.add_bounded_seg(*it);
          }
        }
        m_bounded_segs = temp_bounded_segs;
      }
           
      get_all_arcs_of_hyp(arcs, isolated_points, ext_obj);
    }
  }
      
  template <typename Rational_arc_2,
            typename Isolated_points, 
            typename Ext_obj>
  void get_all_arcs_horizon_seg(std::list<Rational_arc_2>& arcs,
                                Isolated_points& isolated_points,
                                const Ext_obj& ext_obj)
  {
    typename Bounded_segs_vector::iterator horizontal_seg_it;
    Rational_arc_2 t_arc;

    for (horizontal_seg_it = m_horizontal_segment_bounded_segs.begin();
         horizontal_seg_it != m_horizontal_segment_bounded_segs.end();
         ++horizontal_seg_it)
    {
      if ((*horizontal_seg_it).get_min() !=
          (*horizontal_seg_it).get_max())
      {
#if OBJ_ON_ARR_DEBUG
        std::cout <<"Push horizontal segment = ((" 
                  << (*horizontal_seg_it).get_min() << "," 
                  << m_horizontal_segment_y_coordinate << "),(" 
                  << (*horizontal_seg_it).get_max() << "," 
                  << m_horizontal_segment_y_coordinate << "))" 
                  << std::endl;
#endif
               

        Rational source;
        Rational target;
        bool dir_right = true;
        bool dir_left = true;
               
        if ((*horizontal_seg_it).get_min().is_bound())
        {
          source = (*horizontal_seg_it).get_min().bound();
          dir_left = false;
        }

        if ((*horizontal_seg_it).get_max().is_bound())
        {
          if (dir_left)
            source = (*horizontal_seg_it).get_max().bound();
          else
            target = (*horizontal_seg_it).get_max().bound();
          dir_right = false;
        }
               
        if (!dir_left && !dir_right)
        {
          m_traits_2_adapt.create_horizontal_curve_on_plane_arr
            (t_arc, m_horizontal_segment_y_coordinate, source, target,m_creator_segment);
        }
        else if (dir_left && dir_right)
        {
          m_traits_2_adapt.create_horizontal_curve_on_plane_arr
             (t_arc, m_horizontal_segment_y_coordinate,m_creator_segment);
        }
        else
        {
          m_traits_2_adapt.create_horizontal_curve_on_plane_arr
            (t_arc, m_horizontal_segment_y_coordinate, source, dir_right,m_creator_segment);
        }
               
        arcs.push_back (t_arc);
      }
      else
      {
        typedef typename Isolated_points::IP_point_and_line_pair 
          Point_and_line_pair;
        typedef typename Point_and_line_pair::PL_Point Point_2;
        typename Traits_arr_on_plane_2::Point_2 temp_p;
        m_traits_2_adapt.construct_point(temp_p, (*horizontal_seg_it).get_min().bound(),
                                         m_horizontal_segment_y_coordinate);

        isolated_points.add_element
           (Point_and_line_pair(Point_2(temp_p,
                                       (*horizontal_seg_it).get_min().bound(),
                                       m_horizontal_segment_y_coordinate),
                               &ext_obj));
      }
    }
  }

  template <typename Rational_arc_2,
            typename Isolated_points, 
            typename Ext_obj>
  void get_all_arcs_ver_seg(std::list<Rational_arc_2>& arcs,
                            Isolated_points& isolated_points,
                            const Ext_obj& ext_obj)
  {
    typename Bounded_segs_vector::iterator vertical_seg_it;
    Rational_arc_2 t_arc;

    for (vertical_seg_it = m_vertical_segment_bounded_segs.begin();
         vertical_seg_it != m_vertical_segment_bounded_segs.end();
         ++vertical_seg_it)
    {
      if ((*vertical_seg_it).get_min() != 
          (*vertical_seg_it).get_max())
      {
#if OBJ_ON_ARR_DEBUG
        std::cout <<"Push vertical segment = (("
                  << m_vertical_segment_x_coordinate << ","
                  << (*vertical_seg_it).get_min() << "),("
                  << m_vertical_segment_x_coordinate << ","
                  << (*vertical_seg_it).get_max() << "))" 
                  << std::endl;
#endif
               
        bool up_valid = false;
        bool down_valid = false;
        Rational_point_2 up;
        Rational_point_2 down;
               
        if ((*vertical_seg_it).get_min().is_bound())
        {
          down = Rational_point_2(m_vertical_segment_x_coordinate,
                                  (*vertical_seg_it).get_min().bound());
          down_valid = true;
        }
               
        if ((*vertical_seg_it).get_max().is_bound())
        {
          up = Rational_point_2(m_vertical_segment_x_coordinate,
                                (*vertical_seg_it).get_max().bound());
          up_valid = true;
        }
        if (up_valid && down_valid)
        {
          m_traits_2_adapt.create_segment_on_plane_arr (t_arc, down, up,
                                                        m_creator_segment);
        }
        else if (up_valid)
        {
          m_traits_2_adapt.create_vertical_segment_on_plane_arr (
             t_arc,
             up,
             false /* Directed down */);
        }
        else if (down_valid)
        {
          m_traits_2_adapt.create_vertical_segment_on_plane_arr (
             t_arc,
             down,
             true /* Directed down */);
        }
        else
        {
                  
          m_traits_2_adapt.create_vertical_segment_on_plane_arr (
             t_arc,
             Rational_point_2(m_vertical_segment_x_coordinate,
                              Rational(0)));
        }
        arcs.push_back (t_arc);
      }
      else
      {
        typedef typename Isolated_points::IP_point_and_line_pair 
          Point_and_line_pair;
        typedef typename Point_and_line_pair::PL_Point Point_2;
        
        typename Traits_arr_on_plane_2::Point_2 temp_p;
        m_traits_2_adapt.construct_point(temp_p, m_vertical_segment_x_coordinate,
                                         (*vertical_seg_it).get_min().bound());

        isolated_points.add_element
          (Point_and_line_pair(Point_2(temp_p,
                                       m_vertical_segment_x_coordinate,
                                       (*vertical_seg_it).get_min().bound()),
                               &ext_obj));
      }
    }
  }
   
  template <typename Rational_arc_2,
            typename Isolated_points, 
            typename Ext_obj>
  void get_all_arcs_of_hyp(std::list<Rational_arc_2>& arcs,
                           Isolated_points& isolated_points,
                           const Ext_obj& ext_obj)
  {
    typedef typename Isolated_points::IP_point_and_line_pair 
                                                        Point_and_line_pair;
    typedef typename Point_and_line_pair::PL_Point      Point_2;

    Rational_arc_2 t_arc;

    Rational        V[4];
    V[0] = m_a0;
    V[1] = m_b0;
    V[2] = m_a1;
    V[3] = m_b1;

    typename Bounded_segs_vector::iterator it;
    for (it = this->m_bounded_segs.begin(); it != this->m_bounded_segs.end();
         ++it)
    {
      if ((*it).get_max() == (*it).get_min())
      {
        Rational x_coordinate = (*it).get_max().bound();
#if OBJ_ON_ARR_DEBUG
        std::cout<<"Add isolated point. x = " << x_coordinate 
                 << " is valid = " << (*it).is_valid() << std::endl;
#endif
        /* Isolated point */
        Rational y_coordinate = (x_coordinate * this->m_a1 + this->m_a0) / 
          (x_coordinate * this->m_b1 + this->m_b0);
        typename Traits_arr_on_plane_2::Point_2 temp_p;
        m_traits_2_adapt.construct_point(temp_p, x_coordinate, y_coordinate);
        
        isolated_points.add_element
          (Point_and_line_pair(Point_2(temp_p,
                                       x_coordinate,y_coordinate),
                               &ext_obj));
      }
      else
      {
        Rational_arc_2 arc;
        if ((*it).get_min().is_bound() && (*it).get_max().is_bound())
        {
          m_traits_2_adapt.create_curve_on_plane_arr(arc,
                                                     (*it).get_min().bound(),
                                                     (*it).get_max().bound(),
                                                     V,
                                                     m_creator_segment);
        }
        else if ((*it).get_max().is_bound())
        {
          m_traits_2_adapt.create_curve_on_plane_arr(arc,
                                                     (*it).get_max().bound(),
                                                     false, V,
                                                     m_creator_segment);
        }
        else if ((*it).get_min().is_bound())
        {
          m_traits_2_adapt.create_curve_on_plane_arr(arc,
                                                     (*it).get_min().bound(),
                                                     true, V,
                                                     m_creator_segment);
        }
        else
        {
          m_traits_2_adapt.create_curve_on_plane_arr(arc, V,
                                                     m_creator_segment);
        }
               
               
#if OBJ_ON_ARR_DEBUG
        std::cout << arc << std::endl;
#endif
        arcs.push_back(arc);
      }
    }
  }

  /*************************************************************
   * Function description:
   * -------------------- 
   * The following functions checks if S1 and S2 are Perpendiculars.
   *
   *
   *************************************************************/
  bool are_perpendiculars(const Rational_segment_3* S1L,
                          const Rational_segment_3* S2L)
  {
    return (((S1L->target().x() - S1L->source().x()) *
             (S2L->target().x() - S2L->source().x()) +
             (S1L->target().y() - S1L->source().y()) *
             (S2L->target().y() - S2L->source().y()) +
             (S1L->target().z() - S1L->source().z()) *
             (S2L->target().z() - S2L->source().z())) == 0);
  }      
      

  /*************************************************************
   * Function description:
   * -------------------- 
   * The following functions gets 2 intersected segments S1 and S3 and
   * computes the degenrate hyperbola.
   *************************************************************/
  void calc_degenerate_hyperbola(const Rational_segment_3* S1,
                                 const Rational_segment_3* S2,
                                 Bounded_segs_vector&
                                   _horizontal_segment_bounded_segs,
                                 Rational& _vertical_segment_x_coordinate,
                                 Rational& _horizontal_segment_y_coordinate,
                                 bool bound_s1_s2)

  {
    Rational_point_3 ipoint_S2;
    Rational_point_3 ipoint_temp;
         
    CGAL::Object result =
      m_rat_kernel->intersect_3_object()(S1->supporting_line(),
                                         m_creator_segment->supporting_line());
        
    CGAL::assign(ipoint_temp, result);
        
    /* only relelvant when the intersection point is on S1 and CL. */
    Rational S1_t = getL_tFromPointOnSegment(*S1,ipoint_temp);
    Rational_plane_3 PlaneS1CL;
        
    if (m_creator_segment->source() != ipoint_temp)
    {
      PlaneS1CL = Rational_plane_3(S1->source(),S1->target(),
                                   m_creator_segment->source());
    }
    else
    {
      PlaneS1CL = Rational_plane_3(S1->source(),S1->target(),
                                   m_creator_segment->target());
    }
        
    if (m_creator_segment->has_on(ipoint_temp))
    {
      _vertical_segment_x_coordinate = S1_t;
    }
    else
    {
      _vertical_segment_x_coordinate = CGAL_SEGMENT_INVALID;
    }

    /* Find the intersection point of S2 with the plane of S1 and S3 */
    result = m_rat_kernel->intersect_3_object()(PlaneS1CL,
                                                S2->supporting_line());
            
    if (!CGAL::assign(ipoint_S2, result) || 
        (bound_s1_s2 && !S2->has_on(ipoint_S2)))
    { 
      /* _S2 is parallel to PlaneS1CL, it's not at the same plane 
         since this condition was  checked at isOneDegreeOfFreedom(). */
           
      _horizontal_segment_y_coordinate = CGAL_SEGMENT_INVALID;
    }
    else
    {
      /* S1_t will be on the plane created by CL and _S2 */
      Rational S2_t = getL_tFromPointOnSegment(*S2,ipoint_S2);
           
      _horizontal_segment_y_coordinate = S2_t;
           
      /* 
         Find the bounds of _S1.
         In order to find these bounds send segments from the 
         intersection of _S2 with the plane created by
         _S1 and CL, if the segmnet of _S2 intersect with the plane,
         send segments from the intersction point to the ends of
         the segments of _S1 and CL.
         The bounds will be between the segments that intersect 
         with both the segments of _S1 and CL.
      */
           
      if (m_creator_segment->has_on(ipoint_S2))
      {
        if (bound_s1_s2)
        {
           _horizontal_segment_bounded_segs.add_bounded_seg
              (Bounded_seg(Rbound(Rational(1)),
                           Rbound(Rational(0)),
                           true, true));
        }
        else
        {
          _horizontal_segment_bounded_segs.add_bounded_seg
            (Bounded_seg
             (LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY,
              LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY,
              false, false));
                      
        }
        return;
      }
           
      m_arr_g_func.get_bounded_segments_of_point_and_2_lines(
         m_rat_kernel,
         *S1,
         *m_creator_segment,
         m_bound_s1_s2,
         ipoint_S2,
         _horizontal_segment_bounded_segs);
    }
  }

   
  /*************************************************************
   * Function description:
   * -------------------- 
   * The following function gets the following 2 degenerate equations (t = C)
   * of 3 segments, and Returns Hyperbola/Segment that represetns all the
   * segments that passes through these 3 segments.
   * 
   *
   * 1. S1_t * S1_x + S2_t * S2_x = S3_t * S3_x + C_x
   * 2. S1_t * S1_y + S2_t * S2_y = S3_t * S3_y + C_y
   *
   *************************************************************/
  int calc_obj_on_arr_t_is_constant(Rational  S1_x,
                                    Rational  S2_x,
                                    Rational  S3_x,
                                    Rational  C_x,
                                    Rational  S1_y,
                                    Rational  S2_y,
                                    Rational  S3_y,
                                    Rational  C_y)
  {
#if OBJ_ON_ARR_DEBUG    
    std::cout << "S1_t * " << S1_x << " + S2_t * " << S2_x << " = S3_t * "
              << S3_x << " + " << C_x << std::endl;
    std::cout << "S1_t * " << S1_y << " + S2_t * " << S2_y << " = S3_t * "
              << S3_y << " + " << C_y << std::endl;
#endif
    if (S3_y == 0)
    {
      std::swap(S1_x,S1_y);
      std::swap(S2_x,S2_y);
      std::swap(S3_x,S3_y);
      std::swap(C_x,C_y);
    }

    if (S3_x == 0)
    {
      /* The object is a segment at the arrangement. */
      /*  S2_t =  (C_x - S1_t * S1_x)/S2_x */
      CGAL_assertion(S2_x != 0);
      {
        /* S2_t =  (C_x - S1_t * S1_x)/S2_x */
        /*   y = (x * m_a1 + a0) / (x * b1 + m_b0) */
        m_a1 = - (S1_x / S2_x);
        m_a0 = C_x / S2_x;
        m_b0 = 1;
        m_b1 = 0;
        m_obj_on_arr_type = CGAL_SEGMENT;
#if OBJ_ON_ARR_DEBUG
        std::cout << "S3_x == 0" << std::endl;
        std::cout << *this << std::endl;
#endif
        /* S3_t does not make any difference since its coefficients are 
           equal to 0. Compute bounded segment over S1_t */
        if (S3_y == 0)
        {
          if (m_bound_s1_s2)
          {
            m_bounded_segs.add_bounded_seg
              (Bounded_seg(Rbound(Rational(1)),
                           Rbound(Rational(0)),
                           true, true));
          }
          else
          {
            m_bounded_segs.add_bounded_seg
              (Bounded_seg
               (LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY,
                LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY,
                false, false));

          }
        }
        else
        {
          calc_bounded_segments_t_is_constant(S1_y,S2_y,S3_y,C_y);
        }
      }
    }
    else
    {
      /*
       * 1. S1_t * S1_x + S2_t * S2_x = S3_t * S3_x + C_x
       * 2. S1_t * S1_y + S2_t * S2_y = S3_t * S3_y + C_y
       *
       * 1. S3_t = (S1_t * S1_x + S2_t * S2_x - C_x) / S3_x.
       * 2. S3_t = (S1_t * S1_y + S2_t * S2_y - C_y) / S3_y
       *
       *   (S1_t * S1_x + S2_t * S2_x - C_x) / S3_x = (S1_t * S1_y + S2_t *
       *         S2_y - C_y)/S3_y
       *
       *   S1_t * S1_x * S3_y + S2_t * S2_x * S3_y - C_x * S3_y = 
       *          S1_t * S1_y * S3_x + S2_t * S2_y * S3_x - C_y * S3_x
       *
       *   S1_t * (S1_x * S3_y - S1_y * S3_x) + C_y * S3_x - C_x * S3_y = 
       *          S2_t * (S2_y * S3_x - S2_x * S3_y)  
       *
       *   Let B1 = S1_x * S3_y - S1_y * S3_x
       *
       *   Let B2 = C_y * S3_x - C_x * S3_y
       *
       *   Let B3 = S2_y * S3_x - S2_x * S3_y
       *
       *   S1_t * B1 + B2 = S2_t * B3
       * 
       */

      Rational B1 = S1_x * S3_y - S1_y * S3_x;
      Rational B2 = C_y * S3_x - C_x * S3_y;
      Rational B3 = S2_y * S3_x - S2_x * S3_y;
            
      if (B3 == 0)
         return LTS_g_func::CGAL_QUERY_FAILD;
      {
        /*   y = (x * m_a1 + a0) / (x * b1 + b0) */
        /*   S1_t * B1 + B2 = S2_t * B3
         *
         *   S2_t = (S1_t * B1 + B2)/B3
         */
        m_a1 = B1;
        m_a0 = B2;
        m_b0 = B3;
        m_b1 = 0;
        m_obj_on_arr_type = CGAL_SEGMENT;
               
#if OBJ_ON_ARR_DEBUG
        std::cout << "S3_x != 0" << std::endl;
        std::cout << "S1_t * " << B1 <<" + " << B2 << "=  S2_t *" << B3
                  << std::endl;
        std::cout << *this << std::endl;
#endif

        calc_bounded_segments_t_is_constant(S1_y,S2_y,S3_y,C_y);
      }
    }
    return LTS_g_func::CGAL_QUERY_SUCCEED;
  }
            
  /*************************************************************
   * Function description:
   * -------------------- 
   * The following function gets the following 2 degenerate equations
   *   (t = C) of 3 segments,and
   * returns the bouneded segment of the ruling of 0 <= S3_t <=1. 
   * 
   * S1_t * S1_y + S2_t * S2_y = S3_t * S3_y + C_y
   * 
   * 0 <= S3_t <= 1 
   *
   * S1_t * S1_y + S2_t * S2_y = S3_t * S3_y + C_y.
   *
   * S3_t * S3_y = S1_t * S1_y + S2_t * S2_y - C_y.
   *
   * S2_t = ((S1_t * a1 + a0) / b0)
   *
   * S3_t = (S1_t * S1_y + ((S1_t * a1 + a0) / b0) * S2_y - C_y) / S3_y.
   *
   * 1 >= (S1_t * S1_y + ((S1_t * a1 + a0) / b0) * S2_y - C_y) / S3_y >= 0.
   *
   * 1 >= (S1_t * S1_y * b0 + S1_t * a1 * S2_y + a0 * S2_y - b0 * C_y) /
   *     (b0 * S3_y) >= 0.
   *
   * 1 >= (S1_t * (S1_y * b0 + a1 * S2_y) + a0 * S2_y - b0 * C_y) / 
   *    (b0 * S3_y) >= 0.
   *
   * D1 = S1_y * b0 + a1 * S2_y
   *
   * D2 = a0 * S2_y - b0 * C_y
   *
   * D3 = S3_y * b0.
   *
   * 1 >= (S1_t * D1 + D2) / D3 >= 0.
   *
   * 1 - D2 / D3 >= (S1_t * D1 / D3)  >= - D2 / D3.
   *
   *
   *************************************************************/
  void calc_bounded_segments_t_is_constant(const Rational&  S1_y,
                                           const Rational&  S2_y,
                                           const Rational&  S3_y,
                                           const Rational&  C_y)
  {
    Rational D1 = S1_y * m_b0 + m_a1 * S2_y;
    Rational D2 = m_a0 * S2_y - C_y * m_b0;
    Rational D3 = S3_y * m_b0;

#if OBJ_ON_ARR_DEBUG                 
    std::cout << "D1 = " << D1 << " D2 = " << D2 << " D3 = " << D3 
              << std::endl;
#endif
    Rational min_x;
    Rational max_x;
    if ((D3/D1) < 0)
    {
      /* 1 - D2 / D3 <= (S1_t * D1 / D3)  <= - D2 / D3. */
      /* (D3 - D2) / D1 <= S1_t <= - D2 / D1 */
      if (m_bound_s1_s2)
      {
        min_x = max(Rational(0), ((D3 - D2) / D1));
        max_x = min(Rational(1), (-D2/D1));
              
      }
      else
      {
        min_x = ((D3 - D2) / D1);
        max_x = (-D2/D1);
      }    
    }
    else
    {
      /* (D3 - D2) / D1 >= S1_t >= - D2 / D1 */
      if (m_bound_s1_s2)
      {
        min_x = max(Rational(0),(-D2/D1));
        max_x = min(Rational(1),((D3 - D2) / D1));
      }
      else
      {
        min_x = (-D2/D1);
        max_x = ((D3 - D2) / D1);
      }
    }
    m_bounded_segs.add_bounded_seg(Bounded_seg(max_x, min_x, true, true));
  }
      
  /*************************************************************
   * Function description:
   * -------------------- 
   * The following functions gets a Segment L=p1+L_t*v and a point p2 on the
   * segment.
   * The function return the related place on the segment from the point p1.
   * L_t = p2-p1/v
   *
   *
   *************************************************************/
private:
  Rational getL_tFromPointOnSegment(const Rational_segment_3& l,
                                    Rational_point_3& p2)
  {
    if (l.target().x() != l.source().x())
    {
      return (p2.x() - l.source().x())/(l.target().x() - l.source().x());
    }
    else if (l.target().y() != l.source().y())
    {
      return (p2.y() - l.source().y())/(l.target().y() - l.source().y());
    }
    else if (l.target().z() != l.source().z())
    {
      return (p2.z() - l.source().z())/(l.target().z() - l.source().z());
    }
    CGAL_error();
    return Rational();
  }
};
  
} //namespace CGAL

#endif /* LINE_THROUGH_SEGMENTS_ARR_OBJECT_H */
