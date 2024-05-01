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

#ifndef LINES_THROUGH_SEGMENTS_IMPL_H
#define LINES_THROUGH_SEGMENTS_IMPL_H

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>

#include <CGAL/basic.h>
#include <CGAL/Lines_through_segments_arr_object.h>
#include <CGAL/Lines_through_segments_general_functions.h>
#include <CGAL/Lines_through_segments_arr_gen_func.h>
#include <CGAL/Lines_through_segments_isolated_points.h>
#include <CGAL/Lines_through_segments_3/observer.h>
#include <CGAL/Lines_through_segments_arr_ext_dcel.h>
#include <CGAL/Lines_through_segments_arr_plane_faces.h>
#include <CGAL/Lines_through_segments_3_segs_2.h>
#include <CGAL/Lines_through_segments_3/internal.h>
#include <CGAL/Lines_through_segments_3/exceptions.h>
#include <CGAL/Lines_through_segments_point_adapt.h>
#include <CGAL/Lines_through_segments_traits_2_adapt.h>
#include <CGAL/Lines_through_segments_output_obj.h>
#if LTS_DRAW_ARR
#include "Arrangement_general_functions.h"
#endif

#define ONE_SWEEP 1

/*************************************************************
 * This class is the core of the implementation, it holds arrangements on
 * surface at the following types:
 *
 * 1. Arrangement of segment arcs of hyperbolas on the plane.
 * 2. Arrangement of geodesic arcs on sphere.
 *
 *************************************************************/

namespace CGAL {

template <typename Lines_through_segments_traits_3,
          typename With_segments>
class Lines_through_segments_impl {
public:
  using Traits_3 = Lines_through_segments_traits_3;

private:
  using Alg_kernel = typename Traits_3::Alg_kernel;
  using Rational_kernel = typename Traits_3::Rational_kernel;
  using Algebraic = typename Alg_kernel::FT;
  using Rational = typename Rational_kernel::FT;
  using Traits_arr_on_plane_2 = typename Traits_3::Traits_arr_on_plane_2;
  using Traits_arr_on_sphere_2 = typename Traits_3::Traits_arr_on_sphere_2;

  using Alg_point_3 = typename Alg_kernel::Point_3;
  using Alg_line_3 = typename Alg_kernel::Line_3;
  using Alg_segment_3 = typename Alg_kernel::Segment_3;
  using Alg_plane_3 = typename Alg_kernel::Plane_3;
  using Point_on_plane_2 = typename Traits_arr_on_plane_2::Point_2;
  using Point_on_sphere_2 = typename Traits_arr_on_sphere_2::Point_2;

  using Rational_point_3 = typename Rational_kernel::Point_3;
  using Rational_line_3 = typename Rational_kernel::Line_3;
  using Rational_segment_3 = typename Rational_kernel::Segment_3;
  using Rational_plane_3 = typename Rational_kernel::Plane_3;
  using Rational_point_2 = typename Rational_kernel::Point_2;
  using Rational_direction_3 = typename Rational_kernel::Direction_3;

  /***************************************************/
  /*    Arrangement on plane typedefs.               */
  /***************************************************/

  /* Extended each edge with its creator line segment.*/
  using Dcel_on_plane =
    Lines_through_segments_arr_ext_dcel<Traits_arr_on_plane_2,
                                        Rational_segment_3>;
  using Arrangement_on_plane_2 =
    CGAL::Arrangement_2<Traits_arr_on_plane_2, Dcel_on_plane>;

  using Point_2 =
    Lines_through_segments_point_adapt_2<Traits_3, Point_on_plane_2,
                                         Algebraic>;
  using Rational_arc_2 = typename Traits_arr_on_plane_2::Curve_2;

  using Lines_through_segments_arr_observer_on_plane =
    Lines_through_segments_arr_observer<Arrangement_on_plane_2>;

  using Arr_object =
    Lines_through_segments_arr_object<Traits_3, With_segments>;

  /***************************************************/
  /*    Arrangement on sphere typedefs.              */
  /***************************************************/

  using Dcel_geom_traits =
    Lines_through_segments_arr_ext_dcel<Traits_arr_on_sphere_2,
                                        Rational_segment_3>;
  using Topol_traits_2 =
    CGAL::Arr_spherical_topology_traits_2<Traits_arr_on_sphere_2,
                                          Dcel_geom_traits>;

  using Adapted_point_on_sphere_2 =
    Lines_through_segments_point_adapt_3<Traits_3, Point_on_sphere_2,
                                         Rational>;

  using X_monotone_curve_on_Sphere_2 =
    typename Traits_arr_on_sphere_2::X_monotone_curve_2;
  using Curve_2 = typename Traits_arr_on_sphere_2::Curve_2;
  using Curve_2_no_data = typename Traits_arr_on_sphere_2::Base_curve_2;
  using Location_type = typename Adapted_point_on_sphere_2::Location_type;

  using Arrangement_on_sphere_2 =
    CGAL::Arrangement_on_surface_2<Traits_arr_on_sphere_2, Topol_traits_2>;

  using Lines_through_segments_arr_plane_faces =
    class Lines_through_segments_arr_plane_faces
    <Arrangement_on_plane_2, Lines_through_segments_arr_observer_on_plane>;

  using Isolated_points_on_line =
    Lines_through_segments_isolated_points
    <Point_and_segment_pair<Rational_point_3, Rational_segment_3>,
     Compare_points_on_line<Rational_point_3, Rational_segment_3>>;

  using Isolated_points_on_plane =
    Lines_through_segments_isolated_points
    <Point_and_segment_pair<Point_2, Rational_segment_3>,
     Compare_points_on_plane<Point_2, Rational_segment_3>>;

  using Isolated_points_on_sphere =
    Lines_through_segments_isolated_points
    <Point_and_segment_pair<Adapted_point_on_sphere_2, Rational_segment_3>,
     Compare_points_on_sphere <Adapted_point_on_sphere_2, Rational_segment_3>>;

  using Point_on_plane_and_line_pair =
    Point_and_segment_pair<Point_2, Rational_segment_3>;

  using Point_on_sphere_and_line_pair =
    Point_and_segment_pair<Adapted_point_on_sphere_2, Rational_segment_3>;

  using Traits_2_adapt =
    Lines_through_segments_traits_on_plane_adapt<Traits_3>;

  using Through_3 =
    typename Lines_through_segments_output_obj<Traits_3, Rational_segment_3>::Through_3;

public:
  using Lines_through_segments_arr_observer_on_sphere =
    Lines_through_segments_arr_observer<Arrangement_on_sphere_2>;

  //!
  class Arc_end_points {
  private:
    Rational_point_2 source_p;
    Rational_point_2 target_p;

  public:
    Arc_end_points() {}

    Arc_end_points(Rational_point_2& _source, Rational_point_2& _target) {
      source_p = _source;
      target_p = _target;
    }

    Arc_end_points(Rational_point_2 _source, Rational_point_2 _target) {
      source_p = _source;
      target_p = _target;
    }

    Rational_point_2 source() { return source_p; }

    Rational_point_2 target() { return target_p; }

    Rational_point_2& left() {
      if (source_p.x() <= target_p.x()) return source_p;
      else return target_p;
    }

    Rational_point_2& right() {
      if (source_p.x() > target_p.x()) return source_p;
      else return target_p;
    }

    Rational_point_2 top() {
      if (source_p.y() <= target_p.y()) return source_p;
      else return target_p;
    }

    Rational_point_2 bottom() {
      if (source_p.y() > target_p.y()) return source_p;
      else return target_p;
    }
  };

private:
  Lines_through_segments_arr_gen_func<Traits_3, With_segments> m_arr_g_func;
  const Alg_kernel* m_alg_kernel;
  const Rational_kernel* m_rat_kernel;
  boost::shared_ptr<Arrangement_on_plane_2> m_arr_on_plane;

  //! On sphere arrangement
  boost::shared_ptr<Arrangement_on_sphere_2> m_arr_on_sphere;
  const Traits_arr_on_sphere_2* m_traits_on_sphere;

  std::list<Rational_arc_2> arcs_cont;

  Lines_through_segments_arr_observer_on_plane* m_obs_on_plane;
  Lines_through_segments_arr_observer_on_sphere* m_obs_on_sphere;
  Isolated_points_on_plane m_isolated_points_on_plane;
  Isolated_points_on_sphere m_isolated_points_on_sphere;
  const Rational_segment_3* m_S1;
  const Rational_segment_3* m_S2;
  bool m_S1_S2_intersect;
  bool m_S1_S2_on_the_same_plane;
  /* If S1 or S2 are points the arrangement is 1 dimensional. */
  bool m_one_dimensional_arr;
  bool m_S1_S2_on_the_same_line;
  const Rational_segment_3* m_S1_S2_on_the_same_line_S3;
  const Rational_segment_3* m_S1_S2_on_the_same_line_S4;
  bool m_S1_is_a_point;
  bool m_S2_is_a_point;
  Rational_line_3 m_S1_S2_line;
  Rational_plane_3 m_S1_S2_plane;
  Rational_point_3 m_intersection_point_S1S2;
  int m_S1_S2_num_of_common_lines;
  Traits_2_adapt m_traits_2_adapt;
  Lines_through_segments_general_functions<Traits_3>
  m_g_func;

  /* Holds the number of segments that passes through
     m_intersection_point_S1S2  except of S1 and S2. */
  unsigned int m_num_of_intersection;
  const Rational_segment_3* m_S1_S2_concurrent_S3;
  const Rational_segment_3* m_S1_S2_concurrent_S4;

  /* Holds arrangements that were created from 3 lines on the same
   * plane that a plane passes through them.
   * After all lines are inserted, a overlay of all of this arrangement
   * is done.
   * If two faces a overlap the created face represent a plane that passes
   * through 4 lines.
   * All other faces represent a plane that passes through 3 lines.
   * After this process was completed the overlaid arrangement is overlayed
   * with arr_on_plane,
   * Each edge that falls inside the faces of the overlaid arrangement
   * represnt a plane that passes through 4 lines.
   */
  Lines_through_segments_arr_plane_faces m_arr_plane_faces;

public:
  ~Lines_through_segments_impl() {
    if (m_S1_S2_intersect) delete m_obs_on_sphere;
    delete m_obs_on_plane;
  }

  Lines_through_segments_impl(const Rational_segment_3* _s1,
                              const Rational_segment_3* _s2,
                              const Alg_kernel* alg_kernel,
                              const Rational_kernel* rational_kernel) :
    m_arr_g_func(alg_kernel),
    m_alg_kernel(alg_kernel),
    m_rat_kernel(rational_kernel),
    m_arr_on_plane(new Arrangement_on_plane_2()),
    m_arr_on_sphere(new Arrangement_on_sphere_2())
  {
    m_traits_on_sphere = m_arr_on_sphere->geometry_traits();

    m_S1 = _s1;
    m_S2 = _s2;

    m_S1_S2_intersect = false;
    m_num_of_intersection = 0;
    m_S1_S2_on_the_same_plane = false;
    m_S1_S2_on_the_same_line = false;
    m_one_dimensional_arr = false;
    m_S2_is_a_point = m_S2->source() == m_S2->target();
    m_S1_is_a_point = m_S1->source() == m_S1->target();
    m_S1_S2_num_of_common_lines = 0;
    Compare_points_on_plane<Point_2, Rational_segment_3> cp;
    m_isolated_points_on_plane = Isolated_points_on_plane(cp);
    Compare_points_on_sphere<Adapted_point_on_sphere_2, Rational_segment_3> cs;
    m_isolated_points_on_sphere = Isolated_points_on_sphere(cs);

    /* m_S1 is a point */
    if (m_S1_is_a_point) {
#if ARR_ON_SUR_DEBUG
      std::cout << "S1 is a point" << std::endl;
#endif
      line_is_a_point(m_S1->source(),*m_S2);
    }
    else if (m_S2_is_a_point) {
#if ARR_ON_SUR_DEBUG
      std::cout << "S2 is a point" << std::endl;
#endif
      /* In case S2 is a point swap m_S1 and S2. */
      std::swap(m_S1,m_S2);
      m_S1_is_a_point = true;
      m_S2_is_a_point = false;
      line_is_a_point(m_S1->source(),*m_S2);
    }
    else if (m_rat_kernel->do_intersect_3_object()(m_S1->supporting_line(),
                                                   m_S2->supporting_line())) {
      m_S1_S2_on_the_same_plane = true;
      if (m_S1->supporting_line().has_on(m_S2->source()))
        m_S1_S2_plane = Rational_plane_3(m_S1->supporting_line(),
                                         m_S2->target());
      else
        m_S1_S2_plane = Rational_plane_3(m_S1->supporting_line(),
                                         m_S2->source());


      CGAL::Object result =
        m_rat_kernel->intersect_3_object()(m_S1->supporting_line(),
                                           m_S2->supporting_line());
      if (CGAL::assign(m_intersection_point_S1S2, result)) {
        if (m_S1->has_on(m_intersection_point_S1S2) &&
            m_S2->has_on(m_intersection_point_S1S2))
        {
#if ARR_ON_SUR_DEBUG
          std::cout <<"ARR ON SPHERE" << std::endl;
#endif
          m_S1_S2_intersect = true;

          m_obs_on_sphere =
            new Lines_through_segments_arr_observer_on_sphere(*m_arr_on_sphere);
#if ARR_ON_SUR_DEBUG
          std::cout << "S1,S2 intersect" << std::endl;
          std::cout << "m_intersection_point_S1S2 = "
                    << m_intersection_point_S1S2 << std::endl;
#endif
        }
        /* Create Arrangement on the plane from the two lines and
         * the points of the segments that intersect with Plane.
         */
      }
      else {
        m_one_dimensional_arr = true;
        m_S1_S2_on_the_same_line = true;
        m_S1_S2_line = m_S1->supporting_line();

        if (are_overlap(*m_S1,*m_S2))
          throw Lines_through_segments_exp_2_lines_overlap();
      }
    }
    else if (m_rat_kernel->are_parallel_3_object()(m_S1->supporting_line(),
                                                   m_S2->supporting_line()))
    {
      m_S1_S2_on_the_same_plane = true;
      if (m_S1->supporting_line().has_on(m_S2->source()))
        m_S1_S2_plane = Rational_plane_3(m_S1->supporting_line(),
                                         m_S2->target());
      else
        m_S1_S2_plane = Rational_plane_3(m_S1->supporting_line(),
                                         m_S2->source());
    }

#if ARR_ON_SUR_DEBUG
    std::cout <<"ARR ON PLANE" << std::endl;
#endif
    m_obs_on_plane =
      new Lines_through_segments_arr_observer_on_plane(*m_arr_on_plane);
  }

public:
  /*************************************************************
   * The following function adds element to the arrangement.
   * Precondition:
   *        S1 and S2 are on the same plane.
   *
   * Input:
   *      S3
   *************************************************************/

  void add_element_to_arrangement_S1S2_on_the_same_plane
  (const Rational_segment_3& S3) {
    if (S3.source() == S3.target()) {
      if (m_S1_S2_plane.has_on(S3.source())) {
        add_element_to_plane_arrangement(S3.source(), S3);
      }

      if (m_S1_S2_intersect) {
        Rational p0_x_diff = (S3.source().x() - m_intersection_point_S1S2.x());
        Rational p0_y_diff = (S3.source().y() - m_intersection_point_S1S2.y());
        Rational p0_z_diff = (S3.source().z() - m_intersection_point_S1S2.z());

        /* Check if the line passes through the point
           m_intersection_point_S1S2 */
        if (S3.has_on(m_intersection_point_S1S2)) {
#if ARR_ON_SUR_DEBUG
          std::cout << change_color(CGAL_YELLOW,
                                    "num of intersection point ++")
                    << std::endl;
#endif
          if (m_num_of_intersection == 0) m_S1_S2_concurrent_S3 = &S3;
          else if (m_num_of_intersection == 1) m_S1_S2_concurrent_S4 = &S3;
          m_num_of_intersection += 1;

        }
        else {
          auto loc = location(p0_x_diff, p0_y_diff, p0_z_diff);
          Adapted_point_on_sphere_2 pos(p0_x_diff, p0_y_diff, p0_z_diff, loc);
          m_isolated_points_on_sphere.add_element
            (Point_on_sphere_and_line_pair(pos, &S3));
        }
      }
    }
    else {
      CGAL::Object result =
        m_rat_kernel->intersect_3_object()(m_S1_S2_plane,
                                           S3.supporting_line());

      /* Find the intersection of S3 with the Plane of S1 and S2.*/
      Rational_point_3 ipoint;
      Rational_line_3 iline;
      if (CGAL::assign(ipoint, result)) {
        if (S3.has_on(ipoint) &&
            /* Sphere Case - S3 is crossing the intersection point of
               S1 and S2.*/
            !(m_S1_S2_intersect && S3.has_on(m_intersection_point_S1S2)) &&
            /* Sphere Case - Point on Sphere but not on S1.*/
            !(m_S1_S2_intersect &&
              !m_S1->has_on(ipoint) &&
              m_S1->supporting_line().has_on(ipoint)) &&
            /* Sphere Case - Point on Sphere but not on S1.*/
            !(m_S1_S2_intersect &&
              !m_S2->has_on(ipoint) &&
              m_S2->supporting_line().has_on(ipoint)))
        {
          add_element_to_plane_arrangement(ipoint,S3);
        }
      }
      else if (CGAL::assign(iline, result)) {
#if ARR_ON_SUR_DEBUG
        std::cout << "m_S1 = " << *m_S1 << std::endl;
        std::cout << "S2 = " << *m_S2 << std::endl;
        std::cout << "S3 = " << S3 << std::endl;
#endif
        if (!(m_S1_S2_intersect &&
              (m_g_func.has_the_same_supporting_but_not_intersect(*m_S1,S3) ||
               m_g_func.has_the_same_supporting_but_not_intersect(*m_S2,S3))))
        {
          add_element_to_plane_arrangement(S3);
        }
      }
      if (m_S1_S2_intersect) {
        /* Check if the line passes through the point
           m_intersection_point_S1S2 */
        if (S3.has_on(m_intersection_point_S1S2)) {
#if ARR_ON_SUR_DEBUG
          std::cout << change_color(CGAL_YELLOW,
                                    "num of intersection point ++")
                    << std::endl;
#endif
          m_num_of_intersection += 1;
          if (m_num_of_intersection == 0) m_S1_S2_concurrent_S3 = &S3;
          else if (m_num_of_intersection == 1) m_S1_S2_concurrent_S4 = &S3;
        }
        else add_element_to_sphere_arrangement(S3);
      }
    }
  }
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
  void add_element_to_plane_arrangement(const Rational_segment_3& S3) {
    if (! m_S1_S2_intersect &&
        m_S1->supporting_line().has_on(S3.source()) &&
        m_S1->supporting_line().has_on(S3.target())) {
      Rational_point_3 ipoint;
      Rational S2_t;
      CGAL::Object result =
        m_rat_kernel->intersect_3_object()(m_S1->supporting_line(),
                                           m_S2->supporting_line());

      if (CGAL::assign(ipoint, result) && m_S2->has_on(ipoint)) {
        m_g_func.get_scalar_from_point_on_line(ipoint,
                                               m_S2->supporting_line(),S2_t);

        Point_on_plane_2 temp_p;
        m_traits_2_adapt.construct_point(temp_p, Rational(0), S2_t);

        m_isolated_points_on_plane.add_element
          (Point_on_plane_and_line_pair(Point_2(temp_p,
                                                Rational(0), S2_t),
                                        &S3));
      }
    }
    else if (!m_S1_S2_intersect &&
             m_S2->supporting_line().has_on(S3.source()) &&
             m_S2->supporting_line().has_on(S3.target())) {
      Rational_point_3 ipoint;
      Rational S1_t;
      CGAL::Object result =
        m_rat_kernel->intersect_3_object()(m_S1->supporting_line(),
                                           m_S2->supporting_line());

      if (CGAL::assign(ipoint, result) && m_S1->has_on(ipoint)) {
        m_g_func.get_scalar_from_point_on_line(ipoint,
                                               m_S1->supporting_line(),S1_t);

        Point_on_plane_2 temp_p;
        m_traits_2_adapt.construct_point(temp_p, S1_t, Rational(0));

        m_isolated_points_on_plane.add_element
          (Point_on_plane_and_line_pair(Point_2(temp_p,
                                                S1_t, Rational(0)),
                                        &S3));
      }
    }
    else {
      using Lines_through_segments_3_segs_2 =
        Lines_through_segments_3_segs_2<Traits_3,
                                        Lines_through_segments_arr_plane_faces,
                                        Isolated_points_on_plane,
                                        Point_on_plane_and_line_pair,
                                        Arc_end_points,
                                        With_segments>;

      Lines_through_segments_3_segs_2 lines_through_segments_3_segs;

      lines_through_segments_3_segs(m_S1,m_S2,S3,true /* Bound S1 S2 */,
                                    m_rat_kernel,
                                    m_alg_kernel,
                                    m_intersection_point_S1S2,
                                    m_S1_S2_intersect,
                                    m_isolated_points_on_plane,
                                    m_arr_plane_faces);
    }
  }

  /*************************************************************
   * The following function adds isolated points to one dimensional
   * arrangement.
   *
   * Input:
   *      S3
   *************************************************************/

  void
  add_element_to_one_dimensional_arrangement(const Rational_segment_3& S3) {
    m_obs_on_plane->set_is_plane(false);
    /* Get the intersection point of */
    CGAL::Object result =
      m_rat_kernel->intersect_3_object()(S3.supporting_line(), m_S1_S2_plane);

    Rational_point_3 qpoint;
    Rational_line_3 qline;
    if (S3.source() == S3.target() && m_S1_S2_plane.has_on(S3.source())) {
      Rational_line_3 S3_S1_line(m_S1->source(),S3.source());
      Rational_point_3 ipoint;
      CGAL::Object result =
        m_rat_kernel->intersect_3_object()(m_S2->supporting_line(), S3_S1_line);
      if (CGAL::assign(ipoint, result)) {
        if (m_S2->has_on(ipoint)) {
          Rational S2_t;
          m_g_func.get_scalar_from_point_on_line(ipoint, m_S2->supporting_line(),
                                                 S2_t);

          Point_on_plane_2 temp_p;
          m_traits_2_adapt.construct_point(temp_p, Rational(0), S2_t);

          m_isolated_points_on_plane.add_element
            (Point_on_plane_and_line_pair(Point_2(temp_p, Rational(0), S2_t),
                                          &S3));
        }
      }
    }
    else if (CGAL::assign(qpoint, result)) {
      if (S3.has_on(qpoint)) {
        Rational_line_3 line_qpoint_S1(m_S1->source(),qpoint);
        result = m_rat_kernel->intersect_3_object()(m_S2->supporting_line(),
                                                    line_qpoint_S1);
        if (CGAL::assign(qpoint, result)) {
          if (m_S2->has_on(qpoint)) {
            Rational S2_t;
            m_g_func.get_scalar_from_point_on_line(qpoint,
                                                   m_S2->supporting_line(),S2_t);

            Point_on_plane_2 temp_p;
            m_traits_2_adapt.construct_point(temp_p, Rational(0), S2_t);

            m_isolated_points_on_plane.add_element
              (Point_on_plane_and_line_pair(Point_2(temp_p,
                                                    Rational(0),S2_t), &S3));
          }
        }
      }
    }
    else if (CGAL::assign(qline, result)) {
      Lines_through_segments_bounded_segs_vector<Rational>
        ret_bounded_segments;
      Rational_point_3 P1;
#if ARR_ON_SUR_DEBUG
      std::cout <<"qline = " << qline << std::endl;
#endif
      std::list<Rational_arc_2> arcs;
      m_arr_g_func.add_segs_to_arr_S1_is_a_point(*m_S2,
                                                 S3,
                                                 true /* Bound S2 */,
                                                 m_S1->source(),
                                                 arcs,
                                                 m_rat_kernel,
                                                 m_isolated_points_on_plane,
                                                 Rational(0),
                                                 m_arr_on_plane.get(),
                                                 S3);
    }
  }

  /*************************************************************
   * Function description:
   * --------------------
   * The following functions checks if 2 of the segments are on the same
   * supporting line and that the segments overlap.
   *
   *************************************************************/
  bool are_overlap(const Rational_segment_3& seg1,
                   const Rational_segment_3& seg2)
  {
    /* One of the segments is point. */
    if (seg1.source() == seg1.target() ||
        seg2.source() == seg2.target())
      return false;

    bool lines_are_equal = m_rat_kernel->equal_3_object()
      (seg1.supporting_line(),seg2.supporting_line());

    // TODO:
    // change this such that it just used the xyz order of the points
    // The segments overlap if the intervals induced by the ordered points
    // overlap

    if (lines_are_equal) {
      if (seg2.has_on(seg1.source()) && seg2.source() != seg1.source() &&
          seg2.target() != seg1.source())
        return true;
      if (seg2.has_on(seg1.target()) && seg2.source() != seg1.source() &&
          seg2.target() != seg1.source())
        return true;

      if (seg1.has_on(seg2.source()) && seg2.source() != seg1.source() &&
          seg2.target() != seg1.source())
        return true;
      if (seg1.has_on(seg2.target()) && seg2.source() != seg1.source() &&
          seg2.target() != seg1.source())
        return true;

      if (seg1.source() == seg2.source() && seg1.target() == seg2.target())
        return true;
      if (seg1.source() == seg2.target() && seg1.target() == seg2.source())
        return true;
    }

    return false;
  }

  void add_element_to_arrangement(const Rational_segment_3& S3) {
    using Point_on_line_and_segment_pair =
      Point_and_segment_pair<Rational_point_3,Rational_segment_3>;

    if (m_S1_S2_on_the_same_line) {
      CGAL::Object result =
        m_rat_kernel->intersect_3_object()(S3.supporting_line(), m_S1_S2_line);
      Rational_point_3 qpoint;
      Rational_line_3 qline;
      if (CGAL::assign(qpoint, result)) {
        if (S3.has_on(qpoint)) {
          if (m_S1_S2_num_of_common_lines == 0)
            m_S1_S2_on_the_same_line_S3 = &S3;
          else if (m_S1_S2_num_of_common_lines == 1)
            m_S1_S2_on_the_same_line_S4 = &S3;

          m_S1_S2_num_of_common_lines += 1;
        }
      }
      else if (CGAL::assign(qline, result)) {
        if (m_S1_S2_num_of_common_lines == 0)
          m_S1_S2_on_the_same_line_S3 = &S3;
        else if (m_S1_S2_num_of_common_lines == 1)
          m_S1_S2_on_the_same_line_S4 = &S3;

        m_S1_S2_num_of_common_lines += 1;
      }

      if (m_S1_S2_intersect)
        add_element_to_sphere_arrangement(S3);
    }
    else if (m_one_dimensional_arr)
      add_element_to_one_dimensional_arrangement(S3);
    else if (are_overlap(*m_S1,S3)) {
      Rational S1_t_min;
      m_g_func.get_scalar_from_point_on_line(S3.source(),
                                             m_S1->supporting_line(),
                                             S1_t_min);

      Rational S1_t_max;
      m_g_func.get_scalar_from_point_on_line(S3.target(),
                                             m_S1->supporting_line(),
                                             S1_t_max);

      if (S1_t_min > S1_t_max) std::swap(S1_t_max,S1_t_min);
      if (S1_t_min < Rational(0)) S1_t_min = Rational(0);
      if (S1_t_max > Rational(1)) S1_t_max = Rational(1);

      typename std::list<Rational_arc_2> arcs_to_insert;
      Arrangement_on_plane_2 *temp_arr_on_plane =
        new Arrangement_on_plane_2();
      Lines_through_segments_arr_observer_on_plane *temp_obs_on_plane =
        new Lines_through_segments_arr_observer_on_plane(*temp_arr_on_plane);

      Rational_arc_2 t_arc;

      m_traits_2_adapt.
        create_segment_on_plane_arr(t_arc,
                                    Rational_point_2(S1_t_min,Rational(0)),
                                    Rational_point_2(S1_t_max,Rational(0)),
                                    &S3);
      arcs_to_insert.push_back(t_arc);

      m_traits_2_adapt.
        create_segment_on_plane_arr(t_arc,
                                    Rational_point_2(S1_t_min,Rational(1)),
                                    Rational_point_2(S1_t_max,Rational(1)),
                                    &S3);
      arcs_to_insert.push_back(t_arc);

      m_traits_2_adapt.
        create_segment_on_plane_arr(t_arc,
                                    Rational_point_2(S1_t_min,Rational(0)),
                                    Rational_point_2(S1_t_min,Rational(1)),
                                    &S3);
      arcs_to_insert.push_back(t_arc);

      m_traits_2_adapt.
        create_segment_on_plane_arr(t_arc,
                                    Rational_point_2(S1_t_max,Rational(0)),
                                    Rational_point_2(S1_t_max,Rational(1)),
                                    &S3);
      arcs_to_insert.push_back(t_arc);

      temp_obs_on_plane->set_is_plane(true);
      temp_obs_on_plane->set_last_inserted_segment(&S3);

      insert (*temp_arr_on_plane, arcs_to_insert.begin(),
              arcs_to_insert.end());

      m_arr_plane_faces.add_element(temp_arr_on_plane,temp_obs_on_plane);

    }
    else if (are_overlap(*m_S2,S3)) {
      Rational S2_t_min;
      m_g_func.get_scalar_from_point_on_line(S3.source(),
                                             m_S2->supporting_line(),
                                             S2_t_min);

      Rational S2_t_max;
      m_g_func.get_scalar_from_point_on_line(S3.target(),
                                             m_S2->supporting_line(),
                                             S2_t_max);

      if (S2_t_min > S2_t_max) std::swap(S2_t_max,S2_t_min);
      if (S2_t_min < Rational(0)) S2_t_min = 0;
      if (S2_t_max > 1) S2_t_max = 1;

      typename std::list<Rational_arc_2> arcs_to_insert;
      Arrangement_on_plane_2 *temp_arr_on_plane =
        new Arrangement_on_plane_2();
      Lines_through_segments_arr_observer_on_plane *temp_obs_on_plane =
        new Lines_through_segments_arr_observer_on_plane(*temp_arr_on_plane);

      Rational_arc_2 t_arc;

      m_traits_2_adapt.
        create_segment_on_plane_arr(t_arc,
                                    Rational_point_2(Rational(0),S2_t_min),
                                    Rational_point_2(Rational(0),S2_t_max),
                                    &S3);
      arcs_to_insert.push_back(t_arc);

      m_traits_2_adapt.
        create_segment_on_plane_arr(t_arc,
                                    Rational_point_2(Rational(1),S2_t_min),
                                    Rational_point_2(Rational(1),S2_t_max),
                                    &S3);
      arcs_to_insert.push_back(t_arc);

      m_traits_2_adapt.
        create_segment_on_plane_arr(t_arc,
                                    Rational_point_2(Rational(0),S2_t_min),
                                    Rational_point_2(Rational(1),S2_t_min),
                                    &S3);
      arcs_to_insert.push_back(t_arc);

      m_traits_2_adapt.
        create_segment_on_plane_arr(t_arc,
                                    Rational_point_2(Rational(0),S2_t_max),
                                    Rational_point_2(Rational(1),S2_t_max),
                                    &S3);
      arcs_to_insert.push_back(t_arc);

      temp_obs_on_plane->set_is_plane(true);
      temp_obs_on_plane->set_last_inserted_segment(&S3);
      insert (*temp_arr_on_plane, arcs_to_insert.begin(),
              arcs_to_insert.end());

      m_arr_plane_faces.add_element(temp_arr_on_plane,temp_obs_on_plane);

    }
    else if (m_S1_S2_on_the_same_plane)
      add_element_to_arrangement_S1S2_on_the_same_plane(S3);
    else {
      m_obs_on_plane->set_is_plane(false);
      if (S3.source() == S3.target()) {
#if ARR_ON_SUR_DEBUG
        std::cout << "S3 is a point " << std::endl;
#endif
        /* Find the plane that created by S3 and S1.*/
        Rational_plane_3
          PlaneS1CL(m_S1->source(),m_S1->target(),S3.source());

        /* Find the intersection of S2 with the plane PlaneS1CL. */
        CGAL::Object result =
          m_rat_kernel->intersect_3_object()(PlaneS1CL, m_S2->supporting_line());
        Rational_point_3 ipoint_S2;

        if (CGAL::assign(ipoint_S2, result) && m_S2->has_on(ipoint_S2))
        {
          Rational S2_t;
          m_g_func.get_scalar_from_point_on_line(ipoint_S2,
                                                 m_S2->supporting_line(),
                                                 S2_t);
          Rational_line_3 S3_ipoint_S2(S3.source(),
                                       ipoint_S2);

          /* Get the intersection of the line created from S3 and
             ipoint_S2 */
          result = m_rat_kernel->intersect_3_object()(S3_ipoint_S2,
                                                      m_S1->supporting_line());
          Rational_point_3 ipoint_S1;
          if (CGAL::assign(ipoint_S1, result) && m_S1->has_on(ipoint_S1)) {
            Rational S1_t;
            m_g_func.get_scalar_from_point_on_line(ipoint_S1,
                                                   m_S1->supporting_line(),
                                                   S1_t);
            Point_on_plane_2 temp_p;
            m_traits_2_adapt.construct_point(temp_p, S1_t, S2_t),
            m_isolated_points_on_plane.add_element
              (Point_on_plane_and_line_pair(Point_2(temp_p,
                                                    S1_t,S2_t),&S3));
          }
        }
      }
      else if (m_rat_kernel->are_parallel_3_object()(m_S1->supporting_line(),
                                                     S3.supporting_line()))
      {
#if ARR_ON_SUR_DEBUG
        std::cout << "rat_kernel->are_parallel_3_object(S1,S3) = "
                  << m_rat_kernel->are_parallel_3_object()(m_S1->supporting_line(),
                                                           S3.supporting_line())
                  << std::endl;
#endif
        Rational S2_t;
        Rational_point_3 ipoint_S2;


        if (m_arr_g_func.calc_parallel_segments(*m_S1, *m_S2, S3, ipoint_S2,
                                                m_rat_kernel,S2_t,true))
        {
          std::list<Rational_arc_2> arcs;
          m_arr_g_func.add_segs_to_arr_S2_is_a_point(*m_S1,
                                                     S3,
                                                     true /* Bound S1. */,
                                                     ipoint_S2,
                                                     arcs,
                                                     m_rat_kernel,
                                                     m_isolated_points_on_plane,
                                                     S2_t,
                                                     m_arr_on_plane.get(),
                                                     S3);
        }
      }
      else if (m_rat_kernel->are_parallel_3_object()(m_S2->supporting_line(),
                                                     S3.supporting_line()))
      {
#if ARR_ON_SUR_DEBUG
        std::cout << "rat_kernel->are_parallel_3_object(S2,S3) = "
                  << m_rat_kernel->are_parallel_3_object()(m_S2->supporting_line(),
                                                           S3.supporting_line())
                  << std::endl;
#endif
        Rational S1_t;
        Rational_point_3 ipoint_S1;



        if (m_arr_g_func.calc_parallel_segments(*m_S2, *m_S1, S3, ipoint_S1,
                                                m_rat_kernel,S1_t, true))
        {
          std::list<Rational_arc_2> arcs;
          m_arr_g_func.add_segs_to_arr_S1_is_a_point(*m_S2,
                                                     S3,
                                                     true /* Bound S2 */,
                                                     ipoint_S1,
                                                     arcs,
                                                     m_rat_kernel,
                                                     m_isolated_points_on_plane,
                                                     S1_t,
                                                     m_arr_on_plane.get(),
                                                     S3);
        }
      }
      else {
        Arr_object obj(m_S1,m_S2,&S3,m_rat_kernel,true /*Bound S2*/ );
#if ARR_ON_SUR_DEBUG
        std::cout << "*******" << obj << std::endl;
#endif
        add_element_to_plane_arrangement(obj,S3);
      }
    }
  }

  /*************************************************************
   * The following function gets a Line S3:
   * 1. Finds the direction vector between the points of the line and the
   * intersection point of S1 and l2.
   * 2. Normalize the line to length 1.
   * 3. Adds the vector arc between the 2 vector ends point to the
   *    arrangement.
   * Note:
   *     Only arcs on the positive side of the sphere will be added inorder
   *     to prevent duplications.
   *     Some arcs might be splitted to 2 arcs.
   *
   * (x,y,z) - intersection point.
   *
   * (x1-x,y1-y,z1-z) - direction vector.
   * Length = sqrt((x1-x)^2 + (y1-y)^2 + (z1-z)^2).
   * ((x1-x)/Length,(y1-y)/Length,(z1-z)/Length) - normalized direction
   * vector.
   *
   * Input:
   *      S3 - Rational segment.
   *
   * Output:
   *
   *************************************************************/
private:
  //! \todo Efi: It would be better to obtain the location directly
  Location_type location(const Rational& x, const Rational& y,
                         const Rational& z) {
    Rational_direction_3 d(x, y, z);
    auto p = m_traits_on_sphere->construct_point_2_object()(d);
    return p.location();
  }

  //!
  void add_element_to_sphere_arrangement(const Rational_segment_3& S3) {
    m_obs_on_sphere->set_is_plane(false);
    Rational p0_x_diff = (S3.source().x() - m_intersection_point_S1S2.x());
    Rational p0_y_diff = (S3.source().y() - m_intersection_point_S1S2.y());
    Rational p0_z_diff = (S3.source().z() - m_intersection_point_S1S2.z());

    Rational p1_x_diff = (S3.target().x() - m_intersection_point_S1S2.x());
    Rational p1_y_diff = (S3.target().y() - m_intersection_point_S1S2.y());
    Rational p1_z_diff = (S3.target().z() - m_intersection_point_S1S2.z());

    std::list<Curve_2> arcs;

    Rational relation_x, relation_y, relation_z;

    if (S3_is_a_point(p0_x_diff, p0_y_diff, p0_z_diff,
                      p1_x_diff, p1_y_diff, p1_z_diff))
    {
      auto loc = location(p0_x_diff, p0_y_diff, p0_z_diff);
      Adapted_point_on_sphere_2 pos(p0_x_diff, p0_y_diff, p0_z_diff, loc);
#if ARR_ON_SUR_DEBUG
      std::cout << "Point on sphere" << std::endl;
#endif
      m_isolated_points_on_sphere.add_element
        (Point_on_sphere_and_line_pair(pos, &S3));
      return;
    }

    /* Add arcs only on the upper half of the sphere. */
    if (p0_z_diff >= 0 && p1_z_diff >= 0) {
      auto loc1 = location(p0_x_diff, p0_y_diff, p0_z_diff);
      auto loc2 = location(p1_x_diff, p1_y_diff, p1_z_diff);
      Adapted_point_on_sphere_2 pos1(p0_x_diff, p0_y_diff, p0_z_diff, loc1);
      Adapted_point_on_sphere_2 pos2(p1_x_diff, p1_y_diff, p1_z_diff, loc2);
      auto ctr = m_traits_on_sphere->construct_curve_2_object();
      auto cv = ctr(pos1.get_original_point(), pos2.get_original_point());
#if ARR_ON_SUR_DEBUG
      std::cout << "Add arc1 = " << cv < std::endl;
#endif
      arcs.push_back(Curve_2(cv, &S3));
    }
    else if (p0_z_diff <= 0 && p1_z_diff <= 0) {
      auto loc1 = location(-p0_x_diff, -p0_y_diff, -p0_z_diff);
      auto loc2 = location(-p1_x_diff, -p1_y_diff, -p1_z_diff);
      Adapted_point_on_sphere_2 pos1(-p0_x_diff, -p0_y_diff, -p0_z_diff, loc1);
      Adapted_point_on_sphere_2 pos2(-p1_x_diff, -p1_y_diff, -p1_z_diff, loc2);
      auto ctr = m_traits_on_sphere->construct_curve_2_object();
      auto cv = ctr(pos1.get_original_point(), pos2.get_original_point());
#if ARR_ON_SUR_DEBUG
      std::cout << "Add arc2 = " << cv << std::endl;
#endif
      arcs.push_back(Curve_2(cv, &S3));
    }
    else {
      /* find the intersection of the line with the plane x,y,
       * in order to work with the half z-positive sphere.
       */
      Rational_plane_3 x_y_plane(Rational_point_3(0,0,0),
                                 Rational_point_3(1,0,0),
                                 Rational_point_3(0,1,0));

      Rational_point_3 rp1(p0_x_diff, p0_y_diff, p0_z_diff);
      Rational_point_3 rp2(p1_x_diff, p1_y_diff, p1_z_diff);
      CGAL::Object result =
        m_rat_kernel->intersect_3_object()(x_y_plane,
                                           Rational_line_3(rp1, rp2));

      Rational_point_3 z_0_ipoint;
      if (!CGAL::assign(z_0_ipoint, result))
      { CGAL_error_msg("Unepxected error"); }

      if (p0_z_diff < 0) {
        auto loc1 = location(z_0_ipoint.x(), z_0_ipoint.y(), z_0_ipoint.z());
        auto loc2 = location(p1_x_diff, p1_y_diff, p1_z_diff);
        Adapted_point_on_sphere_2 pos1(z_0_ipoint.x(), z_0_ipoint.y(), z_0_ipoint.z(),
                               loc1);
        Adapted_point_on_sphere_2 pos2(p1_x_diff, p1_y_diff, p1_z_diff, loc2);
        auto ctr = m_traits_on_sphere->construct_curve_2_object();
        auto cv1 = ctr(pos1.get_original_point(), pos2.get_original_point());
#if ARR_ON_SUR_DEBUG
        std::cout << "Add arc3 = " << cv1 << std::endl;
#endif
        arcs.push_back(Curve_2(cv1, &S3));

        loc1 = location(-z_0_ipoint.x(), -z_0_ipoint.y(), -z_0_ipoint.z());
        loc2 = location(-p0_x_diff, -p0_y_diff, -p0_z_diff);
        Adapted_point_on_sphere_2 qos1(-z_0_ipoint.x(), -z_0_ipoint.y(),
                               -z_0_ipoint.z(), loc1);
        Adapted_point_on_sphere_2 qos2(-p0_x_diff, -p0_y_diff, -p0_z_diff, loc2);
        auto cv2 = ctr(qos1.get_original_point(), qos2.get_original_point());
#if ARR_ON_SUR_DEBUG
        std::cout << "Add arc4 = " << cv2 << std::endl;
#endif
        arcs.push_back(Curve_2(cv2, &S3));
      }
      else {
        auto loc1 = location(z_0_ipoint.x(), z_0_ipoint.y(), z_0_ipoint.z());
        auto loc2 = location(p0_x_diff, p0_y_diff, p0_z_diff);
        Adapted_point_on_sphere_2 pos1(z_0_ipoint.x(), z_0_ipoint.y(), z_0_ipoint.z(),
                               loc1);
        Adapted_point_on_sphere_2 pos2(p0_x_diff, p0_y_diff, p0_z_diff, loc2);
        auto ctr = m_traits_on_sphere->construct_curve_2_object();
        auto cv1 = ctr(pos1.get_original_point(), pos2.get_original_point());
#if ARR_ON_SUR_DEBUG
        std::cout << "Add arc5 = " << cv1 << std::endl;
#endif
        arcs.push_back(Curve_2(cv1, &S3));

        loc1 = location(-z_0_ipoint.x(), -z_0_ipoint.y(), -z_0_ipoint.z());
        loc2 = location(-p1_x_diff, -p1_y_diff, -p1_z_diff);
        Adapted_point_on_sphere_2 qos1(-z_0_ipoint.x(), -z_0_ipoint.y(),
                               -z_0_ipoint.z(), loc1);
        Adapted_point_on_sphere_2 qos2(-p1_x_diff, -p1_y_diff, -p1_z_diff, loc2);
        auto cv2 = ctr(qos1.get_original_point(), qos2.get_original_point());
#if ARR_ON_SUR_DEBUG
        std::cout << "Add arc6 = " << cv2 << std::endl;
#endif
        const auto* traits = arrangement_on_sphere()->geometry_traits();
        arcs.push_back(Curve_2(cv2, &S3));
      }
    }
    CGAL::insert (*m_arr_on_sphere, arcs.begin(), arcs.end());
  }
  /*************************************************************
   * The following function gets an segment S3, and returns true if S3 is a
   * point on the sphere created by m_S1 and m_S2
   *
   * Input:
   *
   * Output:
   *
   *************************************************************/
  bool S3_is_a_point(Rational p0_x_diff, Rational p0_y_diff, Rational p0_z_diff,
                     Rational p1_x_diff, Rational p1_y_diff, Rational p1_z_diff)
  {
    if (p0_x_diff == p1_x_diff &&
        p0_y_diff == p1_y_diff &&
        p0_z_diff == p1_z_diff)
      return true;

    Rational_line_3 l(Rational_point_3(p0_x_diff, p0_y_diff, p0_z_diff),
                      Rational_point_3(p1_x_diff, p1_y_diff, p1_z_diff));

    if (l.has_on(Rational_point_3(0,0,0))) return true;

    return false;
  }

private:
  /*************************************************************
   * The following function gets an hyperbola H1 and add all of its relevant
   * arcs to the arrangement.
   * If one of the arcs degenerated to a point the arc won't be added to the
   * arrangement and the point
   *
   * Input:
   *      H1 - Hyperbola.
   *
   * Output:
   *
   *************************************************************/
  void add_element_to_plane_arrangement(Arr_object& obj,
                                        const Rational_segment_3& S3)
  {

    /*  y = (a0 + x * a1)/(b1 * x + b0) */

    /*
     * Since we are using segments we only need the parts of the hyperbola
     * in the bounded segment [0,1].
     *We will find the intersection of the hyperbola with the lines S2.t = 1
     * and S2.t = 0, if the intersection points are in this segment we will
     * split the hyperbola into two parts or 3 parts.
     * Since the hyperbola will intersect with this lines maximum 2 times it
     * won't affect the asymptotic run time.
     */

#if ARR_ON_SUR_DEBUG
    std::cout <<"get_all_arcs_in_positive_unit_square" <<std::endl;
#endif


#if ONE_SWEEP
    obj.get_all_arcs_in_positive_unit_square(arcs_cont,
                                             m_isolated_points_on_plane);
#else
    typename std::list<Rational_arc_2>  arcs;
    obj.get_all_arcs_in_positive_unit_square(arcs,
                                             m_isolated_points_on_plane);

//     typename std::list<Rational_arc_2>::iterator it;
//     for (it = arcs.begin(); it != arcs.end(); ++it)
//     {
// //       insert (*m_arr_on_plane, *it);
// //       std::cout << *it << std::endl;
//    }

    insert (*m_arr_on_plane, arcs.begin(), arcs.end());
#endif
    // std::cout << "Arr on plane arrangement size:"
    //   << "   V = " << m_arr_on_plane->number_of_vertices()
    //   << ",  E = " << m_arr_on_plane->number_of_edges()
    //   << ",  F = " << m_arr_on_plane->number_of_faces() << std::endl;

  }

  /*************************************************************
   * The following function gets a point and adds a segment which represents
   * all the lines that passes through the point and S1 and S2.
   *
   * Input:
   *      p - point.
   *      S3 - Segment.
   *
   * Output:
   *
   *************************************************************/
  void add_element_to_plane_arrangement(const Rational_point_3& qpoint,
                                        const Rational_segment_3& S3)
  {
    std::list<Arc_end_points > ret_end_points;
#if ARR_ON_SUR_DEBUG
    std::cout<<"Add point to plane arrangement = " << qpoint << std::endl;
#endif

    /*
     * Since we are using segments we only need the parts of the segment in
     * the bounded segment [0,1].
     */

    m_obs_on_plane->set_is_plane(false);
#if ONE_SWEEP
    m_arr_g_func.
      get_all_lines_through_point_and_2_lines(*m_S1,
                                              *m_S2, S3, true /* bound s1 s2*/, qpoint,
                                              *m_rat_kernel, m_intersection_point_S1S2,
                                              m_isolated_points_on_plane, std::back_inserter(arcs_cont),
                                              false, &ret_end_points,
                                              m_S1_S2_intersect);
#else
    std::list<Rational_arc_2> arcs;
    m_arr_g_func.
      get_all_lines_through_point_and_2_lines(*m_S1,
                                              *m_S2, S3, true /* bound s1 s2*/, qpoint,
                                              *m_rat_kernel, m_intersection_point_S1S2,
                                              m_isolated_points_on_plane, std::back_inserter(arcs),
                                              false, &ret_end_points,
                                              m_S1_S2_intersect);

//     typename std::list<Rational_arc_2>::iterator it;
//     for (it = arcs.begin(); it != arcs.end(); ++it)
//     {
// //       insert (*m_arr_on_plane, *it);
// //       std::cout << *it << std::endl;
//     }
    insert (*m_arr_on_plane, arcs.begin(), arcs.end());
#endif
  }

   /*************************************************************
   * The following functor checks that the vertex was created from two edges
   * of two unique hyperbolas. This is important in the case of degenerate
   * hyperbola that has been created from 2 line segements.
   **************************************************************/
  template<typename Arrangement_2>
  class Created_from_2_unique_lines {
  public:

    using Halfedge_around_vertex_const_circulator =
      typename Arrangement_2::Halfedge_around_vertex_const_circulator;
    bool operator()(Halfedge_around_vertex_const_circulator first,
                    bool output_s3, /* When true S3 is also part of the output. */
                    const typename Arrangement_2::Dcel::Ext_obj** obj3,
                    const typename Arrangement_2::Dcel::Ext_obj** obj4) {
      Halfedge_around_vertex_const_circulator curr;

      /* Its sufficient to look only on the first originating segment
         since if there are two distinct segments the entire edge is
         the output.
      */
      const typename Arrangement_2::Dcel::Ext_obj* temp_obj =
        *(first->segs_begin());
      if (output_s3) *obj3 = temp_obj;
      else temp_obj = *obj3 ;
      curr = first;

      do {
        const typename Arrangement_2::Dcel::Ext_obj* obj =
          *curr->segs_begin();

        if (obj != temp_obj) {
          *obj4 = obj;
          return true;
        }
        curr++;
      } while (curr != first);

      return false;
    }
  };

public:
  /*************************************************************
   * The following function runs over all the vertices at the arrangement,
   * for each vertex of degree > 1, it finds the line that passes through the
   * creator lines and the lines which represents the edges of the hyperbola.
   **************************************************************/
  template <typename OutputIterator>
  void find_all_lines(bool rational_output, OutputIterator out) {
#if ONE_SWEEP
    if (arcs_cont.size() != 0)
    {
       // for (typename std::list<Rational_arc_2>::iterator
       //      it = arcs_cont.begin(); it != arcs_cont.end(); ++it)
       // {
       //    std::cout << *it << std::endl;
       //    std::cout << *(it->data()) << std::endl;
       // }
       // std::cout << "Afer insert" << std::endl;

       // for (typename std::list<Rational_arc_2>::iterator
       //      it = arcs_cont.begin(); it != arcs_cont.end(); ++it)
       // {
       //    std::cout << *it << std::endl;
       //    std::cout << *(it->data()) << std::endl;
       // }

       insert (*m_arr_on_plane, arcs_cont.begin(), arcs_cont.end());
    }
#endif

    if (m_S1_S2_on_the_same_line) {
      if (m_S1_S2_num_of_common_lines >= 2) {
        LTS::insert_transversal(out,
                                m_S1_S2_line,
                                m_S1,m_S2,
                                m_S1_S2_on_the_same_line_S3,
                                m_S1_S2_on_the_same_line_S4, With_segments());
      }

      if (m_S1_S2_intersect) find_all_lines_sphere(out);
    }
    else {
      if (m_S1_S2_intersect) find_all_lines_sphere(out);

      /* Merge all of the arrangements at arr_vector to one arrangement. */
      if (m_arr_plane_faces.size() != 0) {
        m_arr_on_plane.
          reset(m_arr_plane_faces.get_overlaid_arr(m_arr_on_plane.get()));
      }

      Created_from_2_unique_lines<Arrangement_on_plane_2> valid_vertex;

      m_arr_g_func.find_all_lines_plane(m_arr_on_plane,
                                        out,
                                        *m_S1,*m_S2,m_S1_S2_intersect,
                                        m_isolated_points_on_plane,
                                        valid_vertex,
                                        m_intersection_point_S1S2,
                                        rational_output);
    }
#if LINES_DEBUG
    /* Print the arrangement size. */
    std::cout << *this << std::endl;
#endif

  }

private:
  template<typename OutputIterator>
  void find_all_lines_sphere(OutputIterator out) {
    using Vertex_it = typename Arrangement_on_sphere_2::Vertex_iterator;
    Vertex_it vit;
    typename Arrangement_on_sphere_2::Edge_iterator   eit;

    if (m_num_of_intersection >= 2) {
      Through_3 through_3(m_intersection_point_S1S2);
      LTS::insert_transversal(out,
                              through_3,
                              m_S1,m_S2,
                              m_S1_S2_concurrent_S3,
                              m_S1_S2_concurrent_S4, With_segments());
    }
    else {
      for (eit = m_arr_on_sphere->edges_begin();
           eit != m_arr_on_sphere->edges_end();
           ++eit)
      {
        /* Add Edges that created by two or more identical curves. */
        if (!eit->get_added_to_output() &&
            (eit->num_of_segments() >=
             (2 - m_num_of_intersection)))
        {
#if ARR_ON_SUR_DEBUG
          std::cout << change_color(CGAL_RED,"ADD PLANE") << std::endl;
          std::cout << eit->curve() << std::endl;
#endif

          Rational_segment_3
            segment(Rational_point_3(m_intersection_point_S1S2.x() +
                                     eit->curve().source().dx(),
                                     m_intersection_point_S1S2.y() +
                                     eit->curve().source().dy(),
                                     m_intersection_point_S1S2.z() +
                                     eit->curve().source().dz()),
                    Rational_point_3(m_intersection_point_S1S2.x() +
                                     eit->curve().target().dx(),
                                     m_intersection_point_S1S2.y() +
                                     eit->curve().target().dy(),
                                     m_intersection_point_S1S2.z() +
                                     eit->curve().target().dz()));

          if (m_num_of_intersection == 1) {
            Through_3 through_3(segment,m_intersection_point_S1S2);
            LTS::insert_transversal(out,
                                    through_3,
                                    m_S1,m_S2,
                                    m_S1_S2_concurrent_S3,
                                    *eit->segs_begin(), With_segments());
          }
          else {
            using Edge = typename Arrangement_on_sphere_2::Halfedge;
            typename Arrangement_on_sphere_2::Dcel::const_iterator it =
              eit->segs_begin();
            typename Arrangement_on_sphere_2::Dcel::const_iterator next_it =
              eit->segs_begin();
            next_it++;

            Through_3 through_3(segment,m_intersection_point_S1S2);
            LTS::insert_transversal(out, through_3, m_S1,m_S2, *it,
                                    *next_it, With_segments());
          }

          eit->set_added_to_output(true);
          eit->twin()->set_added_to_output(true);
          eit->source()->set_added_to_output(true);
          eit->target()->set_added_to_output(true);
        }
      }

      for (vit = m_arr_on_sphere->vertices_begin();
           vit != m_arr_on_sphere->vertices_end(); ++vit)
      {
#if ARR_ON_SUR_DEBUG
        std::cout << "On Sphere vit->degree() = " << vit->degree()
                  << std::endl;
#endif
        Created_from_2_unique_lines<Arrangement_on_sphere_2>
          created_from_2_unique_lines;
        const Rational_segment_3* S3;
        const Rational_segment_3* S4;

        /* Check only vertices that are intersection of 2 hyperbolas. */
        if (!vit->get_added_to_output() &&
            vit->degree() >= (2 - m_num_of_intersection) &&
            created_from_2_unique_lines(vit->incident_halfedges(),true,&S3,&S4))
        {
#if ARR_ON_SUR_DEBUG
          std::cout << "\n\n*************************************"
                    << "NEW LINE SPHERE*****************\n\n"
                    << std::endl;
#endif
          Adapted_point_on_sphere_2 intersection_point = vit->point();
#if ARR_ON_SUR_DEBUG
          std::cout << "vit->degree() = " <<vit->degree() << std::endl;
#endif

          Rational_line_3 transversal
            (m_intersection_point_S1S2,
             Rational_point_3(intersection_point.dx() +
                              m_intersection_point_S1S2.x(),
                              intersection_point.dy() +
                              m_intersection_point_S1S2.y(),
                              intersection_point.dz() +
                              m_intersection_point_S1S2.z()));

          LTS::insert_transversal(out,
                                  transversal,
                                  m_S1,m_S2,S3,S4, With_segments());

#if CGAL_DEBUG_OUTPUT
          Rational_line_3 common_line;
          Rational x = intersection_point.dx();
          Rational y = intersection_point.dy();
          Rational z = intersection_point.dz();

#if ARR_ON_SUR_DEBUG
          std::cout << "x = " << x << std::endl;
          std::cout << "y = " << y << std::endl;
          std::cout << "z = " << z << std::endl;
#endif

          m_g_func.get_line_from_intersection_point(m_intersection_point_S1S2,
                                                    x,y,z,
                                                    common_line);
#if ARR_ON_SUR_DEBUG
          Rational_point_3 source_p(common_line.point(0));
          Rational_point_3 target_p(common_line.point(1));
          Rational temp_scale =
            ((Rational(-20) - source_p.z())/target_p.z());
          if (source_p.z() != Rational(-20)) {
            source_p = Rational_point_3(common_line.point(0).x() + temp_scale *
                                        (common_line.point(1).x() -
                                         common_line.point(0).x()),
                                        common_line.point(0).y() + temp_scale *
                                        (common_line.point(1).y() -
                                         common_line.point(0).y()),
                                        common_line.point(0).z() + temp_scale *
                                        (common_line.point(1).z() -
                                         common_line.point(0).z()));
          }

          temp_scale =
            ((Rational(68) - common_line.point(0).z())/target_p.z());

          if (target_p.z() != Rational(68))
          {
            target_p = Rational_point_3(common_line.point(0).x() + temp_scale *
                                        (common_line.point(1).x() -
                                         common_line.point(0).x()),
                                        common_line.point(0).y() + temp_scale *
                                        (common_line.point(1).y() -
                                         common_line.point(0).y()),
                                        common_line.point(0).z() + temp_scale *
                                        (common_line.point(1).z() -
                                         common_line.point(0).z()));

          }

          std::cout << source_p.x() << " "  << source_p.y()
                    << " " << source_p.z() << ", "
                    << target_p.x() << " "  << target_p.y() << " "
                    << target_p.z() << ","
                    << std::endl;
#endif

#endif
        }
      }

      if (m_isolated_points_on_sphere.size() > 0) {
        Created_from_2_unique_lines<Arrangement_on_sphere_2> valid_vertex;
        using Naive_pl =
          CGAL::Arr_naive_point_location<Arrangement_on_sphere_2>;
        Naive_pl naive_pl;

        m_arr_g_func.add_isolated_points(m_isolated_points_on_sphere,
                                         m_arr_on_sphere, naive_pl,
                                         valid_vertex, out,
                                         *m_S1, *m_S2,
                                         m_intersection_point_S1S2);
      }
    }
  }

private:

  /*************************************************************************
   * The following function handles the case where S1/S2 is a point.
   * It saves the plane/line created by S1 and S2 and mark the arr type as
   * arr_on_plane if S2 is a line and otherwise as arr_on_line.
   *************************************************************************/
  void line_is_a_point(const Rational_point_3& P1,
                       const Rational_segment_3& S2) {
    if (S2.source() == S2.target()) {
      m_S1_S2_line = Rational_line_3(P1,S2.source());
      m_S1_S2_on_the_same_line = true;
    }
    else if (S2.has_on(P1)) {
      m_S1_S2_intersect = true;
      m_intersection_point_S1S2 = P1;
      m_S1_S2_on_the_same_line = true;
      if (P1 != S2.source())
        m_S1_S2_line = Rational_line_3(P1,S2.source());
      else
        m_S1_S2_line = Rational_line_3(P1,S2.target());
      m_obs_on_sphere =
        new Lines_through_segments_arr_observer_on_sphere(*m_arr_on_sphere);
    }
    else if (S2.supporting_line().has_on(P1)) {
      m_S1_S2_line = Rational_line_3(P1,S2.source());
      m_S1_S2_on_the_same_line = true;
    }
    else {
      m_S1_S2_on_the_same_plane = true;
      m_S1_S2_plane = Rational_plane_3(S2.source(),S2.target(),P1);
      m_one_dimensional_arr = true;
    }
  }

public:
  std::string to_string() {
    std::ostringstream o;
    if (m_S1_S2_intersect) {
      o << "Arr on sphere arrangement size:"
        << "   V = " << m_arr_on_sphere->number_of_vertices()
        << ",  E = " << m_arr_on_sphere->number_of_edges()
        << ",  F = " << m_arr_on_sphere->number_of_faces() << std::endl;
    }
    o << "Arr on plane arrangement size:"
      << "   V = " << m_arr_on_plane->number_of_vertices()
      << ",  E = " << m_arr_on_plane->number_of_edges()
      << ",  F = " << m_arr_on_plane->number_of_faces() << std::endl;

    return o.str();
  }

public:
   boost::shared_ptr<Arrangement_on_plane_2> arrangement_on_plane()
   { return m_arr_on_plane; }

   boost::shared_ptr<Arrangement_on_sphere_2> arrangement_on_sphere()
   { return m_arr_on_sphere; }

#if LTS_DRAW_ARR
  void draw_arr() {
    using Arrangement_draw =
      Arrangement_general_functions<Rational_kernel, Alg_kernel,
                                    Arrangement_on_plane_2,
                                    Traits_2_adapt, Point_2>;
    Arrangement_draw arr_draw;

    arr_draw(*m_arr_on_plane);
  }
#endif

};

} //namespace CGAL

#endif /*LINES_THROUGH_SEGMENTS_IMPL_H*/
