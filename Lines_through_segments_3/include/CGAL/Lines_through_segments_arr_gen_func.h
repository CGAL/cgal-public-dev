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

#ifndef LINES_THROUGH_SEGMENTS_ARR_GEN_FUNCTIONS_H
#define LINES_THROUGH_SEGMENTS_ARR_GEN_FUNCTIONS_H

#include <string>

#include <CGAL/basic.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Lines_through_segments_3/internal.h>
#include <CGAL/Lines_through_segments_bounded_seg.h>
#include <CGAL/Lines_through_segments_bounded_segs_vector.h>
#include <CGAL/Lines_through_segments_traits_2_adapt.h>
#include <CGAL/Lines_through_segments_point_adapt.h>
#include <CGAL/Lines_through_segments_isolated_points.h>
#include <CGAL/Lines_through_segments_general_functions.h>
#include <CGAL/Lines_through_segments_output_obj.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

namespace CGAL {

template <typename Lines_through_segments_traits_3,
          typename With_segments>
class Lines_through_segments_arr_gen_func
{
public:
  typedef Lines_through_segments_traits_3 Traits_3;

private:
  typedef typename Traits_3::Alg_kernel      Alg_kernel;
  typedef typename Traits_3::Rational_kernel Rational_kernel;
  typedef typename Traits_3::Traits_arr_on_plane_2  Traits_arr_on_plane_2;
  typedef typename Traits_3::Traits_arr_on_sphere_2 Traits_arr_on_sphere_2;

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

  typedef typename Traits_arr_on_plane_2::Curve_2         Rational_arc_2;
  typedef CGAL::Polyhedron_3<Rational_kernel>             Rational_polyhedron_3;

  typedef typename Traits_arr_on_plane_2::Point_2 Point_2_on_arr;
  typedef Lines_through_segments_point_adapt_2<
    Traits_3,
    typename Traits_arr_on_plane_2::Point_2,Algebraic> Point_2;

  typedef Lines_through_segments_point_adapt_3<
    Traits_3,
    typename Traits_arr_on_sphere_2::Point_2,
    Rational> Point_on_sphere_2;

  Lines_through_segments_general_functions<
    Traits_3> m_g_func;

  typedef Lines_through_segments_mapped_2<
    Traits_3> Mapped_2;


  typedef Lines_through_segments_rbound_unbound_union<Rational>  LTS_rbound;

  const Alg_kernel* m_alg_kernel;

  template <typename Point,typename Ext_obj>
  class Point_and_two_objs
  {
  public:
    Point point;
    const Ext_obj* obj1;
    const Ext_obj* obj2;

    Point_and_two_objs(){ };

    Point_and_two_objs(Point _point,
                       const Ext_obj* _obj1,
                       const Ext_obj* _obj2)
      : point(_point), obj1(_obj1), obj2(_obj2)
    { }

    Point_and_two_objs(const Point_and_two_objs& p2)
    {
      point = p2.point;
      obj1 = p2.obj1;
      obj2 = p2.obj2;
    }
  };

public:

  Lines_through_segments_arr_gen_func(const Alg_kernel* alg_kernel)
  {
    m_alg_kernel = alg_kernel;
  }

  // XXX this is doomed to break and used from arr_object
  Lines_through_segments_arr_gen_func() : m_alg_kernel(NULL)
  {}

   template <typename Halfedge_around_vertex_const_circulator>
   bool is_degenerate_hyp(Halfedge_around_vertex_const_circulator first,
                          Rational_point_2& output_p)
   {
      Lines_through_segments_traits_on_plane_adapt<
      Traits_3> traits_2_adapt;

      Halfedge_around_vertex_const_circulator curr = first;
      Rational y,x;
      do {
         if (traits_2_adapt.is_horizontal(curr->curve())) /* Horizontal line segment*/
         {
            y = traits_2_adapt.get_y_val_of_horizontal_curve(curr->curve());
            curr++;
            curr++;
            if (!traits_2_adapt.has_x_value_at_y(curr->curve(), y))
               return false;
            x = traits_2_adapt.get_x_val(curr->curve(),y);
            output_p = Rational_point_2(x,y);
            return true;
         }
         else if (traits_2_adapt.is_vertical( curr->curve())) /* Vertical line segment*/
         {
            x = traits_2_adapt.get_x_val_of_vertical_cruve(curr->curve());
            curr++;
            curr++;
            if (!traits_2_adapt.has_y_value_at_x(curr->curve(), x))
               return false;

            y = traits_2_adapt.get_y_val_ratioal(curr->curve(), x);
            output_p = Rational_point_2(x,y);
            return true;
         }


         curr++;
      } while (curr != first);
      return false;
   }


  /*************************************************************
   * The following function Iterates over all of the vertics of a face.
   * The vetices represents the bounded lines of a plane that passes through
   * 4 lines.
   **************************************************************/
  template <typename Arrangement_2,
            typename Ccb_halfedge_circulator,
            typename OutputIterator,
            typename Ext_obj>
  void find_plane_faces(boost::shared_ptr<Arrangement_2> arr,
                        Ccb_halfedge_circulator circ,
                        OutputIterator out,
                        const Rational_segment_3& s1,
                        const Rational_segment_3& s2,
                        const Ext_obj& s3,
                        const Ext_obj& s4)
  {
    Ccb_halfedge_circulator curr = circ;
    typedef typename Traits_arr_on_plane_2::X_monotone_curve_2
       X_monotone_curve_2;

    std::list<X_monotone_curve_2> arcs;
    do {
      arcs.push_back(curr->curve());
      curr->source()->set_added_to_output(true);
      curr->target()->set_added_to_output(true);

      curr->set_added_to_output(true);
      curr->twin()->set_added_to_output(true);

    } while (++curr != circ);

    Mapped_2 OAP(arcs.begin(),arcs.end(),s1,s2);
    LTS::insert_transversal(out,
                            OAP,
                            &s1,&s2,&s3,&s4, With_segments());

  }



  /*************************************************************
   * The following function runs over all the vertices at the arrangement,
   * for each vertex of degree > 1, it finds the line that passes through the
   * creator lines and the lines which represents the edges of the hyperbola.
   **************************************************************/
  template <typename Arr_on_plane,
            typename OutputIterator,
            typename Isolated_points_on_plane,
            typename Vertex_valid_functor>
  void find_all_lines_plane(boost::shared_ptr<Arr_on_plane> arr_on_plane,
                            OutputIterator out,
                            const Rational_segment_3& s1,
                            const Rational_segment_3& s2,
                            bool s1_s2_intersect,
                            Isolated_points_on_plane& isolated_points_on_plane,
                            Vertex_valid_functor& validate_vertex,
                            const Rational_point_3& intersection_point_S1S2,
                            bool rational_output)
  {
    /* Run over all of the faces, if face i represents a plane then:
     * 1. Run over all of the half edges, if half edge that doen't not
     * represent a plane,
     *    It a half edge that splits a face which represnts plane.
     *    Hence, all the points on this half edge are valid lines that
     *    passes through 4 segments.
     *    All the lines together represnt half plane.
     * 2. Run over all of the holls at the face, each holl also represent
     *    a plane that crosses 4 segments.
     */
    typename Arr_on_plane::Face_iterator   fit;
    for (fit = arr_on_plane->faces_begin(); fit !=
           arr_on_plane->faces_end();
         ++fit)
    {
      if (fit->num_of_overlap_plane_faces() > 1)
      {
        typename Arr_on_plane::Dcel::const_iterator it = fit->segs_begin();
        typename Arr_on_plane::Dcel::const_iterator next_it = fit->segs_begin();
        ++next_it;
        find_plane_faces(arr_on_plane,
                         fit->outer_ccb(),
                         out,s1,s2,
                         **it,**next_it);
      }
    }

    typename Arr_on_plane::Edge_iterator   eit;
    for (eit = arr_on_plane->edges_begin();
         eit != arr_on_plane->edges_end();
         ++eit)
    {
       /* Add Edges that created by two or more identical curves. */
      if (!eit->get_added_to_output() &&
          ((eit->twin()->face()->num_of_overlap_plane_faces() == 1 &&
            eit->face()->num_of_overlap_plane_faces() == 1 &&
            *(eit->face()->segs_begin()) !=
            *(eit->twin()->face()->segs_begin())) ||
           /* An edge added inside a face with num_of_overlap_plane_faces = 1. */
           (eit->twin()->face()->num_of_overlap_plane_faces() == 1 &&
            eit->face()->num_of_overlap_plane_faces() == 1 &&
            *(eit->face()->segs_begin()) !=
            *eit->segs_begin()) ||
           eit->num_of_segments() >= 2))
      {


#if ARR_ON_SUR_DEBUG
        std::cout << change_color(CGAL_RED,"ADD PLANE") << std::endl;
        std::cout << eit->curve() << std::endl;
        if (eit->twin()->face()->num_of_overlap_plane_faces() == 1 &&
            eit->face()->num_of_overlap_plane_faces() == 1 &&
            *(eit->face()->segs_begin()) !=
            *eit->segs_begin())
        {
          std::cout << "*eit->face()->segs_begin() = "
                    << *eit->face()->segs_begin() << std::endl;
          std::cout << "*eit->segs_begin()" << *eit->segs_begin() << std::endl;
        }
#endif
        const typename Arr_on_plane::Dcel::Ext_obj *S4;
        if (eit->num_of_segments() >= 2)
        {
          typename Arr_on_plane::Dcel::const_iterator next_it =
             eit->segs_begin();
          next_it++;
          S4 = *next_it;
        }
        else
        {
          S4 = *(eit->face()->segs_begin());
        }

        Mapped_2 output_curve(eit->curve(), s1, s2);
        LTS::insert_transversal(out,
                                output_curve,
                                &s1, &s2,
                                *eit->segs_begin(),
                                S4, With_segments());

        eit->source()->set_added_to_output(true);
        eit->target()->set_added_to_output(true);
        eit->set_added_to_output(true);
        eit->twin()->set_added_to_output(true);
      }
    }

    typedef typename Arr_on_plane::Vertex_iterator Vertex_it;
    Vertex_it   vit;
    for (vit = arr_on_plane->vertices_begin();
         vit != arr_on_plane->vertices_end(); ++vit)
    {
#if ARR_ON_SUR_DEBUG
      std::cout<<"On Plane vit->degree() = "<<vit->degree()<<std::endl;
#endif
      /* Check only vertices that are intersection of 2 hyperbolas. */
      if (vit->degree() >= 2 &&
          !vit->get_added_to_output()) /* Was handled earlier. */
      {
        /* Check that the vertex was created from two edges of two
           unique hyperbolas.
           This is important in the case of degenerate hyperbola
           that has been created from
           2 line segements.

        */
        const typename Arr_on_plane::Dcel::Ext_obj *S3,*S4;
        if (validate_vertex( vit->incident_halfedges(),true,&S3,&S4))
        {
          /* In case the point represent the intersection point
             of S1 and S2. */
          Alg_line_3 common_line;
          if (!s1_s2_intersect || m_g_func.add_line_to_output(s1,s2,s1_s2_intersect,vit->point(),
                                                              common_line))
          {
#if ARR_ON_SUR_DEBUG

            std::cout << "S1 = " << s1 << std::endl;
            std::cout << "S2 = " << s2 << std::endl;

            std::cout
              <<   common_line.point(-10)
              << ", "
              <<   common_line.point(10000)
              << ","
              << std::endl;
#endif
           Rational_point_2 rp;

            bool is_rational = (rational_output &&
                                is_degenerate_hyp(vit->incident_halfedges(),
                                                  rp));
            if (is_rational)
            {
              Rational_line_3 output_line;
              m_g_func.get_line_from_intersection_point(rp.x(), rp.y(), s1, s2,
                                                        output_line);
              LTS::insert_transversal(out,
                                      output_line,
                                      &s1, &s2, S3, S4, With_segments());

            }
            else
            {
              Mapped_2 output_point(vit->point(), s1, s2);
              LTS::insert_transversal(out,
                                      output_point,
                                      &s1, &s2, S3, S4, With_segments());
            }
          }

#if CGAL_DEBUG_OUTPUT
          debug_output<Arr_on_plane,Algebraic>(s1, s2, common_line, vit);
#endif

          vit->set_added_to_output(true);
        }
      }
    }

    if (isolated_points_on_plane.size() > 0)
    {
      typedef CGAL::Arr_naive_point_location<Arr_on_plane> Naive_pl;
      Naive_pl     naive_pl;

      add_isolated_points(isolated_points_on_plane, arr_on_plane,
                          naive_pl, validate_vertex,
                          out, s1, s2,
                          intersection_point_S1S2);
    }
  }

  /*************************************************************
   * The following function adds isolated points to the arrangement.
   * If the isolated point falls on edge or on vertex, it returns the line
   * through the 4 segments.
   *
   * Input:
   *      H1 - Hyperbola that will be added to the arrangement.
   *
   * Output:
   *      arr - The arrangemnet.
   *      isolated_points - A data sturcture that contains all
   * of the isolated points.
   *
   *************************************************************/
  template <typename Isolated_points,
            typename Arrangement_2,
            typename Point_location,
            typename Validate_vertex,
            typename OutputIterator>
  void add_isolated_points(Isolated_points& isolated_points,
                           boost::shared_ptr<Arrangement_2> arr,
                           Point_location& pl,
                           Validate_vertex& validate_vertex,
                           OutputIterator out,
                           const Rational_segment_3& S1,
                           const Rational_segment_3& S2,
                           const Rational_point_3& intersection_point_S1S2)
  {

    typename Isolated_points::iterator it;
    typedef typename
      Isolated_points::IP_point_and_line_pair::PL_Point Point;
    typedef typename Arrangement_2::Dcel::Ext_obj Ext_obj;

    /* Each duplicated point represents a line in 3D, that intersect
       4 lines. */
    typedef Point_and_two_objs<Point,Ext_obj> List_element;

    std::list<List_element> duplicated_points;
    typename std::list<List_element>::iterator dup_it;

    isolated_points.sort_points();
    isolated_points.remove_duplicated_points(duplicated_points);

    /* for each duplicated point add a line */
    for (dup_it = duplicated_points.begin();
         dup_it != duplicated_points.end();
         ++dup_it)
    {
       get_3D_point_from_point_on_arr(
          arr,
          dup_it->point,
          out,
          intersection_point_S1S2, S1, S2 ,
          *dup_it->obj1 ,
          *dup_it->obj2);
    }

    /* Perform point location query on each point, if the point is on
     * an edge than it represents a valid line that goes through 4 segments.
     */

    pl.attach(*arr);

    for (it = isolated_points.begin(); it != isolated_points.end(); ++it)
    {
      add_isolated_point_to_arr(
         arr, ((*it).get_point()),S1, S2,
           *((*it).get_segment()),intersection_point_S1S2,
           pl, validate_vertex, out);
    }
  }

  template <typename OutputIterator,
            typename Arrangement_2,
            typename Ext_obj>
  void get_3D_point_from_point_on_arr(boost::shared_ptr<Arrangement_2> arr,
                                      const Point_on_sphere_2& pt,
                                      OutputIterator out,
                                      const Rational_point_3& intersection_point_S1S2,
                                      const Rational_segment_3& S1,
                                      const Rational_segment_3& S2,
                                      const Ext_obj& S3,
                                      const Ext_obj& S4)
  {
#if CGAL_DEBUG_OUTPUT
    Rational x = pt.dx();
    Rational y = pt.dy();
    Rational z = pt.dz();

    Rational_line_3 common_line;
    m_g_func.get_line_from_intersection_point(intersection_point_S1S2,
                                              x,y,z, common_line);
#endif


    Rational_line_3 output_line(
       intersection_point_S1S2,
       Rational_point_3(intersection_point_S1S2.x() + pt.dx(),
                        intersection_point_S1S2.y() + pt.dy(),
                        intersection_point_S1S2.z() + pt.dz()));

    LTS::insert_transversal(out,
                            output_line,
                            &S1,&S2,&S3,&S4,With_segments());
  }

  template <typename OutputIterator,
            typename Arrangement_2,
            typename Ext_obj>
  void get_3D_point_from_point_on_arr(
     boost::shared_ptr<Arrangement_2> arr,
     const Point_2& pt,
     OutputIterator out,
     const Rational_point_3& intersection_point_S1S2,
     const Rational_segment_3& S1,
     const Rational_segment_3& S2,
     const Ext_obj& S3,
     const Ext_obj& S4)
  {
    typedef typename Arrangement_2::Point_2 Point_L;

    Point_L ptl;
    pt.get_original_point(ptl);

    Mapped_2 output_point(ptl, S1, S2);
    LTS::insert_transversal(out,
                            output_point,
                            &S1, &S2, &S3, &S4, With_segments());

#if CGAL_DEBUG_OUTPUT
    Alg_line_3 common_line;
    Algebraic S1_t;
    Algebraic S2_t;
    Lines_through_segments_get_algebraic_number_adapt<
      Traits_3>   get_algebraic_number_adapt;

    S1_t = get_algebraic_number_adapt(pt.x());
    S2_t = get_algebraic_number_adapt(pt.y());

    int status =
      m_g_func.get_line_from_intersection_point(S1_t, S2_t, S1, S2,
                                                common_line);

    if (!m_g_func.do_intersect_line_segment(common_line, S1,m_alg_kernel))
    {
      CGAL_error_msg("Debug Error The line does not intersect with S1!");
    }
    if (!m_g_func.do_intersect_line_segment(common_line, S2,m_alg_kernel))
    {
      CGAL_error_msg("Debug Error The line does not intersect with S2!");
    }

    if (!m_g_func.do_intersect_line_segment(common_line, S3,m_alg_kernel))
    {
      CGAL_error_msg("Debug Error The line does not intersect with S3!");
    }

    if (!m_g_func.do_intersect_line_segment(common_line, S4,m_alg_kernel))
    {
      CGAL_error_msg("Debug Error The line does not intersect with S4!");
    }
#endif

  }

  template <typename Point,
            typename Ext_obj,
            typename Arrangement_2,
            typename Point_location,
            typename Validate_vertex,
            typename OutputIterator>
  void add_isolated_point_to_arr(boost::shared_ptr<Arrangement_2> arr,
                                 const Point& pt,
                                 const Rational_segment_3& S1,
                                 const Rational_segment_3& S2,
                                 const Ext_obj& S3,
                                 const Rational_point_3& intersection_point_S1S2,
                                 Point_location& pl,
                                 Validate_vertex& created_from_2_unique_lines,
                                 OutputIterator out)
  {
    typedef typename Arrangement_2::Point_2 L_point_2;
    typedef typename Arrangement_2::Vertex_const_handle Vertex;
#if ARR_ON_SUR_DEBUG
    std::cout << change_color(CGAL_BLUE,"add_isolated_point_to_arr  = ",
                              pt, "   ",S3) << std::endl;
#endif

    Vertex    v;
    typename Arrangement_2::Halfedge_const_handle  e;
    typename Arrangement_2::Face_const_handle      f;


    L_point_2 ptl;
    pt.get_original_point(ptl);

    CGAL::Object obj = pl.locate(ptl);
    const Ext_obj* S4 = NULL;

    typename Arrangement_2::Halfedge_around_vertex_const_circulator hec;
    if (CGAL::assign(v,obj))
    {
      hec = v->incident_halfedges();
      if (!v->get_added_to_output())
      {
        const Ext_obj* S3_ptr = &S3;
        if (created_from_2_unique_lines(hec, false, &S3_ptr ,&S4))
        {
#if ARR_ON_SUR_DEBUG
          std::cout << change_color(CGAL_GREEN,
                                    "add_isolated_point_to_arr  = ",
                                    pt , "   ", S3)
                    << std::endl;
#endif

          // typedef typename Arrangement_2::Dcel::Ext_obj Ext_obj;
          get_3D_point_from_point_on_arr(
             arr,
             pt,
             out,intersection_point_S1S2,
             S1, S2, S3, *S4);
        }
      }
    }
    else if ( CGAL::assign(e,obj))
    {
#if ARR_ON_SUR_DEBUG
      std::cout << change_color(CGAL_GREEN,
                                "add_isolated_point_to_arr  = ", pt,
                                "   ", S3)
                << std::endl;
#endif
      typename Arrangement_2::Dcel::const_iterator it;
      for (it = e->segs_begin(); it != e->segs_end(); ++it)
      {
        if (&S3 != *it)
        {
          // typedef typename Arrangement_2::Dcel::Ext_obj Ext_obj;
          get_3D_point_from_point_on_arr(arr, pt,
                                         out,
                                         intersection_point_S1S2,
                                         S1, S2, S3, **it);
          return;
        }
      }
    }
    else if ( CGAL::assign(f,obj) && f->num_of_overlap_plane_faces() >= 1)
    {
#if ARR_ON_SUR_DEBUG
      std::cout << change_color(CGAL_GREEN,
                                "add_isolated_point_to_arr  = ", pt ,
                                "   ", S3)
                << std::endl;
#endif
      typename Arrangement_2::Dcel::const_iterator it;
      for (it = f->segs_begin(); it != f->segs_end(); ++it)
      {
        if (&S3 != *it)
        {
          // typedef typename Arrangement_2::Dcel::Ext_obj Ext_obj;
          get_3D_point_from_point_on_arr(arr, pt,
                                         out,
                                         intersection_point_S1S2,
                                         S1, S2, S3, **it);
        }
      }
    }
  }

  template <typename Arrangement_2,typename Algebraic>
  void debug_output(const Rational_segment_3& s1,
                    const Rational_segment_3& s2,
                    const Alg_line_3& common_line,
                    const typename Arrangement_2::Vertex_handle& vit)
  {
    if (!m_g_func.do_intersect_line_segment(common_line,s1,
                                            m_alg_kernel))
    {
      CGAL_error_msg( "Debug Error The line does not intersect with S1!");
    }
    if (!m_g_func.do_intersect_line_segment(common_line, s2,
                                            m_alg_kernel))
    {
      CGAL_error_msg( "Debug Error The line does not intersect with S2!");
    }

    typename Arrangement_2::Halfedge_around_vertex_circulator first,
      curr, next;

    first = curr = vit->incident_halfedges();

    int num_of_segs = 3;

    do {
      CGAL_assertion(curr->num_of_segments() == 1);
      const typename Arrangement_2::Dcel::Ext_obj* obj = *curr->segs_begin();
#if ARR_ON_SUR_DEBUG
      std::cout << "debug output S3 = " << *obj << std::endl;
#endif

      if (!m_g_func.do_intersect_line_segment(common_line, *obj,m_alg_kernel))
      {
#if CGAL_DEBUG_OUTPUT
        std::cout << change_color(CGAL_RED,"OBJ = ",*obj) << std::endl;
#endif
        // CGAL_error_msg(
        //    "Debug Error The line does not intersect with S3!");
        // exit(0);
      }
      if (*curr->segs_begin() != *first->segs_begin())
        num_of_segs++;
      curr++;
    } while (curr != first);

    if ( num_of_segs < 4)
    {
      CGAL_error_msg("Debug Error The line does not intersect with 4 segs!");
    }
  }
  /*************************************************************
   * The following function gets 2 segments and a point on the plane and returns
   * the bounds at the lines of all the line that passes through the point and
   * the 2 lines.
   * The return value is represented as vector of Bounded segments on S1.
   *
   *************************************************************/
  void get_bounded_segments_of_point_and_2_lines(
     const Rational_kernel* rat_kernel,
     const Rational_segment_3& S1,
     const Rational_segment_3& S2,
     bool bound_s1,
     const Rational_point_3& qpoint,
     Lines_through_segments_bounded_segs_vector<Rational>&
     ret_bounded_segments)
  {
    typedef Lines_through_segments_bounded_seg<Rational> Bounded_seg;
    Rational_point_3 ipoint_temp;

    /* Get the intersection points of the 4 lines that passes through the
       segments end points with qpoint. */
    std::list<Rational> end_points_sorted_list;

    Rational_line_3 temp_line;
    CGAL::Object result;

    if (bound_s1)
    {
      /* S1(0) */
      temp_line = Rational_line_3(qpoint,S1.source());
      result =
        rat_kernel->intersect_3_object()(temp_line,S2.supporting_line());

      if (CGAL::assign(ipoint_temp, result))
      {
        if (S2.has_on(ipoint_temp))
        {
          end_points_sorted_list.push_back(Rational(0));
        }
      }

      /* S1(1) */
      temp_line = Rational_line_3(qpoint,S1.target());
      result = rat_kernel->intersect_3_object()(temp_line,S2.supporting_line());
      if (CGAL::assign(ipoint_temp, result))
      {
        if (S2.has_on(ipoint_temp))
        {
          end_points_sorted_list.push_back(Rational(1));

        }
      }
    }

    if (S2.has_on(qpoint))
    {
      if (bound_s1)
      {
        ret_bounded_segments.add_bounded_seg(Bounded_seg(
                                                Rational(1),Rational(0),
                                                true,true));
      }
      else
      {
        ret_bounded_segments.add_bounded_seg(
           Bounded_seg(
              LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY,
              LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY,false,false));
      }
      return;
    }


    /* S2(0) */
    temp_line = Rational_line_3(qpoint,S2.source());
    result = rat_kernel->intersect_3_object()(temp_line,
                                              S1.supporting_line());
    if (CGAL::assign(ipoint_temp, result))
    {
      if (!bound_s1 || S1.has_on(ipoint_temp))
      {
        Rational S1_t;
        m_g_func.get_scalar_from_point_on_line(ipoint_temp,
                                               S1.supporting_line(),S1_t);
        end_points_sorted_list.push_back(S1_t);
      }
    }

    /* S2(1) */
    temp_line = Rational_line_3(qpoint,S2.target());
    result = rat_kernel->intersect_3_object()(temp_line,S1.supporting_line());
    if (CGAL::assign(ipoint_temp, result))
    {
      if (!bound_s1 || S1.has_on(ipoint_temp))
      {
        Rational S1_t;
        m_g_func.get_scalar_from_point_on_line(ipoint_temp,S1.supporting_line(),
                                               S1_t);
        end_points_sorted_list.push_back(S1_t);

      }
    }

    /* Sort the points. */
    end_points_sorted_list.sort();

    /*  Remove duplicated points. */
    typename std::list<Rational>::iterator it =
      end_points_sorted_list.begin();
    typename std::list<Rational>::iterator next_it;

    while (it != end_points_sorted_list.end())
    {
      next_it = it;
      next_it++;
      if (next_it != end_points_sorted_list.end() &&
          *it == *next_it)
      {
        end_points_sorted_list.erase(next_it);
      }
      else
      {
        it++;
      }
    }


    for (it = end_points_sorted_list.begin();
         it != end_points_sorted_list.end();++it)
    {
      next_it = it;
      next_it++;
      if (next_it != end_points_sorted_list.end())
      {
        Rational_point_3 temp_point(S1.source().x() +
                                    (*it + (*next_it - *it)/2) *
                                    (S1.target().x() - S1.source().x()),
                                    S1.source().y() +
                                    (*it + (*next_it - *it)/2) *
                                    (S1.target().y() - S1.source().y()),
                                    S1.source().z() +
                                    (*it + (*next_it - *it)/2) *
                                    (S1.target().z() - S1.source().z()));

        Rational_line_3 temp_line(qpoint,temp_point);

        /* Segment */
        if (m_g_func.do_intersect_line_segment(temp_line,S2,rat_kernel))
        {
          ret_bounded_segments.add_bounded_seg(Bounded_seg(*next_it,*it,
                                                           true,true));
          it++;
        }
        else /* Point */
        {
          if (bound_s1)
          {
            ret_bounded_segments.add_bounded_seg(Bounded_seg(*it,*it,true,true));
          }
          else
          {
            /* If S1 is not bounded the bounded segments is split
               into two parts. The split point is when the line that
               passes through S3 is parallel to S1.*/
            ret_bounded_segments.add_bounded_seg(
               Bounded_seg(
                  LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY,
                  *next_it,false,true));

            ret_bounded_segments.add_bounded_seg(
               Bounded_seg(*it,LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY,
                           true,false));
            it++;
          }
        }
      }
      else /* Point */
      {
        if (bound_s1)
        {
          ret_bounded_segments.add_bounded_seg(Bounded_seg(*it,*it,true,true));
        }
        else
        {
          /* If S1 is not bounded the bounded segments is either
             a point, or a half segment. */
          Rational_point_3 temp_point_plus(S1.source().x() +
                                           (*it + (Rational(1)/Rational(2))) *
                                           (S1.target().x() - S1.source().x()),
                                           S1.source().y() +
                                           (*it + (Rational(1)/Rational(2))) *
                                           (S1.target().y() - S1.source().y()),
                                           S1.source().z() +
                                           (*it + (Rational(1)/Rational(2))) *
                                           (S1.target().z() - S1.source().z()));

          Rational_line_3 temp_line_plus(qpoint,temp_point_plus);

          if (m_g_func.do_intersect_line_segment(temp_line_plus,S2,
                                                 rat_kernel))
          {
            ret_bounded_segments.add_bounded_seg(
               Bounded_seg(LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY,*it,
                           false, true));
          }
          else
          {
            Rational_point_3 temp_point_minus(
               S1.source().x() +
               (*it - (Rational(1)/Rational(2))) *
               (S1.target().x() - S1.source().x()),
               S1.source().y() +
               (*it - (Rational(1)/Rational(2))) *
               (S1.target().y() - S1.source().y()),
               S1.source().z() +
               (*it - (Rational(1)/Rational(2))) *
               (S1.target().z() - S1.source().z()));

            Rational_line_3 temp_line_minus(qpoint,temp_point_minus);
            if (m_g_func.do_intersect_line_segment(temp_line_minus,
                                                   S2,rat_kernel))
            {
              ret_bounded_segments.add_bounded_seg(
                 Bounded_seg(*it,LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY,
                             true, false));
            }
            else
            {
              ret_bounded_segments.add_bounded_seg(Bounded_seg(*it,*it,
                                                               true,true));
            }
          }
        }
      }
    }
  }


  template <typename Isolated_points_on_plane,
            typename Arr_on_plane,
            typename Ext_obj>
  void add_segs_to_arr_S2_is_a_point(
     const Rational_segment_3& S1,
     const Rational_segment_3& S3,
     bool bound_s1,
     const Rational_point_3& P2,
     std::list<Rational_arc_2> & arcs,
     const Rational_kernel* rat_kernel,
     Isolated_points_on_plane& isolated_points_on_plane,
     const Rational& S2_t,
     Arr_on_plane* arr_on_plane,
     const Ext_obj &ext_point_obj)
  {
    Lines_through_segments_traits_on_plane_adapt<
      Traits_3> traits_2_adapt;

    typedef typename Isolated_points_on_plane::IP_point_and_line_pair
      Point_on_plane_and_line_pair;

    Lines_through_segments_bounded_segs_vector<Rational>
      ret_bounded_segments;

    this->get_bounded_segments_of_point_and_2_lines(rat_kernel,
                                                    S1,
                                                    S3,
                                                    bound_s1,
                                                    P2,
                                                    ret_bounded_segments);

    typename Lines_through_segments_bounded_segs_vector<Rational>::iterator
      seg_it;

#if LINES_DEBUG
    std::cout <<"ret bounded segments = "
              << ret_bounded_segments << std::endl;
#endif

    for (seg_it = ret_bounded_segments.begin();
         seg_it != ret_bounded_segments.end();
         ++seg_it)
    {
      /* Add point to the arrangement */
      if ((*seg_it).get_min() == (*seg_it).get_max())
      {
        typedef typename Point_on_plane_and_line_pair::PL_Point Point_2;
        typename Traits_arr_on_plane_2::Point_2 temp_p;
        traits_2_adapt.construct_point(temp_p, ((*seg_it).get_min().bound()),S2_t);
        isolated_points_on_plane.add_element(
           Point_on_plane_and_line_pair(
              Point_2(
                 temp_p,
                 ((*seg_it).get_min().bound()),S2_t),
              &ext_point_obj));
      }
      else
      {
        Rational_arc_2 t_arc;

        if (((*seg_it).get_min() ==
             LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY) &&
            ((*seg_it).get_max() ==
             LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY))
        {
          traits_2_adapt.create_horizontal_curve_on_plane_arr(t_arc,
                                                              S2_t,
                                                              &S3);
        }
        else if ((*seg_it).get_min() ==
                 LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY)
        {
          traits_2_adapt.create_horizontal_curve_on_plane_arr(
             t_arc,
             S2_t,
             ((*seg_it).get_max().bound()),
             false, &S3);
        }
        else if ((*seg_it).get_max() ==
                 LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY)
        {
          traits_2_adapt.create_horizontal_curve_on_plane_arr(
             t_arc,
             S2_t,
             ((*seg_it).get_min().bound()),
             true, &S3);
        }
        else
        {
          traits_2_adapt.create_segment_on_plane_arr(
             t_arc,
             Rational_point_2(((*seg_it).get_min().bound()),S2_t),
             Rational_point_2(((*seg_it).get_max().bound()),S2_t),
             &S3);
        }

#if LINES_DEBUG
        std::cout << t_arc << std::endl;
#endif
        arcs.push_back(t_arc);

      }
      insert (*arr_on_plane, arcs.begin(), arcs.end());
    }
  }
  template <typename Isolated_points_on_plane,
            typename Arr_on_plane,
            typename Ext_obj>
  void add_segs_to_arr_S1_is_a_point(
     const Rational_segment_3& S2,
     const Rational_segment_3& S3,
     bool bound_s2,
     const Rational_point_3& P1,
     std::list<Rational_arc_2>&  arcs,
     const Rational_kernel* rat_kernel,
     Isolated_points_on_plane& isolated_points_on_plane,
     const Rational& S1_t,
     Arr_on_plane* arr_on_plane,
     const Ext_obj &ext_point_obj)
  {
    Lines_through_segments_traits_on_plane_adapt<
      Traits_3> traits_2_adapt;

    typedef typename Isolated_points_on_plane::IP_point_and_line_pair
      Point_on_plane_and_line_pair;

    Lines_through_segments_bounded_segs_vector<Rational>
      ret_bounded_segments;

    this->get_bounded_segments_of_point_and_2_lines(rat_kernel,S2, S3, bound_s2,
                                                    P1, ret_bounded_segments);

    typename Lines_through_segments_bounded_segs_vector<Rational>::iterator
      seg_it;

#if LINES_DEBUG
    std::cout <<"ret bounded segments = "
              << ret_bounded_segments << std::endl;
#endif

    for (seg_it = ret_bounded_segments.begin();
         seg_it != ret_bounded_segments.end();
         ++seg_it)
    {
      /* Add point to the arrangement */
      if (((*seg_it).get_min()) == ((*seg_it).get_max()))
      {
        typedef typename Point_on_plane_and_line_pair::PL_Point Point_2;
        typename Traits_arr_on_plane_2::Point_2 temp_p;
        traits_2_adapt.construct_point(temp_p, S1_t, (*seg_it).get_min().bound());

        isolated_points_on_plane.add_element(
           Point_on_plane_and_line_pair(Point_2(
                                           temp_p,
                                           S1_t, (*seg_it).get_min().bound()),
                                        &ext_point_obj));
      }
      else
      {
        Rational_arc_2 arc;
        if ((*seg_it).get_min() == LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY &&
            (*seg_it).get_max() == LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY)
        {
          traits_2_adapt.create_vertical_segment_on_plane_arr(
             arc,
             Rational_point_2(S1_t,Rational(0)));
        }
        else if ((*seg_it).get_min() ==
                 LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY)
        {
          traits_2_adapt.create_vertical_segment_on_plane_arr(
             arc,
             Rational_point_2(S1_t, ((*seg_it).get_max().bound())),
             false /* Directed down */);

        }
        else if ((*seg_it).get_max() ==
                 LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY)
        {
          traits_2_adapt.create_vertical_segment_on_plane_arr(
             arc,
             Rational_point_2(S1_t, ((*seg_it).get_min().bound())),
             true /* Directed down */);
        }
        else
        {
          traits_2_adapt.create_segment_on_plane_arr(
             arc,
             Rational_point_2(S1_t, ((*seg_it).get_min().bound())),
             Rational_point_2(S1_t, ((*seg_it).get_max().bound())),
             &S3);
        }

#if LINES_DEBUG
        std::cout << arc << std::endl;
        std::cout << "(" << (*seg_it).get_min()
                  << (*seg_it).get_max() << ")" << std::endl;
#endif
        arcs.push_back(arc);
      }
      insert (*arr_on_plane, arcs.begin(), arcs.end());
    }
  }
  /*************************************************************
   * The following function gets 2 segments and a point on the plane and
   * returns all the lines that passes through the point and the 2 lines.
   * The return value is represented as vector of segement on S1 and S2.
   *
   * The segment is represented by a line where the segement is from the first
   * point on the line and the second point on the line
   *************************************************************/
  template <typename Isolated_points_on_plane,
            typename OutputIterator1,
            typename Arc_end_points>
  void get_all_lines_through_point_and_2_lines(
     const Rational_segment_3& S1,
     const Rational_segment_3& S2,
     const Rational_segment_3& S3,
     bool bound_s1_s2,
     const Rational_point_3& qpoint,
     const Rational_kernel& rat_kernel,
     const Rational_point_3& intersection_point_S1S2,
     Isolated_points_on_plane& ret_isolated_points_on_plane,
     OutputIterator1 ret_arcs,
     bool ret_end_points,
     std::list<Arc_end_points >* ret_end_points_list,
     bool S1_S2_intersect)

  {
    typedef typename Isolated_points_on_plane::IP_point_and_line_pair
      Point_on_plane_and_line_pair;

    Lines_through_segments_bounded_segs_vector<Rational>
      ret_bounded_segments;

    Lines_through_segments_traits_on_plane_adapt<
      Traits_3> traits_2_adapt;

    Rational_arc_2 te_arc;

#if LINES_DEBUG
    std::cout <<"qpoint = " << qpoint << std::endl;
#endif
    /* The lines are concurrent, handled at the sphere arrangement */
    if (S1.has_on(qpoint) && S2.has_on(qpoint))
      return;

    if (S2.supporting_line().has_on(qpoint))
    {
      if (S2.has_on(qpoint))
      {
        Rational S2_t;
        m_g_func.get_scalar_from_point_on_line(qpoint,
                                               S2.supporting_line(),
                                               S2_t);

        traits_2_adapt.create_segment_on_plane_arr(
           te_arc,
           Rational_point_2(Rational(0),S2_t),
           Rational_point_2(Rational(1),S2_t),
             &S3);

#if LINES_DEBUG
        std::cout << " push horizontal line " << te_arc << std::endl;
#endif
        *ret_arcs++ = te_arc;

        if (ret_end_points)
        {
          ret_end_points_list->push_back
            (Arc_end_points(Rational_point_2(Rational(0), S2_t),
                            Rational_point_2(Rational(1), S2_t)));
        }
      }
      /* Add isolated point on S2_t = 1 and S1_t = interesection point
       * on S1 and S2.
       * This point represents the line that contains S2 and
       * passes through S1.
       */

      CGAL::Object result =
        rat_kernel.intersect_3_object()(S1.supporting_line(),
                                        S2.supporting_line());
      Rational_point_3 ipoint;
      if (CGAL::assign(ipoint, result) && S1.has_on(ipoint))
      {
        typedef typename Point_on_plane_and_line_pair::PL_Point Point_2;
        Rational S1_t;
        m_g_func.get_scalar_from_point_on_line(ipoint,
                                               S1.supporting_line(),S1_t);
        traits_2_adapt.create_segment_on_plane_arr(
           te_arc,
           Rational_point_2(S1_t,Rational(0)),
           Rational_point_2(S1_t,Rational(1)),
             &S3);

#if LINES_DEBUG
        std::cout << " push vertical line " << te_arc << std::endl;
#endif
        *ret_arcs++ = te_arc;

        if (ret_end_points)
        {
          ret_end_points_list->push_back
            (Arc_end_points(Rational_point_2(S1_t,Rational(0)),
                            Rational_point_2(S1_t,Rational(1))));
        }
      }

      return;
    }

    if (S1.supporting_line().has_on(qpoint))
    {
      if (S1.has_on(qpoint))
      {
        Rational S1_t;
        m_g_func.get_scalar_from_point_on_line(qpoint,
                                               S1.supporting_line(),S1_t);

        traits_2_adapt.create_segment_on_plane_arr(
           te_arc,
           Rational_point_2(S1_t,Rational(0)),
           Rational_point_2(S1_t,Rational(1)),
             &S3);
#if LINES_DEBUG
        std::cout <<" push vertical line " << te_arc << std::endl;
#endif
        *ret_arcs++ = te_arc;

        if (ret_end_points)
        {
          ret_end_points_list->push_back(Arc_end_points
                                         (Rational_point_2(S1_t, Rational(0)),
                                          Rational_point_2(S1_t,Rational(1))));
        }
      }

      /* Add isolated point on S1_t = 1 and S2_t = interesection point on
       * S1 and S2.
       * This isolated point represents the line that contains S1 and
       * passes through S2.
       */
      CGAL::Object result =
        rat_kernel.intersect_3_object()(S1.supporting_line(),
                                        S2.supporting_line());
      Rational_point_3 ipoint;
      if (CGAL::assign(ipoint, result) && S2.has_on(ipoint))
      {
        typedef typename Point_on_plane_and_line_pair::PL_Point Point_2;

        Rational S2_t;
        m_g_func.get_scalar_from_point_on_line(ipoint,
                                               S2.supporting_line(),
                                               S2_t);
        traits_2_adapt.create_segment_on_plane_arr(
           te_arc,
           Rational_point_2(Rational(0),S2_t),
           Rational_point_2(Rational(1),S2_t),
             &S3);

#if LINES_DEBUG
        std::cout << " push horizontal line " << te_arc << std::endl;
#endif
        *ret_arcs++ = te_arc;

        if (ret_end_points)
        {
          ret_end_points_list->push_back
            (Arc_end_points(Rational_point_2(Rational(0), S2_t),
                            Rational_point_2(Rational(1), S2_t)));
        }
      }
      return;
    }


    this->get_bounded_segments_of_point_and_2_lines(&rat_kernel,
                                                    S1,
                                                    S2,
                                                    bound_s1_s2,
                                                    qpoint,
                                                    ret_bounded_segments);

    typename Lines_through_segments_bounded_segs_vector<Rational>::iterator
      seg_it;

#if LINES_DEBUG
    std::cout <<"ret bounded segments = " << ret_bounded_segments
              << std::endl;
#endif

    for (seg_it = ret_bounded_segments.begin();
         seg_it != ret_bounded_segments.end();
         ++seg_it)
    {
      /* Add point to the arrangement */
      if (((*seg_it).get_min().bound()) == ((*seg_it).get_max().bound()))
      {
        Rational y_coord;
        Rational_point_3 point_on_s1(
           S1.source().x() + ((*seg_it).get_max().bound()) *
           (S1.target().x() - S1.source().x()),
           S1.source().y() + ((*seg_it).get_max().bound()) *
           (S1.target().y() - S1.source().y()),
           S1.source().z() + ((*seg_it).get_max().bound()) *
           (S1.target().z() - S1.source().z()));

        /* The intersection point of S1 and S2 is handled at the
         *  arrangement on sphere case.
         */
        if (!S1_S2_intersect ||
            point_on_s1 != intersection_point_S1S2)
        {
          Rational_line_3 line_qpoint_to_s1(point_on_s1,qpoint);
          CGAL::Object result =
            rat_kernel.intersect_3_object()(S2.supporting_line(),
                                            line_qpoint_to_s1);

          Rational_point_3 ipoint_temp;
          if (CGAL::assign(ipoint_temp, result))
          {
            typedef typename
              Point_on_plane_and_line_pair::PL_Point Point_2;
            m_g_func.get_scalar_from_point_on_line(ipoint_temp,
                                                   S2.supporting_line(),
                                                   y_coord);
            typename Traits_arr_on_plane_2::Point_2 temp_p;
            traits_2_adapt.construct_point(temp_p, ((*seg_it).get_min().bound()),
                                             y_coord);

            ret_isolated_points_on_plane.add_element
               (Point_on_plane_and_line_pair(Point_2(
                                                temp_p,
                                                ((*seg_it).get_min().bound()),
                                                y_coord), &S3));
          }
          else
          {
            CGAL_error_msg("Unexpected error should have been handled earlier");
          }
        }
      }
      else
      {
        Rational x_coord[5];
        Rational y_coord[5];
        x_coord[0] = ((*seg_it).get_min().bound());
        x_coord[4] = ((*seg_it).get_max().bound());
        x_coord[1] =
          (x_coord[4]-x_coord[0])*0.25 + ((*seg_it).get_min().bound());
        x_coord[2] =
          (x_coord[4]-x_coord[0])*0.5 + ((*seg_it).get_min().bound());
        x_coord[3] =
          (x_coord[4]-x_coord[0])*0.75 + ((*seg_it).get_min().bound());

        for (int ii = 0; ii < 5; ++ii)
        {
          Rational_point_3 temp_p(S1.source().x() +  x_coord[ii] *
                                  (S1.target().x() - S1.source().x()),
                                  S1.source().y() +  x_coord[ii] *
                                  (S1.target().y() - S1.source().y()),
                                  S1.source().z() +  x_coord[ii] *
                                  (S1.target().z() - S1.source().z()));

          Rational_line_3 temp_line(temp_p,qpoint);

          CGAL::Object result =
            rat_kernel.intersect_3_object()(S2.supporting_line(),temp_line);

          Rational_point_3 ipoint_temp;
          if (CGAL::assign(ipoint_temp, result))
          {
            m_g_func.get_scalar_from_point_on_line(ipoint_temp,
                                                   S2.supporting_line(),
                                                   y_coord[ii]);
          }
          else
          {
            CGAL_error_msg("Unexpected error");
          }
        }

        /* Add horizontal line to the arrangement. */
        if (y_coord[0] == y_coord[2] &&
            y_coord[0] == y_coord[4])
        {

          traits_2_adapt.create_segment_on_plane_arr
            (te_arc,
             Rational_point_2(x_coord[0],y_coord[0]),
             Rational_point_2(x_coord[4],y_coord[4]),
             &S3);
          *ret_arcs++ = te_arc;

          if (ret_end_points)
          {
            ret_end_points_list->push_back(Arc_end_points
                                           (Rational_point_2(x_coord[0],
                                                             y_coord[0]),
                                            Rational_point_2(x_coord[4],
                                                             y_coord[4])));
          }



#if LINES_DEBUG
          std::cout <<" push Segment " << te_arc << std::endl;
#endif

        }
        /* Add line segment to the arrangement. */
        else if (((x_coord[4]-x_coord[2])/(x_coord[2]-x_coord[0])) ==
                 ((y_coord[4]-y_coord[2])/(y_coord[2]-y_coord[0])))
        {


          traits_2_adapt.create_segment_on_plane_arr
            (te_arc,
             Rational_point_2(x_coord[0], y_coord[0]),
             Rational_point_2(x_coord[4], y_coord[4]),
             &S3);
          *ret_arcs++ = te_arc;

          if (ret_end_points)
          {
            ret_end_points_list->push_back(Arc_end_points
                                           (Rational_point_2(x_coord[0],
                                                             y_coord[0]),
                                            Rational_point_2(x_coord[4],
                                                             y_coord[4])));
          }


#if LINES_DEBUG
          std::cout <<" push Segment " << te_arc << std::endl;
#endif
        }
        /* Add hyperbola to the arrangement. */
        else
        {
          Rational_arc_2 arc;

          traits_2_adapt.
            create_curve_on_plane_arr(arc,
                                      Rational_point_2(x_coord[0],y_coord[0]),
                                      Rational_point_2(x_coord[1],y_coord[1]),
                                      Rational_point_2(x_coord[2],y_coord[2]),
                                      Rational_point_2(x_coord[3],y_coord[3]),
                                      Rational_point_2(x_coord[4],y_coord[4]),
                                      &S3);
          *ret_arcs++ = arc;

          if (ret_end_points)
          {
            ret_end_points_list->push_back(Arc_end_points
                                           (Rational_point_2(x_coord[0],
                                                             y_coord[0]),
                                            Rational_point_2(x_coord[4],
                                                             y_coord[4])));
          }

#if LINES_DEBUG
          std::cout <<" push Arc " << arc << std::endl;
#endif
        }
      }
    }
  }

  /*************************************************************
   * Function description:
   * --------------------
   * The following functions gets 2 parallel lines S1 S3 and computes all
   * lines through S2 and these lines.
   * The function returns the scalar on S2 which represent the intersection
   * point of S2 and the plane of S1 and S3.
   *************************************************************/
  bool calc_parallel_segments(const Rational_segment_3& _S1,
                              const Rational_segment_3& _S2,
                              const Rational_segment_3& _S3,
                              Rational_point_3& ipoint_S2,
                              const Rational_kernel* rat_kernel,
                              Rational& res,
                              bool bound_s2)
  {
    typedef typename Rational_kernel::Plane_3 Rational_plane_3;

    Rational_plane_3 PlaneS1CL(_S1.source(),_S1.target(),_S3.source());

    CGAL::Object result =
      rat_kernel->intersect_3_object()(PlaneS1CL, _S2.supporting_line());
    if (CGAL::assign(ipoint_S2, result))
    {
      if (!bound_s2 || (bound_s2 &&  _S2.has_on(ipoint_S2)))
      {
        m_g_func.get_scalar_from_point_on_line(ipoint_S2,_S2,res);
        return true;
      }
    }

    return false;
  }
};

} //namespace CGAL

#endif //LINES_THROUGH_SEGMENTS_ARR_GEN_FUNCTIONS_H
