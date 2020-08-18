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

#ifndef LINE_THROUGH_POLYTOPES_IMPL_H
#define LINE_THROUGH_POLYTOPES_IMPL_H

#include <CGAL/Lines_through_segments_general_functions.h>
#include <CGAL/Lines_through_segments_arr_gen_func.h>
#include <CGAL/Lines_through_segments_isolated_points.h>
#include <CGAL/Lines_through_segments_point_adapt.h>
#include <CGAL/Lines_through_segments_arr_object.h>
#include <CGAL/Lines_through_segments_arr_ext_dcel.h>
#include <CGAL/Lines_through_segments_3/observer.h>

#include <CGAL/Envelope_diagram_1.h>
#include <CGAL/envelope_2.h>
#include <CGAL/Lines_through_polytopes_3/Lines_through_polytopes_con_component.h>
#include <CGAL/Lines_through_polytopes_3/Lines_through_polytopes_intersection.h>
#include <CGAL/Lines_through_polytopes_3/Lines_through_polytopes_overlay.h>

namespace CGAL {
   
const bool bound_s2_true = true;
const bool bound_s1_s2_false = false;
const bool insert_to_arr_false = false;

template <typename Lines_through_segments_traits_3, typename Insert_iterator>
class Lines_through_polytopes_impl {      
  typedef typename Lines_through_segments_traits_3::Traits_arr_on_plane_2  Traits_arr_on_plane_2;
  typedef typename Lines_through_segments_traits_3::Traits_arr_on_sphere_2 Traits_arr_on_sphere_2;
  typedef typename Lines_through_segments_traits_3::Alg_kernel             Alg_kernel;
  typedef typename Lines_through_segments_traits_3::Rational_kernel        Rational_kernel;
  typedef typename Alg_kernel::FT                                          Algebraic;
  typedef typename Rational_kernel::FT                                     Rational;

  typedef typename Alg_kernel::Point_3                    Alg_point_3;
  typedef typename Alg_kernel::Line_3                     Alg_line_3;
  typedef typename Alg_kernel::Segment_3                  Alg_segment_3;
  typedef typename Alg_kernel::Plane_3                    Alg_plane_3;
  typedef typename Alg_kernel::Point_2                    Alg_point_2;
    
  typedef typename Rational_kernel::Point_3               Rational_point_3;
  typedef typename Rational_kernel::Line_3                Rational_line_3;
  typedef typename Rational_kernel::Segment_3             Rational_segment_3;
  typedef typename Rational_kernel::Plane_3               Rational_plane_3;
  typedef typename Rational_kernel::Point_2               Rational_point_2;

  typedef CGAL::Polyhedron_3<Rational_kernel>             Rational_polyhedron_3;
  typedef typename Rational_polyhedron_3::Halfedge_handle Halfedge_handle;
  typedef typename Rational_polyhedron_3::Edge_const_iterator   Edge_iterator;
  typedef typename Rational_polyhedron_3::Vertex_handle   Vertex_handle;
  typedef typename Rational_polyhedron_3::Facet_iterator  Facet_iterator;
         
  typedef typename Rational_polyhedron_3::Halfedge_around_vertex_const_circulator
  Halfedge_around_vertex_const_circulator;

  typedef typename Rational_polyhedron_3::Halfedge_around_facet_const_circulator
  Halfedge_around_facet_const_circulator;
      
         
  typedef typename Traits_arr_on_plane_2::Curve_2            Rational_arc_2;
  typedef typename Traits_arr_on_plane_2::X_monotone_curve_2 Mon_rat_arc_2;

  typedef Lines_through_segments_arr_ext_dcel<Traits_arr_on_plane_2,
                                              Rational_polyhedron_3>
  Dcel_on_plane;

  typedef CGAL::Arrangement_2<Traits_arr_on_plane_2,
                              Dcel_on_plane>  Arrangement_on_plane_2;
   
  typedef Lines_through_segments_arr_observer<Rational_polyhedron_3,
                                              Arrangement_on_plane_2> 
  Lines_through_segments_arr_observer_on_plane;
      
      
  typedef Lines_through_segments_point_adapt_2<
    Lines_through_segments_traits_3,
    typename Traits_arr_on_plane_2::Point_2,Algebraic> Point_2;

  typedef Point_and_segment_pair<Point_2,Rational_polyhedron_3>
  Point_on_plane_and_segment_pair;
      
  typedef Lines_through_segments_isolated_points<
    Point_on_plane_and_segment_pair, 
    Compare_points_on_plane<Point_2, Rational_polyhedron_3> >
  Isolated_points_on_plane;
      
  typedef Lines_through_segments_arr_object<Lines_through_segments_traits_3,
                                            boost::false_type,
                                            boost::false_type>
  Arr_object;

  typedef Lines_through_segments_traits_on_plane_adapt<Lines_through_segments_traits_3>
  Traits_2_adapt;
      
  typedef CGAL::Envelope_diagram_1<Traits_arr_on_plane_2> Diagram_1;
      
  typedef Lines_through_polytopes_con_component<
    Point_2,Rational_arc_2,Traits_2_adapt> LTP_con_comp;
      
  typedef typename LTP_con_comp::Lines_through_polytopes_component LTP_comp;
      
  typedef Lines_through_polytopes_intersection<Lines_through_segments_traits_3>
  LTP_intersection;
private:
      
  
  /* Member variables. */
  Rational_segment_3 m_s1;
  Rational_segment_3 m_s2;
  Edge_iterator m_e1;
  Edge_iterator m_e2;
  const Rational_kernel *m_rational_kernel;
  Traits_2_adapt m_traits_2_adapt;
      
  Arrangement_on_plane_2 m_arr_on_plane;
  Arrangement_on_plane_2 m_s1_s2_rest_arr; /* Holds the invalid faces
                                              generated by the polytopes P1 and P2.
                                           */
  Lines_through_segments_arr_observer_on_plane m_obs_on_plane;
  Lines_through_segments_arr_observer_on_plane m_obs_s1_s2_rest_arr;
  Lines_through_segments_general_functions<Lines_through_segments_traits_3>
  m_g_func;
      
   Lines_through_segments_arr_gen_func<Lines_through_segments_traits_3,
                                       boost::false_type,
                                       boost::false_type>
  m_arr_g_func;
      
  Insert_iterator* m_insert_itertor;
  Isolated_points_on_plane m_isolated_points_on_plane;
      
  /* True if the arrangment is valid, I.e there exists lines
     that passes through s1 and s2, and do not intersect 
     _p1 and _p2.
  */
  bool m_valid;
      
public:      
  Lines_through_polytopes_impl(Rational_polyhedron_3 &_p1,
                               Edge_iterator _e1,
                               Rational_polyhedron_3 &_p2,
                               Edge_iterator _e2,
                               Insert_iterator *insert_it,
                               const Alg_kernel *alg_kernel,
                               const Rational_kernel *rational_kernel)
  {
    m_valid = false;
    m_arr_g_func = 
       Lines_through_segments_arr_gen_func<Lines_through_segments_traits_3,
        boost::false_type,
        boost::false_type>(alg_kernel);
    m_obs_on_plane.attach(m_arr_on_plane);
    m_obs_s1_s2_rest_arr.attach(m_s1_s2_rest_arr);
    m_insert_itertor = insert_it;
    m_s1 = Rational_segment_3(_e1->vertex()->point(),
                              _e1->opposite()->vertex()->point());

    m_s2 = Rational_segment_3(_e2->vertex()->point(),
                              _e2->opposite()->vertex()->point());
         
    // if (m_s1 != Rational_segment_3(Rational_point_3(1,1,-1),Rational_point_3(1,1,1)))
    //   return;
    m_s1 = Rational_segment_3(Rational_point_3(3,3,-3),Rational_point_3(3,3,3));
    m_s2 = Rational_segment_3(Rational_point_3(5,4,1),Rational_point_3(4,4,3));
//     m_s2 = Rational_segment_3(Rational_point_3(2,0,0),Rational_point_3(2.5,0.5,0.5));

    // #if CGAL_DEBUG_SWAP
    //          std::swap(m_s1,m_s2);
    //          std::swap(_e1,_e2);
    //          std::swap(_p1,_p2);
    // #endif
    //          std::cout << "S1  == " << m_s1 << std::endl;
    //          std::cout << "S2  == " << m_s2 << std::endl;
         
    //          {
    // m_s1 = Rational_segment_3(Rational_point_3(-14, 18, -22),Rational_point_3(30, 16, 21));
    // m_s2 = Rational_segment_3(Rational_point_3(-14, 13, -12),Rational_point_3(20, 13, -13));
                m_valid = true;
                return;
    //          }
         
    m_e1 = _e1;
    m_e2 = _e2;
    m_rational_kernel = rational_kernel;
         
    if (m_rational_kernel->do_intersect_3_object()(m_s1.supporting_line(),
                                                   m_s2.supporting_line()))
    {
      CGAL_error_msg("Not supported m_s1 and m_s2 are intersecting m_s1 = ");
    }
    else if (m_rational_kernel->are_parallel_3_object()(m_s1.supporting_line(),m_s2.supporting_line()))
    {
      CGAL_error_msg("Not supported m_s1 and m_s2 are co-planer");
    }
    else
    {
      /* Add to the arrangement the restricttion on M_S1 and M_S2,
         that derives from the polytopes of m_p1 and m_p2. */
      add_p1_p2_restrictions_to_arr(_p1,_p2);
            
      set_arr_is_valid();
    }
  }

  ~Lines_through_polytopes_impl()
  {
  }
      
  void add_polytope(const Rational_polyhedron_3 &_p3)
  {
    Arrangement_on_plane_2 arr_on_plane;
         
    Lines_through_polytopes_overlay_add_edges_in_box<Arrangement_on_plane_2,
      Lines_through_segments_arr_observer_on_plane> overlay_object;
    Edge_iterator unused;
    Valid_return_false ret_false;
         
    /* Find all the arcs of the envelope inside the bounded box. */
    if (m_s1_s2_rest_arr.number_of_vertices() != 0)
    {
      find_polytope_arcs(m_s1_s2_rest_arr, _p3,m_s1,m_s2, 
                         arr_on_plane,ret_false,overlay_object);
    }
    else
    {
      /* Add the bounded box ([0,0], [1,1]) to the arrangment. */
      Arrangement_on_plane_2 bounded_box_arr;
      Lines_through_segments_arr_observer_on_plane 
        obs_bounded_box_arr(bounded_box_arr);
            
      obs_bounded_box_arr.set_last_inserted_segment(NULL,false);
      add_bounded_box(bounded_box_arr);

      find_polytope_arcs(bounded_box_arr, _p3,m_s1,m_s2, 
                         arr_on_plane,ret_false,overlay_object);

    }
         

    //          std::cout << "arr_on_plane size:" 
    //                    << "   V = " << arr_on_plane.number_of_vertices()
    //                    << ",  E = " << arr_on_plane.number_of_edges()
    //                    << ",  F = " << arr_on_plane.number_of_faces() << std::endl;



    /* Add all of the curve in arr_on_plane to m_arr_on_plane. */
    std::list<Rational_arc_2> temp_list;
    typename Arrangement_on_plane_2::Edge_iterator   eit;
    for (eit = arr_on_plane.edges_begin();
         eit != arr_on_plane.edges_end();
         ++eit)
    {
      Rational_arc_2 temp_arc;
      m_traits_2_adapt.create_curve_on_plane_arr(temp_arc,eit->curve());
      temp_list.push_back(temp_arc);
    }
         
    m_obs_on_plane.set_last_inserted_segment(&_p3,false);
    insert(m_arr_on_plane, temp_list.begin(), temp_list.end());
         
    //          std::cout << "m_arr_on_plane size:" 
    //                    << "   V = " << m_arr_on_plane.number_of_vertices()
    //                    << ",  E = " << m_arr_on_plane.number_of_edges()
    //                    << ",  F = " << m_arr_on_plane.number_of_faces() << std::endl;

  }

  void find_all_lines()
  {
    Created_from_2_unique_lines<Arrangement_on_plane_2> vertex_validator;
         
    /* The segments do not intersect. */
    Rational_point_3 point_unused;
         
    m_arr_g_func.find_all_lines_plane(m_arr_on_plane,
                                      m_insert_itertor,
                                      m_s1,
                                      m_s2,
                                      false,/* s1_s2_intersect */
                                      m_isolated_points_on_plane,
                                      vertex_validator,
                                      point_unused);
  }

#if LTS_DRAW_ARR      
  void draw_arr()
  {
    typedef Arrangement_general_functions<Rational_kernel,Alg_kernel, Arrangement_on_plane_2, Traits_2_adapt, Point_2>
      Arrangement_draw;
    Arrangement_draw arr_draw;
         
    arr_draw(m_arr_on_plane);
  }

  void draw_rest_arr()
  {
    typedef Arrangement_general_functions<Rational_kernel,Alg_kernel, Arrangement_on_plane_2, Traits_2_adapt, Point_2>
      Arrangement_draw;
    Arrangement_draw arr_draw;
         
    arr_draw(m_s1_s2_rest_arr);
  }
#endif
         
  bool is_valid()
  {
    return m_valid;
  }
  
private:
  template<typename Arrangement_2>
  class Created_from_2_unique_lines
  {
  public:
     bool operator()
     (typename Arrangement_2::Halfedge_around_vertex_const_circulator first,
      bool output_s3, /* When true S3 is also part of the output. */
      const typename Arrangement_2::Dcel::Ext_obj** obj3,
      const typename Arrangement_2::Dcel::Ext_obj** obj4)
     {
      typename Arrangement_2::Halfedge_around_vertex_const_circulator 
        curr;
            
      /* Its sufficient to look only on the first originating segment
         since if there are two distinct segments the entire edge is 
         the output.
      */
      const typename Arrangement_2::Dcel::Ext_obj *temp_obj = 
        *(first->segs_begin());
      if (output_s3)
        *obj3 = temp_obj;
      else
        temp_obj = *obj3 ;
      curr = first;
            
      do {
        const typename Arrangement_2::Dcel::Ext_obj* obj = 
          *curr->segs_begin();
               
        if (obj != temp_obj)
        {
          *obj4 = obj;
          return true;
        }
        curr++;
      } while (curr != first);
                  
      return false;
    }

  };

  void add_p1_p2_restrictions_to_arr(Rational_polyhedron_3 &_p1,
                                     Rational_polyhedron_3 &_p2)
  {
#if 0// LTS_POLY_IMPL_DEBUG
    LTP_intersection intersection_obj(m_rational_kernel);           

    /* DEBUG start */

    Rational_line_3 l2(m_s2.source(),m_s2.target());
    Rational_line_3 l1(m_s1.source(),m_s1.target());
            
    for (double s1_t = -5; s1_t <= 5; s1_t+=0.5)
    {
               
      std::cout << "s1_t = " << s1_t << std::endl;
              
      Rational_line_3 l3((l1.point(0) + s1_t * (l1.point(1) - l1.point(0))),
                         (l2.point(0) + 10 * (l2.point(1) - l2.point(0))));

      //                std::cout << change_color(CGAL_YELLOW,
      //                "       ",
      //                CGAL::to_double(l3.point(-1).x())," ",
      //                CGAL::to_double(l3.point(-1).y())," ",
      //                CGAL::to_double(l3.point(-1).z()),", ",
      //                CGAL::to_double(l3.point(1).x())," ",
      //                CGAL::to_double(l3.point(1).y())," ",
      //                CGAL::to_double(l3.point(1).z()),",")  << std::endl;

      switch (intersection_obj(_p1,l3))
      {
       case CGAL::ON_BOUNDED_SIDE:
        std::cout <<change_color(CGAL_GREEN,"CGAL::ON_BOUNDED_SIDE") << std::endl;
        break;
       case CGAL::ON_BOUNDARY:
        std::cout <<change_color(CGAL_BLUE,"CGAL::ON_BOUNDARY") << std::endl;
        break;
       case CGAL::ON_UNBOUNDED_SIDE:
        std::cout <<"CGAL::ON_UNBOUNDED_SIDE" << std::endl;
        break;
      }
               
      std::cout << std::endl << std::endl<< std::endl;
               
      l3 = Rational_line_3((l1.point(0) + s1_t * (l1.point(1) - l1.point(0))),
                           (l2.point(0) + 15 * (l2.point(1) - l2.point(0))));

      switch (intersection_obj(_p1,l3))
      {
       case CGAL::ON_BOUNDED_SIDE:
        std::cout <<change_color(CGAL_GREEN,"CGAL::ON_BOUNDED_SIDE") << std::endl;
        break;
       case CGAL::ON_BOUNDARY:
        std::cout <<change_color(CGAL_BLUE,"CGAL::ON_BOUNDARY") << std::endl;
        break;
       case CGAL::ON_UNBOUNDED_SIDE:
        std::cout <<"CGAL::ON_UNBOUNDED_SIDE" << std::endl;
        break;
      }
      std::cout << std::endl << std::endl<< std::endl;
               
    }

    //             std::cout <<"S1 = " << m_s1 << std::endl;
    /* DEBUG end */
    std::cout <<"S2 = " << m_s2 << std::endl;
#endif            

    // #if !CGAL_DEBUG_SWAP
    Arrangement_on_plane_2 arr_on_plane_s1_rest;

    add_p_restrictions_to_arr(_p1,m_s1,m_s2,m_e1,arr_on_plane_s1_rest);
            


    // #endif            
    /* Add the restriction of S2 . */
    // #if CGAL_DEBUG_SWAP
    Arrangement_on_plane_2 arr_on_plane_s2_rest;
    add_p_restrictions_to_arr(_p2,m_s1,m_s2,m_e2,arr_on_plane_s2_rest);


    // #endif
    /* Add all of the curve in arr_on_plane_s1_rest and in arr_on_plane_s2_rest
       to m_s1_s2_rest_arr. */
            
    std::list<Rational_arc_2> temp_list;
    typename Arrangement_on_plane_2::Edge_iterator   eit;
    for (eit = arr_on_plane_s1_rest.edges_begin();
         eit != arr_on_plane_s1_rest.edges_end();
         ++eit)
    {
      Rational_arc_2 temp_arc;
      m_traits_2_adapt.create_curve_on_plane_arr(temp_arc,eit->curve());
      temp_list.push_back(temp_arc);
    }
    m_obs_s1_s2_rest_arr.set_last_inserted_segment(&_p1,false);
    insert(m_s1_s2_rest_arr, temp_list.begin(), temp_list.end());
    temp_list.clear();
            
    for (eit = arr_on_plane_s2_rest.edges_begin();
         eit != arr_on_plane_s2_rest.edges_end();
         ++eit)
    {
      Rational_arc_2 temp_arc;
      m_traits_2_adapt.create_curve_on_plane_arr(temp_arc,eit->curve());
      temp_list.push_back(temp_arc);
    }
            
    m_obs_s1_s2_rest_arr.set_last_inserted_segment(&_p2,false);
    insert(m_s1_s2_rest_arr, temp_list.begin(), temp_list.end());
  }
  
  void add_p_restrictions_to_arr(Rational_polyhedron_3 &p,
                                 Rational_segment_3 &seg1,
                                 Rational_segment_3 &seg2,
                                 Edge_iterator& seg_edge,
                                 Arrangement_on_plane_2& arr_on_plane)
  {
    Lines_through_polytopes_overlay_add_bounded_faces<Arrangement_on_plane_2,
      Lines_through_segments_arr_observer_on_plane> overlay_object;

    Are_equal_he<Edge_iterator> are_equal_he(seg_edge);

    /* Add the bounded box ([0,0], [1,1]) to the arrangment. */
    Arrangement_on_plane_2 bounded_box_arr;
    Lines_through_segments_arr_observer_on_plane 
      obs_bounded_box_arr(bounded_box_arr);
            
    obs_bounded_box_arr.set_last_inserted_segment(NULL,false);
    add_bounded_box(bounded_box_arr);
                     

    find_polytope_arcs(bounded_box_arr, p, seg1,seg2,arr_on_plane,
                       are_equal_he,overlay_object);
            
    //             typedef Arrangement_general_functions<Rational_kernel,Alg_kernel,Arrangement_on_plane_2,
    //                Traits_2_adapt, Point_2>
    //                Arrangement_draw;
    //             Arrangement_draw arr_draw2;

    //             arr_draw2(arr_on_plane);

  }
         
  /* The following function adds find only the relevent segments
     of the polytope and add them to arr_on_plane.
  */
            
  template <typename Overlay_object,typename Valid_edge>
  void find_polytope_arcs(const Arrangement_on_plane_2& restrictions_arr,
                          const Rational_polyhedron_3 &p,
                          const Rational_segment_3 &seg1,
                          const Rational_segment_3 &seg2,
                          Arrangement_on_plane_2& arr_result,
                          Valid_edge& dont_use_edge,
                          Overlay_object& overlay_object)
  {
    /* Run over all of the faces that contains eiher S1 or
       S1 vertices and find the bounding hyperbolas. */
    typename std::list<Rational_arc_2>  arcs;
    Isolated_points_on_plane isolated_points_on_plane;
            

    /* Run over all of the edges and for each edje compute the 
       projection into 2D arc and add it to arc list. */
    for (Edge_iterator e_p = p.edges_begin(); 
         e_p != p.edges_end();
         ++e_p)
    {
      if (!dont_use_edge(e_p))
      {
        Rational_segment_3 s3(e_p->vertex()->point(),
                              e_p->opposite()->vertex()->point());
#if LTS_POLY_IMPL_DEBUG
        std::cout << "S1  == " << m_s1 << std::endl;
        std::cout << "S2  == " << m_s2 << std::endl;
        std::cout << change_color(CGAL_GREEN,"S3 = ",s3) << std::endl;
#endif
        add_arc_to_arcs_list(arcs, isolated_points_on_plane, seg1, seg2,
                             s3, p, bound_s1_s2_false);
      }
            
    }
            
    /* Since we are looking for faces the isolated points
       can be cleared.
    */
    isolated_points_on_plane.clear();

    if (arcs.size() > 0)
    {
      typename LTP_con_comp::iterator
        con_comp_it;
      LTP_con_comp con_components(arcs.begin(),arcs.end());
#if 1//LTS_POLY_IMPL_DEBUG
      std::cout << con_components << std::endl;
#endif
      for (con_comp_it = con_components.begin();
           con_comp_it != con_components.end();
           ++con_comp_it)
      {
        Diagram_1 min_diag;
        Diagram_1 max_diag;
        Arrangement_on_plane_2 arr_on_plane_polytope;
        Lines_through_segments_arr_observer_on_plane 
          obs_on_plane_polytope(arr_on_plane_polytope);
        obs_on_plane_polytope.set_last_inserted_segment(&p,false);
                  
        if ((*con_comp_it)->begin() != (*con_comp_it)->end())
        {
          /* The envelope procedure does not support vertical asymptote,
             hence, we cut them at the biggest/smallest y, and add a short horizontal
             segment to close the envelope.
          */
#if LTS_POLY_IMPL_DEBUG
          std::cout << change_color(CGAL_BLUE,"eliminate_vertical_asymptotes") << std::endl;
#endif
          con_components.eliminate_asymptotes(**con_comp_it);

#if LTS_POLY_IMPL_DEBUG
          std::cout << change_color(CGAL_BLUE,"compute_envelope") << std::endl;
#endif
          compute_envelope(**con_comp_it,isolated_points_on_plane,min_diag,max_diag);
                     
          typename std::list<Rational_arc_2>  arcs_to_insert;
          set_diagram_edges(min_diag,arcs_to_insert);
          set_diagram_edges(max_diag,arcs_to_insert);
                     
          /* Connect the left most top and bottom vertices. */
          connect_vertex(arcs_to_insert,
                         min_diag.rightmost()->left()->point(),
                         max_diag.rightmost()->left()->point());
                     
          connect_vertex(arcs_to_insert,
                         min_diag.leftmost()->right()->point(),
                         max_diag.leftmost()->right()->point());
          insert(arr_on_plane_polytope,arcs_to_insert.begin(),arcs_to_insert.end());

          insert(arr_on_plane_polytope,
                 con_components.vertical_edges_begin(),
                 con_components.vertical_edges_end());
          {
            typedef Arrangement_general_functions<Rational_kernel,Alg_kernel,Arrangement_on_plane_2,
              Traits_2_adapt, Point_2>
              Arrangement_draw;
            Arrangement_draw arr_draw;
            arr_draw(arr_on_plane_polytope);
          }


          std::cout << "arr_on_plane_polytope size:" 
                    << "   V = " << arr_on_plane_polytope.number_of_vertices()
                    << ",  E = " << arr_on_plane_polytope.number_of_edges()
                    << ",  F = " << arr_on_plane_polytope.number_of_faces() << std::endl;
                     

          //                      typename Arrangement_on_plane_2::Edge_iterator   eit;
          //                      for (eit = arr_on_plane_polytope.edges_begin();
          //                           eit != arr_on_plane_polytope.edges_end();
          //                           ++eit)
          //                      {
          //                         std::cout << eit->curve() << std::endl;
          //                      }
          overlay_object(m_traits_2_adapt,restrictions_arr,
                         arr_on_plane_polytope, arr_result ,p);
        }
      }
    }
  }
         
  template <typename Point_2>
  void connect_vertex(typename std::list<Rational_arc_2>& arcs_to_insert,
                      Point_2& bottom,
                      Point_2& top)
  {
    if (!(bottom == top) && (bottom.x() == top.x()))
    {
      Rational_arc_2 arc;
      m_traits_2_adapt.create_segment_on_plane_arr(arc,bottom,top);
      arcs_to_insert.push_back(arc);
    }
  }
         
  template <typename He>
  class Are_equal_he
  {
  private:
    He m_he;
               
  public:
    Are_equal_he(He he)
    {
      m_he = he;
    }
               
    template <typename He1> 
    bool operator()(He1 e1)
    {
      if (e1->vertex() == m_he->vertex() &&
          e1->opposite()->vertex() == m_he->opposite()->vertex())
      {
        return true;
      }
                  
      if (e1->opposite()->vertex() == m_he->vertex() &&
          e1->vertex() == m_he->opposite()->vertex())
      {
        return true;
      }
                  
      return false;
    }
  };
         
  class Valid_return_false
  {
  public:
    Valid_return_false()
    {
    }
               
    template <typename He1> 
    bool operator()(He1 e1)
    {
      return false;
    }
  };
         
  template <typename He1, typename He2> 
  bool has_common_vertex(He1 e1,He2 e2)
  {
    //             std::cout << "has_common_vertex:" << e1->vertex()->point() <<
    //                " " << e1->opposite()->vertex()->point() << std::endl;

    //             std::cout << "has_common_vertex:" << e2->vertex()->point() <<
    //                " " << e2->opposite()->vertex()->point() << std::endl;
            
    if (e1->vertex() == e2->vertex() ||
        e1->vertex() == e2->opposite()->vertex())
    {
      return true;
    }

    if (e1->opposite()->vertex() == e2->vertex() ||
        e1->opposite()->vertex() == e2->opposite()->vertex())
    {
      return true;
    }

    return false;
            
  }
         
  //          void insert_digram_to_arr(
  //             Arrangement_on_plane_2& arr_on_plane_polytope,
  //             const Diagram_1& diag)
  //          {
  //             typename Diagram_1::Edge_const_handle     e = diag.leftmost();
  //             typename Diagram_1::Vertex_const_handle   v;
  //             typename std::list<Rational_arc_2> arcs;
  //             while (e != diag.rightmost())
  //             {
  //                if (! e->is_empty())
  //                {
  //                   /* Cast the arc to Curve_2 from X monotone curve. */
  //                   Rational_arc_2 arc;
  //                   m_traits_2_adapt.create_curve_on_plane_arr(arc,e->curve());
  //                   arcs.push_back(arc);
  //                }
               
  //                v = e->right();
  //                e = v->right();
  //             }
  //             CGAL_assertion (e->is_empty());
  //             insert(arr_on_plane_polytope,arcs.begin(),arcs.end());
  //             return;
  
  //          }
         
  /*! Print the given envelope diagram. */
  void print_diagram (const Diagram_1& diag)
  {
    typename Diagram_1::Edge_const_handle     e = diag.leftmost();
    typename Diagram_1::Vertex_const_handle   v;
            
    while (e != diag.rightmost())
    {
      std::cout << "Edge: ";
      if (! e->is_empty())
      {
        std::cout << e->curve() << std::endl;
      }
      else
        std::cout << " [empty]" << std::endl;
               
      v = e->right();
      std::cout << "Vertex x = (" << CGAL::to_double(v->point().x()) 
        //                          << ' '
        //                          << CGAL::to_double(m_traits_2_adapt.get_y_val(e->curve(),v->point().x())) 
                << ')' 
                << std::endl;
               
      e = v->right();
    }
    CGAL_assertion (e->is_empty());
    std::cout << "Edge: [empty]" << std::endl;
            
    return;
  }

  void set_diagram_edges (const Diagram_1& diag,
                          std::list<Rational_arc_2>& arcs)
  {
    typename Diagram_1::Edge_const_handle     e = diag.leftmost();
    typename Diagram_1::Vertex_const_handle   v;
    typename Diagram_1::Vertex_const_handle left_ver,right_ver;
            
    Rational_arc_2 prev_arc;/* The right vertex of the last inserted edge. */
    bool first = true;
            
    while (e != diag.rightmost())
    {
      if (! e->is_empty())
      {
        right_ver = e->right();
        left_ver = e->left();

        Rational_arc_2 new_arc;
        m_traits_2_adapt.create_curve_on_plane_arr(new_arc,
                                                   left_ver->point().x(),
                                                   right_ver->point().x(),
                                                   e->curve());
        arcs.push_back(new_arc);

        /* In case the previous right vertex is not equal to
           the left vertex the vertices are on a vertical line. */
        if (!first && !(prev_arc.right() == new_arc.left()))
        {
          /* Add a vertical line. */
          Rational_arc_2 arc;
          m_traits_2_adapt.create_segment_on_plane_arr(arc,
                                                       new_arc.left(),
                                                       prev_arc.right());
          arcs.push_back(arc);
        }
        prev_arc = new_arc;
                  
        first = false;
      }

      v = e->right();
      e = v->right();
    }
    CGAL_assertion (e->is_empty());
            
    return;
  }

  void insert_envelope_arcs_to_plane_arr(const Diagram_1& diag,
                                         typename std::list<Rational_arc_2>& arcs,
                                         const Rational_polyhedron_3& p3)
  {
    typename Diagram_1::Edge_const_handle     e = diag.leftmost();
    typename Diagram_1::Vertex_const_handle   v;
            
            
    while (e != diag.rightmost())
    {
      Point_2 source(e->curve().source());
      Point_2 target(e->curve().target());

      if (! e->is_empty() &&
          !(((source.y() > Rational(1) && target.y() > Rational(1)) ||
             (source.y() < Rational(0) && target.y() < Rational(0)))))
      {
        Rational_arc_2 arc;
        Rational coefficients[4];
        coefficients[3] = e->curve().t();
        coefficients[2] = e->curve().u();
        coefficients[1] = e->curve().v();
        coefficients[0] = e->curve().w();

        /* Add only the intersection of the cureve with the 
           square [0,0],[1,1]*/
        if (e->number_of_curves() > 1)
        {
          std::cout << change_color(CGAL_RED,
                                    "TODO more than one curve") << std::endl;
        }

        if (source.y() <= Rational(1) &&
            source.y() >= Rational(0) &&
            target.y() <= Rational(1) &&
            target.y() >= Rational(0))
        {
          /* The entire curve is inside the square. */
          m_traits_2_adapt.create_curve_on_plane_arr(arc, source, target,
                                                     coefficients);
          arcs.push_back(arc);

        }
        else 
        {
          typedef typename Traits_arr_on_plane_2::Point_2 Arr_point_2;
          Arr_point_2 source_on_arr,target_on_arr;
                     
          if (source.y() > Rational(1))
          {
            Rational numerator = - (e->curve().w() + e->curve().v());
            Rational denominator = (e->curve().t() + e->curve().u());
            Rational x = (numerator / denominator);
            source_on_arr = Arr_point_2(x,Rational(1));
          }
          else if (source.y() < Rational(0))
          {
            Rational numerator = - (e->curve().w());
            Rational denominator = (e->curve().u());

            Rational x = (numerator / denominator);
            (e->curve().w())/ (e->curve().u());
            source_on_arr = Arr_point_2(x,Rational(0));

          }
          else
          {
            source_on_arr = e->curve().source();
          }

          if (target.y() > Rational(1))
          {
            /* Only the right vertex is inside the square. */
            Rational numerator = - (e->curve().w() + e->curve().v());
            Rational denominator = (e->curve().t() + e->curve().u());
            Rational x = (numerator / denominator);
            target_on_arr = Arr_point_2(x,Rational(1));
          }
          else if (target.y() < Rational(0))
          {
            Rational numerator = - (e->curve().w());
            Rational denominator = (e->curve().u());

            Rational x = (numerator / denominator);
            target_on_arr = Arr_point_2(x,Rational(0));
          }
          else
          {
            target_on_arr = e->curve().target();
          }
                     
          if (source_on_arr == target_on_arr)
          {
            m_isolated_points_on_plane.
              add_element(Point_on_plane_and_segment_pair(source_on_arr, &p3));    
          }
          else
          {
            m_traits_2_adapt.create_curve_on_plane_arr(arc, source_on_arr,
                                                       target_on_arr,
                                                       coefficients);
            arcs.push_back(arc);
          }
                     
        }
      }
               
      v = e->right();
      e = v->right();
    }
            
    return;
  }

  void set_arr_is_valid()
  {
    typedef typename Arrangement_on_plane_2::Face_const_handle     
      Face_const_handle;
    typedef typename Arrangement_on_plane_2::Hole_const_iterator     
      Hole_const_iterator;

    /* The arr is not valid if the entire square ([0,0],[1,1])
       is not valid.
       Set to false if the unbounded face contains the vertices
       (0,0), (0,1), (1,0) and (1,1).
       There is only one unbounded face.
    */
    if (m_s1_s2_rest_arr.number_of_faces() >= 2)
    {
      Face_const_handle unbounded_face = m_s1_s2_rest_arr.unbounded_face();
      Hole_const_iterator hit = unbounded_face->holes_begin();
               
      if (hit == unbounded_face->holes_end())
      {
        m_valid = true;
        return;
      }
      typedef typename Arrangement_on_plane_2::Ccb_halfedge_const_circulator
        Ccb_halfedge_const_circulator;

      Ccb_halfedge_const_circulator circ = *hit;
      Ccb_halfedge_const_circulator curr = circ;
               
      typename Arrangement_on_plane_2::Point_2 p1(Rational(0),Rational(0));
      typename Arrangement_on_plane_2::Point_2 p2(Rational(0),Rational(1));
      typename Arrangement_on_plane_2::Point_2 p3(Rational(1),Rational(0));
      typename Arrangement_on_plane_2::Point_2 p4(Rational(1),Rational(1));
      bool p1_found = false;
      bool p2_found = false;
      bool p3_found = false;
      bool p4_found = false;
                
      do {
        if (curr->source()->point() == p1)
          p1_found = true;
        if (curr->source()->point() == p2)
          p2_found = true;
        if (curr->source()->point() == p3)
          p3_found = true;
        if (curr->source()->point() == p4)
          p4_found = true;
      } while (++curr != circ);
               
      if (p1_found && p2_found && p3_found && p4_found)
      {
        m_valid = false;
        return;
      }
               
               
    }
    m_valid = true;
  }
  
  void compute_envelope(const LTP_comp& comp,
                        Isolated_points_on_plane& isolated_points_on_plane,
                        Diagram_1&  min_diag,
                        Diagram_1&  max_diag)
  {
    {
      Arrangement_on_plane_2 temp_arr;
      insert(temp_arr,comp.begin(), comp.end());
      typedef Arrangement_general_functions<Rational_kernel,Alg_kernel,Arrangement_on_plane_2, Traits_2_adapt, Point_2>
        Arrangement_draw;
      Arrangement_draw arr_draw;
      arr_draw(temp_arr);
    }
         
    if (isolated_points_on_plane.size() > 0)
    {
      CGAL_error_msg("TODO");
      std::list<Rational_arc_2> arcs_plus_points;
      typename LTP_comp::const_iterator it_arcs_list;
             
      for (it_arcs_list = comp.begin();
           it_arcs_list != comp.end();
           it_arcs_list++)
      {
        arcs_plus_points.push_back(*it_arcs_list);
      }
             

      typename Isolated_points_on_plane::iterator it;
      for (it = isolated_points_on_plane.begin();
           it != isolated_points_on_plane.end();
           it++)
      {
        if (it->get_point().y() <= Rational(1))
        {
          Rational_arc_2 ver_seg;

          m_traits_2_adapt.create_segment_on_plane_arr
            (ver_seg,
             it->get_point().get_rational_point(),
             Rational_point_2(it->get_point().get_rational_point().x(),
                              Rational(2)));

          arcs_plus_points.push_back(ver_seg);
        }
      }
            
      lower_envelope_2(arcs_plus_points.begin(), arcs_plus_points.end(),
                       min_diag);
    }
    else
    {
      lower_envelope_2 (comp.begin(), comp.end(),
                        min_diag);
    }
    //          std::cout << "LOWER ENVELOPE" << std::endl;
         
    //          print_diagram (min_diag);
         
    /* Compute the maximization diagram that represents the upper envelope. */
         
         
    //          std::cout << "UPPER ENVELOPE" << std::endl;
    if (isolated_points_on_plane.size() > 0)
    {
      std::list<Rational_arc_2> arcs_plus_points;
      typename LTP_comp::const_iterator it_arcs_list;

      for (it_arcs_list = comp.begin();
           it_arcs_list != comp.end(); 
           it_arcs_list++)
      {
        arcs_plus_points.push_back(*it_arcs_list);
      }

      typename Isolated_points_on_plane::iterator it;
      for (it = isolated_points_on_plane.begin();
           it != isolated_points_on_plane.end();
           it++)
      {
        if (it->get_point().y() >= Rational(0))
        {
          Rational_arc_2 ver_seg;
          m_traits_2_adapt.create_segment_on_plane_arr
            (ver_seg,
             it->get_point().get_rational_point(),
             Rational_point_2(it->get_point().get_rational_point().x(),
                              Rational(-1)));

          arcs_plus_points.push_back(ver_seg);
        }
      }
            
      lower_envelope_2(arcs_plus_points.begin(), arcs_plus_points.end(),
                       max_diag);
    }
    else
    {
      upper_envelope_2 (comp.begin(), comp.end(),
                        max_diag);
    }
#if LTS_POLY_IMPL_DEBUG         
    print_diagram (min_diag);
    print_diagram (max_diag);
#endif
  }
         

  void add_arc_to_arcs_list(typename std::list<Rational_arc_2>&  arcs,
                            Isolated_points_on_plane& isolated_points_on_plane,
                            const Rational_segment_3& s1,
                            const Rational_segment_3& s2,
                            const Rational_segment_3& s3,
                            const Rational_polyhedron_3& p3,
                            bool bound_s1_s2)
  {
    /* Only bound S1 and S3, S2 will be bounded on later stage. */

    /* Bound e1 with the adjacnt facets, and than bound the source and
       target with the faces they belong two except of the adjacnt facets
       of e1.
    */
    //         std::cout << "Add the following arc " << s3 << std::endl;
        
    if (m_rational_kernel->are_parallel_3_object()(s1.supporting_line(),
                                                   s3.supporting_line()))
    {
#if LTS_POLY_IMPL_DEBUG
      std::cout << change_color(CGAL_BLUE,"rational_kernel->are_parallel_3_object(S1,S3) = ")
                << std::endl;
#endif
      Rational S2_t;
      Rational_point_3 ipoint_S2;
           
      if (m_arr_g_func.calc_parallel_segments(s1, s2, s3, ipoint_S2,
                                              m_rational_kernel,S2_t,
                                              bound_s1_s2))
      {
        m_arr_g_func.add_segs_to_arr_S2_is_a_point(s1, s3, false /* bound S1*/,
                                                   ipoint_S2, 
                                                   arcs,
                                                   m_rational_kernel,
                                                   isolated_points_on_plane,
                                                   S2_t,
                                                   &m_arr_on_plane,
                                                   insert_to_arr_false, p3);
      }
    }
    else if (m_rational_kernel->are_parallel_3_object()(s2.supporting_line(),
                                                        s3.supporting_line()))
    {
#if LTS_POLY_IMPL_DEBUG
      std::cout << "rational_kernel->are_parallel_3_object(S2,S3) = "
                << std::endl;
#endif
      Rational S1_t;
      Rational_point_3 ipoint_S1;
           
      if (m_arr_g_func.calc_parallel_segments(s2, s1, s3, ipoint_S1,
                                              m_rational_kernel,S1_t,
                                              bound_s1_s2))
      {
        m_arr_g_func.add_segs_to_arr_S1_is_a_point(s2, s3, false /* bound S2*/,
                                                   ipoint_S1, 
                                                   arcs,m_rational_kernel,
                                                   isolated_points_on_plane,
                                                   S1_t,&m_arr_on_plane,
                                                   insert_to_arr_false, p3);
      }
    }
    else 
    {
      Arr_object obj(&s1, &s2, &s3, m_rational_kernel, false);
      obj.get_all_arcs(arcs, isolated_points_on_plane, p3);
#if LTS_POLY_IMPL_DEBUG
      std::cout << obj << std::endl;
#endif
    }
  }

  void add_bounded_box(Arrangement_on_plane_2& bounded_box_arr)
  {
    typename std::list<Rational_arc_2>  arcs_to_insert;
    Rational_arc_2 arc;
        
    /* Add the bounded box ([0,0], [1,1]) to the arrangment. */
    /* Push the 4 lines that bound the square [0,0], [1,1]. */
    m_traits_2_adapt.create_segment_on_plane_arr(arc,
                                                 Rational_point_2(Rational(0),
                                                                  Rational(0)),
                                                 Rational_point_2(Rational(0),
                                                                  Rational(1)));
    arcs_to_insert.push_back(arc);

    m_traits_2_adapt.create_segment_on_plane_arr(arc,
                                                 Rational_point_2(Rational(0),
                                                                  Rational(1)),
                                                 Rational_point_2(Rational(1),
                                                                  Rational(1)));
    arcs_to_insert.push_back(arc);
        
    m_traits_2_adapt.create_segment_on_plane_arr(arc,
                                                 Rational_point_2(Rational(1),
                                                                  Rational(1)),
                                                 Rational_point_2(Rational(1),
                                                                  Rational(0)));
    arcs_to_insert.push_back(arc);
        
    m_traits_2_adapt.create_segment_on_plane_arr(arc,
                                                 Rational_point_2(Rational(1),
                                                                  Rational(0)),
                                                 Rational_point_2(Rational(0),
                                                                  Rational(0)));
    arcs_to_insert.push_back(arc);
        
    insert(bounded_box_arr,arcs_to_insert.begin(),arcs_to_insert.end());
  }
};
   
} //namespace CGAL

#endif /*LINE_THROUGH_POLYTOPES_IMPL_H*/
