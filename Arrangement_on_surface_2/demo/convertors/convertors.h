//Function templates
#ifndef files_convertors_h
#define files_convertors_h
#include <iostream>
#include <vector>
/*<<<<<<< HEAD
=======*/

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_2.h>

#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Combinatorial_map_operations.h>
/*
>>>>>>> fa8bd1acbcc9d226f464224fe8b956d0035fabb2
 */
/*
#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_constructors.h>
#include <CGAL/Combinatorial_map_operations.h>
 */

namespace CGAL {
  /** @file convertors.h
   *
  */
/*
//Function template that transfer Arrangement to Linear_cell_complex
template <class Arrangement, class LCC>
typename LCC::Dart_handle arr2lcc(const Arrangement &arr, LCC& lcc)
{
  lcc.clear();
  typename LCC::Dart_handle th;
  typename Arrangement::Face_const_iterator fit;
  //vectors
  std::vector<typename LCC::Point> po;
  std::vector<typename LCC::Dart_handle> dart;
  std::vector<typename Arrangement::Halfedge_const_handle> hch;
  typename LCC::Dart_handle dh;
  //Iteration on the edge to get the coordinates of points
  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    if (!fit->is_unbounded()) {
      //Get the outer boundary
      typename Arrangement::Ccb_halfedge_const_circulator circ =
      fit->outer_ccb();
      typename Arrangement::Ccb_halfedge_const_circulator curr = circ;
      typename Arrangement::Ccb_halfedge_const_circulator curr2 = circ;
      //Store the last edge in this facet
      --curr2;
      //Store the first point obtained
      typename Arrangement::Traits_2::Point_2 origin;
      int size = hch.size();
      //Dart_handle for sew function
      typename LCC::Dart_handle d0;
      bool first = true;
      do {
        typename Arrangement::Halfedge_const_handle he = curr;
        //push the Dart_handle into the vector
        hch.push_back(he);
        //Consider the conditiono that this edge is the first one obtained
        //from the facet
        if (first) {
          typename Arrangement::Traits_2::Point_2 p1 = he->source()->point();
          origin = p1;
          typename Arrangement::Traits_2::Point_2 p2 = he->target()->point();
          ++curr;
          he = curr;
          hch.push_back(he);
          typename Arrangement::Traits_2::Point_2 p3 = he->target()->point();
          typename LCC::Point v1 =typename LCC::Point(p1.x(), p1.y());
          typename LCC::Point v2 =typename LCC::Point(p2.x(), p2.y());
          typename LCC::Point v3 =typename LCC::Point(p3.x(), p3.y());
          dh = lcc.make_triangle(v1, v2, v3);
          dart.push_back(dh);
          dh = dh->beta(1);
          dart.push_back(dh);
          d0 = dh->beta(1);
          first = false;
        }
        //Consider the condition that this edge is the second last one
        else if (he != curr2) {
          typename Arrangement::Traits_2::Point_2 p1 = he->source()->point();
          typename Arrangement::Traits_2::Point_2 p2 = he->target()->point();
          dh = lcc.make_triangle(origin, p1, p2);
          dh = dh->beta(1);
          dart.push_back(dh);
          lcc.template sew<2>(dh->beta(0), d0);
          // typename LCC::Dart_handle d0 = dh->beta(1)->beta(1);
          CGAL::remove_cell<LCC, 1>(lcc,d0);
          d0 = dh->beta(1);
        }
        //If the edge is the last one obtained from the facet
        else dart.push_back(d0);
      } while (++curr != circ);
      curr = circ;
      //Iteration over the edge to see whether this facet is adjacent to
      //another facet
      int pos = 0;
      do {
        typename Arrangement::Halfedge_const_handle he = curr;
        for (int i = 0; i < size; i++) {
          d0 = d0->beta(1);
          if (hch[i] == he->twin()) {
            lcc.template sew<2>(dart[i], dart[size+pos]);
          }
        }
        ++pos;
      } while (++curr != circ);
    }
  }
  CGAL_assertion(dh!=NULL);
  return dh;
}
 */
  
  template <class Arrangement, class LCC>
  typename LCC::Dart_handle arr2lcc(const Arrangement &arr, LCC& lcc)
  {
    lcc.clear();
    typename LCC::Dart_handle dh;
    //Vectors that are used to store vertex_attribute_handle
    //and points
    std::vector<typename LCC::Vertex_attribute_handle> vah;
    std::vector<typename LCC::Point> vert;
    //Iteration over the vertices and create associate vertex_attribute
    typename Arrangement::Vertex_const_iterator it;
    for (it = arr.vertices_begin(); it != arr.vertices_end(); ++it) {
      typename Arrangement::Traits_2::Point_2 p = it->point();
      typename LCC::Point v = typename LCC::Point(p.x(), p.y());
      typename LCC::Vertex_attribute_handle h = lcc.create_vertex_attribute(v);
      vert.push_back(v);
      vah.push_back(h);
    }
    //Iteration over the edges and create associate darts and links
    typename Arrangement::Edge_const_iterator eit;
    for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
      typename Arrangement::Halfedge_const_handle he = eit;
      typename Arrangement::Traits_2::Point_2 p1 = he->source()->point();
      typename Arrangement::Traits_2::Point_2 p2 = he->target()->point();
      typename LCC::Point v1 = typename LCC::Point(p1.x(), p1.y());
      typename LCC::Point v2 = typename LCC::Point(p2.x(), p2.y());
      int m = 0, n = 0;
      bool first = false, second = false;
      for (int i = 0; i < vert.size(); ++i) {
        if (vert[i] == v1) {
          m = i;
          first = true;
        }
        else if (vert[i] == v2) {
          n = i;
          second = true;
        }
        if (first && second) break;
      }
      dh = lcc.create_dart(vah[m]);
      typename LCC::Dart_handle dh2 = lcc.create_dart(vah[n]);
      lcc.template link_beta<2>(dh, dh2);
    }
    CGAL_assertion(dh!=NULL);
    return dh;
  }
    

  //Function template that transfer Linear_cell_complex to Arrangement
  template <class LCC, class Arrangement>
  typename Arrangement::Halfedge_handle lcc2arr(const LCC &lcc,
                                                Arrangement& arr)
  {
    arr.clear();
    typename Arrangement::Halfedge_handle he;
    //Vectors that used to store the vertex_handle and points
    std::vector<typename Arrangement::Vertex_handle> vvh;
    std::vector<typename LCC::Point> lp;
    //Generate all the points and insert associate Vertex_handle
    for (typename LCC::Vertex_attribute_range::const_iterator it =
        lcc.vertex_attributes().begin(), itend=lcc.vertex_attributes().end();
        it!=itend; ++it) {
      typename LCC::Point p = it->point();
      typename Arrangement::Traits_2::Point_2 p1 = typename
        Arrangement::Traits_2::Point_2(p.x(), p.y());
      typename Arrangement::Vertex_handle vh =
        arr.insert_in_face_interior(p1, arr.unbounded_face());
      vvh.push_back(vh);
      lp.push_back(p);
    }
    for (typename LCC::template One_dart_per_cell_range<1>::
        const_iterator it=lcc.template one_dart_per_cell<1>().begin(),
        itend=lcc.template one_dart_per_cell<1>().end();
        it!=itend; ++it) {
      //Get the two endpoints of each edge
      typename LCC::Point v1 = LCC::point(it);
      typename LCC::Point v2 = LCC::point(it->other_extremity());
      //Create associate points in Arrangement
      typename Arrangement::Traits_2::Point_2 p1 = typename
        Arrangement::Traits_2::Point_2(v1.x(), v1.y());
      typename Arrangement::Traits_2::Point_2 p2 = typename
        Arrangement::Traits_2::Point_2(v2.x(), v2.y());
      //Create this edge using the two points
      typename Arrangement::Traits_2::X_monotone_curve_2 c =
        typename Arrangement::Traits_2::X_monotone_curve_2(p1 ,p2);
      //Find the corresponding Vertex_handle
      int m = 0, n = 0;
      bool first = false, second = false;
      for (int i = 0; i < lp.size(); ++i) {
        if (lp[i] == v1) {
          m = i;
          first = true;
        }
        else if (lp[i] == v2) {
          n = i;
          second = true;
        }
        if (first && second) break;
      }
      //Insert the edge into the Arrangement
      he = arr.insert_at_vertices(c, vvh[m], vvh[n]);
    }
    CGAL_assertion(he!=NULL);
    return arr.edges_begin();
  }

}//namespace
#endif
