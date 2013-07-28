//Function templates
#ifndef files_convertors_h
#define files_convertors_h
#include <iostream>
#include <vector>

#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Combinatorial_map_operations.h>

namespace CGAL {
  /** @file convertors.h
   *
  */
  //Function template that transfer Arrangement to Linear_cell_complex
  template <class Arrangement, class LCC>
  typename LCC::Dart_handle arr2lcc(const Arrangement &arr, LCC& lcc)
  {
    lcc.clear();
    typename LCC::Dart_handle dh;
    //Vectors for future referfence
    std::vector<typename LCC::Vertex_attribute_handle> vah;
    std::vector<typename Arrangement::Traits_2::Point_2> arrp;
    std::vector<typename LCC::Point> vert;
    std::vector<typename Arrangement::Halfedge_const_handle> heh;
    std::vector<typename LCC::Dart_handle> dhd;
    //Iteration over the vertices and create associate vertex_attribute
    typename Arrangement::Vertex_const_iterator it;
    for (it = arr.vertices_begin(); it != arr.vertices_end(); ++it) {
      typename Arrangement::Traits_2::Point_2 p = it->point();
      typename LCC::Point v = typename LCC::Point(p.x(), p.y());
      typename LCC::Vertex_attribute_handle h = lcc.create_vertex_attribute(v);
      vert.push_back(v);
      vah.push_back(h);
      arrp.push_back(p);
    }
    //Iteraton over each edge to create all the connected darts with beta<1>
    typename Arrangement::Edge_const_iterator eit;
    for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
      typename Arrangement::Halfedge_const_handle he = eit;
      bool is_in = false, twin = false;
      //check whether this halfedge_handle has been reached
      for (int i = 0; i < heh.size(); ++i) {
        if (he == heh[i]) {
          is_in = true;
          break;
        }
      }
      //If it has been reached, check its twin halfedge
      if (is_in) {
        typename Arrangement::Halfedge_const_handle th = he->twin();
        for (int i = 0; i < heh.size(); ++i) {
          if (th == heh[i]) {
            twin = true;
            break;
          }
        }
        //If the halfedge has not been reached, iterate all the halfedge that
        //could be reached by next()
        if (!twin) {
          typename Arrangement::Halfedge_const_handle th2 = th;
          typename LCC::Dart_handle fp;
          do {
            heh.push_back(th);
            int pos = 0;
            for (int i = 0; i < arrp.size(); ++i) {
              if (th->source()->point() == arrp[i]) {
                pos = i;
                break;
              }
            }
            //Consider it is the first one of this loop
            if (th == th2) {
              dh = lcc.create_dart(vah[pos]);
              fp = dh;
              dhd.push_back(dh);
            }
            else {
              typename LCC::Dart_handle dh2 = lcc.create_dart(vah[pos]);
              dhd.push_back(dh2);
              //link by beta<1>
              lcc.template link_beta<1>(dh, dh2);
              dh = dh2;
            }
            th = th->next();
          } while (th != th2);
          //link the first one and the last one
          lcc.template link_beta<1>(dh, fp);
        }
      }
      //The halfedge has not been reached
      else {
        typename Arrangement::Halfedge_const_handle he2 = he;
        typename LCC::Dart_handle fp;
        do {
          heh.push_back(he);
          int pos = 0;
          for (int i = 0; i < arrp.size(); ++i) {
            if (he->source()->point() == arrp[i]) {
              pos = i;
              break;
            }
          }
          if (he == he2) {
            dh = lcc.create_dart(vah[pos]);
            fp = dh;
            dhd.push_back(dh);
          }
          else {
            typename LCC::Dart_handle dh2 = lcc.create_dart(vah[pos]);
            dhd.push_back(dh2);
            lcc.template link_beta<1>(dh, dh2);
            dh = dh2;
          }
          he = he->next();
        } while (he != he2);
        lcc.template link_beta<1>(dh, fp);
      }
    }
    //Iteration over all the edges, and link the dart by beta<2>
    for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
      typename Arrangement::Halfedge_const_handle he = eit;
      typename Arrangement::Halfedge_const_handle ht = he->twin();
      int m = 0, n = 0;
      bool first = false, second = false;
      //find out the associate positions of the two Dart_handle
      for (int i = 0; i < heh.size(); ++i) {
        if (heh[i] == he) {
          m = i;
          first = true;
        }
        else if (heh[i] == ht) {
          n = i;
          second = true;
        }
        if (first && second) break;
      }
      //link by beta<2>
      lcc.template link_beta<2>(dhd[m], dhd[n]);
    }
    CGAL_assertion(dh != NULL);
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
