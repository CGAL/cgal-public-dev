//Function templates

#ifndef files_convertors_h
#define files_convertors_h
#include <iostream>
#include <vector>
//#include <CGAL/Combinatorial_map_operations.h>

namespace CGAL {
    /** @file convertors.h
     * 
     */

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
    
    
    
    
  //Iteration on the edge to get the coordinates of points
  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    if (!fit->is_unbounded()) {
      typename Arrangement::Ccb_halfedge_const_circulator circ =
        fit->outer_ccb();
      typename Arrangement::Ccb_halfedge_const_circulator curr = circ;
        
        typename LCC::Dart_handle dh;
        typename Arrangement::Traits_2::Point_2 origin;
        
      bool first = true;
        
      do {
        typename Arrangement::Halfedge_const_handle he = curr;
          hch.push_back(he);
    
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
            typename LCC::Dart_handle d0 = lcc.make_triangle(v1, v2, v3);
            dh = d0->beta(1)->beta(1);
            first = false;
          } else {
              typename Arrangement::Traits_2::Point_2 p1 = he->source()->point();
              typename Arrangement::Traits_2::Point_2 p2 = he->target()->point();
              typename LCC::Dart_handle d0 = lcc.make_triangle(origin, p1, p2);
              lcc.template sew<2>(dh, d0);
              dh = d0->beta(1)->beta(1);
              CGAL::template remove_cell<LCC, 1>(lcc,d0);
          }
      } while (++curr != circ);
    }
  }
    
    return th;
   // return dart[0];
};


template <class LCC, class Arrangement>
typename Arrangement::Halfedge_handle lcc2arr(const LCC &lcc, Arrangement& arr)
{
  typename Arrangement::Halfedge_handle he;
  //Iteration throuth the edges of the Linear_cell_complex
  for (typename LCC::template One_dart_per_cell_range<1>::
         const_iterator it=lcc.template one_dart_per_cell<1>
         ().begin(),
         itend=lcc.template one_dart_per_cell<1>
         ().end();
         it!=itend; ++it)    {
    //Get the two endpoints of this edge
    typename LCC::Point p1 = LCC::point(it);
    typename LCC::Point p2 = LCC::point(it->other_extremity());
    //Get the coordinates of the two endpoints
    typename LCC::FT sx1 = p1.x();
    typename LCC::FT sy1 = p1.y();
    typename LCC::FT sx2 = p2.x();
    typename LCC::FT sy2 = p2.y();
    //Create the associate points of Arrangement
    typename Arrangement::Traits_2::Point_2 v1 =
      typename Arrangement::Traits_2::Point_2(sx1, sy1);
    typename Arrangement::Traits_2::Point_2 v2 =
      typename Arrangement::Traits_2::Point_2(sx2, sy2);
    //Create the segment
    typename Arrangement::Traits_2::Segment_2 s =
      typename Arrangement::Traits_2::Segment_2(v1, v2);
    insert(arr, s);
  }
  return he;
}
}


#endif
