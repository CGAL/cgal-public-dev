//Function templates

#ifndef files_convertors_h
#define files_convertors_h
#include <iostream>
#include <vector>
//Function template that transfer Arrangement to Linear_cell_complex
template <class Arrangement, class LCC>
typename LCC::Dart_handle arr2lcc(const Arrangement &arr, LCC& lcc)
{
  typename LCC::Dart_handle dh;
  typename Arrangement::Edge_const_iterator he;
    
    std::vector<typename Arrangement::Traits_2::Point_2> po;
    
    std::vector<typename LCC::Dart_handle> dart;
    
  //Iteration on the edge to get the coordinates of points
  for (he = arr.edges_begin(); he != arr.edges_end(); ++he) {
    typename Arrangement::Traits_2::Point_2 p1 = he->source()->point();
    typename Arrangement::Traits_2::Point_2 p2 = he->target()->point();
      
      po.push_back(p1);
      po.push_back(p2);
      
    //Get coordinates of the endpoints of the edge
    typename Arrangement::Traits_2::FT sx1 = p1.x();
    typename Arrangement::Traits_2::FT sy1 = p1.y();
    typename Arrangement::Traits_2::FT sx2 = p2.x();
    typename Arrangement::Traits_2::FT sy2 = p2.y();
    //Create associate points of Linear_cell_complex
    typename LCC::Point v1 =typename LCC::Point(sx1, sy1);
    typename LCC::Point v2 =typename LCC::Point(sx2, sy2);
    //Add the segment into the Linear_cell_complex object
    dh = lcc.make_segment(v1, v2);
      /*
      dart.push_back(dh);
      for (int i=0; i<po.size(); i++) {
          if (po[i]==p1 || po[i]==p2) {
              if (i % 2==0) {
                  lcc.template sew<1>(dart[i/2],dh);
              }
              else {
                  lcc.template sew<1>(dart[(i-1)/2], dh);
              }
              break;
          }
      }
       */
      
  }
  return dh;
}


template <class LCC, class Arrangement>
void lcc2arr(const LCC &lcc, Arrangement& arr)
{
  typename Arrangement::Halfedge_handle he;
      
  typename LCC::Dart_handle dart;
   // typename LCC::Dart d;
  std::vector<typename Arrangement::Traits_2::Point_2> p_vec;
  int count = 0;

    for(typename LCC::template One_dart_per_cell_range<1>::iterator it=lcc.template one_dart_per_cell_range<1>(dart).begin(), itend =lcc.template one_dart_per_cell_range<1>(dart).end(); it!=itend; ++it)
    {
        
        typename LCC::Point temp = typename LCC::Point(1,1);//it->point();
       // std::cout<<temp<<std::endl;
    ++count;
       // std::cout<<"count is "<<count<<std::endl;
    typename Arrangement::Traits_2::Point_2 p =
      typename Arrangement::Traits_2::Point_2(temp.x(), temp.y());
    p_vec.push_back(p);
      
    if (count %2 == 0) {
        typename Arrangement::Traits_2::Segment_2 s = typename Arrangement::Traits_2::Segment_2(p_vec[0], p_vec[1]);
      insert(arr, s);
      p_vec.clear();
    }
  }
}


#endif
