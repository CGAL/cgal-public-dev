#ifndef TYPEDEFS_H
#define TYPEDEFS_H
#include<CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Aff_transformation_3.h>
#include<CGAL/Point_3.h>
#include<cassert>
#include <CGAL/Bbox_3.h>
#define DEBUG 0
typedef struct dartinfo{

//all the darts corresponding to the same tet will have the same cell_no. This is to make it easier to assign a parametrization matrix, which can be determined from any dart of the tet.
  double cell_no;
  bool flipped;
//a variable to store if the particular dart is singular.
  bool singular;
  int singular_edges;
//parameters is a point corresponding to the point of gven dart in input mesh in the parametrized space 
  CGAL::Point_3<CGAL::Exact_predicates_inexact_constructions_kernel> parameters;
}dart_info;

struct myItem
{
  template < class CMap >
  struct Dart_wrapper
  {
    typedef dart_info Dart_info; 
    //typedef CGAL::Cell_attribute<GMap, int> Edge_attrib;
   // typedef CGAL::cpp11::tuple<void,Edge_attrib> Attributes;
    typedef CGAL::Cell_attribute_with_point< CMap, int, CGAL::Tag_true>   Vertex_attribute;
    
    typedef CGAL::cpp11::tuple<Vertex_attribute> Attributes;
  };
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel                    K; 
typedef CGAL::Linear_cell_complex_traits<3, K>	                               Traits;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3, 3, Traits, myItem>  LCC_3;
typedef LCC_3::Dart_handle                                                     Dart_handle;
typedef LCC_3::Dart_const_handle                                               Dart_const_handle;
typedef CGAL::Vector_3<K>                                                      Vector_3;
typedef CGAL::Direction_3<K>                                                   Direction;
typedef CGAL::Aff_transformation_3<K>                                          Aff_transformation;
typedef LCC_3::Point                                                           Point;  	
typedef CGAL::Point_3<K>			                               Point_3;
typedef CGAL::Tetrahedron_3<K>				                       Tetrahedron_3;
typedef CGAL::Bbox_3                                                           Bbox_3;

namespace std{ //TODO: is this needed? 
  /*int dart_count = 0;
  template<>
  struct hash<Face_handle>{
    std::size_t operator()(const Face_handle& fh) const{
      return (fh.enumeration)%1000;
    } 
  };*/
  template<>
  struct hash<Point_3>{
    std::size_t operator()(const Point_3& p) const{
      return (p.x()+p.y()+p.z())/p.x();
    } 
  };
}


#endif
