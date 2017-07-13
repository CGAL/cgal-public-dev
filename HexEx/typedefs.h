#ifndef TYPEDEFS_H
#define TYPEDEFS_H
#include<CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Aff_transformation_3.h>
#include<CGAL/Point_3.h>
#include<cassert>


struct myItem
{
  template < class GMap >
  struct Dart_wrapper
  {
    typedef double Dart_info;
    //typedef CGAL::Cell_attribute<GMap, int> Edge_attrib;
   // typedef CGAL::cpp11::tuple<void,Edge_attrib> Attributes;
    typedef CGAL::Cell_attribute_with_point< GMap, int, CGAL::Tag_true>   Vertex_attribute;
    
    typedef CGAL::cpp11::tuple<Vertex_attribute> Attributes;
  };
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel                    K; 
typedef CGAL::Linear_cell_complex_traits<3, K>	                               Traits;
typedef CGAL::Linear_cell_complex_for_generalized_map<3, 3, Traits, myItem>    LCC_3;
typedef LCC_3::Dart_handle                                                     Dart_handle;
typedef LCC_3::Dart_const_handle                                               Dart_const_handle;
typedef CGAL::Vector_3<K>                                                      Vector_3;
typedef CGAL::Direction_3<K>                                                   Direction;
typedef CGAL::Aff_transformation_3<K>                                          Aff_transformation;
typedef LCC_3::Point                                                           Point;  	
typedef CGAL::Point_3<K>			                               Point_3;
typedef CGAL::Tetrahedron_3<K>				                       Tetrahedron_3;
#endif
