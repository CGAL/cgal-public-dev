#ifndef TYPEDEFS_H
#define TYPEDEFS_H
#include<CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Aff_transformation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K; 
typedef CGAL::Linear_cell_complex_for_generalized_map<3>    LCC_3;
typedef LCC_3::Dart_handle                                  Dart_handle;
typedef LCC_3::Dart_const_handle                            Dart_const_handle;
typedef CGAL::Vector_3<K>                                   Vector_3;
typedef CGAL::Direction_3<K>                                Direction;
typedef CGAL::Aff_transformation_3<K>                       Aff_transformation;
typedef LCC_3::Point                                        Point;  					
#endif
