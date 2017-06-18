#ifndef HFAT_H
#define HFAT_H
#pragma once
#include<CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Aff_transformation_3.h>
#include<cstdlib>
namespace HexEx{

class Half_face_and_transition{
public:
  Half_face_and_transition(){};
  Half_face_and_transition(CGAL::Linear_cell_complex_for_generalized_map<3>::Dart_const_handle dh,
                           CGAL::Aff_transformation_3<CGAL::Exact_predicates_inexact_constructions_kernel> &tr){
    dart_handle = dh;

    min_transformation = tr;
  }
  CGAL::Linear_cell_complex_for_generalized_map<3>::Dart_const_handle dart_handle;
  CGAL::Aff_transformation_3<CGAL::Exact_predicates_inexact_constructions_kernel> min_transformation;  
};
}//namespace HexEx
#endif
