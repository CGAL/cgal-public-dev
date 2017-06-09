#ifndef HEXEXTR_H
#define HEXEXTR_H
#pragma once 
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Cell_attribute.h>
#include <CGAL/Aff_transformation_3.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Direction_3.h>
#include "triangulation_to_LCC.h"
#include "Half_face_and_transition.h"
#include "func.h"
#include <cmath>
#include <CGAL/Generalized_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K; 
typedef CGAL::Vector_3<K>                                   Vector_3;
typedef CGAL::Direction_3<K>                                Direction;
typedef CGAL::Linear_cell_complex_for_generalized_map<3>    LCC_3;
typedef LCC_3::Dart_handle                                  Dart_handle;
typedef LCC_3::Point                                        Point;
typedef LCC_3::Vertex_attribute_handle                      Vertex_attribute_handle;
typedef CGAL::Aff_transformation_3<K>                       Transformation;
typedef LCC_3::size_type				    size_type;		
typedef CGAL::Point_3<K>				    Point_3;
typedef HexEx::Half_face_and_transition                     Half_face_and_transition;
namespace HexEx{
class HexExtr{
  public:
    HexExtr();
    HexExtr(std::string infilename): identity(1,0,0,0,1,0,0,0,1,1){
      input_tet_mesh = load_off_to_LCC(infilename); //tetmesh to lcc
      directions.push_back(Direction(1,0,0));
      directions.push_back(Direction(0,1,0));
      directions.push_back(Direction(0,0,1));
      directions.push_back(Direction(-1,0,0));
      directions.push_back(Direction(0,-1,0));
      directions.push_back(Direction(0,0,-1));
      for(int i = 0;i < 6;i++) 
        for(int j = 0;j < 6;j++) 
          for(int k = 0;k < 6;k++) 
            if(CGAL::cross_product(directions[i].vector(),directions[j].vector()) == directions[k].vector()) 
              G.push_back(Transformation(directions[i].dx(), directions[i].dy(), directions[i].dz(), directions[j].dx(), directions[j].dy(), directions[j].dz(), directions[k].dx(), directions[k].dy(), directions[k].dz(), 1));

      for(LCC_3::One_dart_per_cell_range<3>::iterator it = input_tet_mesh.one_dart_per_cell<3>().begin(), itend = input_tet_mesh.one_dart_per_cell<3>().end(); it != itend; it++){
		all_faces_with_transitions.push_back(extract_transition_function(*it, input_tet_mesh, G));
      }

//Sanitization 
    
    }

    std::vector<Half_face_and_transition> all_faces_with_transitions;
    std::vector<std::vector<Point>> points_in_each_cell;
    std::vector<Direction> directions;
    Transformation identity;//(1,0,0,0,1,0,0,0,1,1);
    LCC_3 input_tet_mesh;
    std::vector<Transformation> G; //chiral cubical symmetry group

};
}

#endif
