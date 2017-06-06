#ifndef EACHCELL_H
#define EACHCELL_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Cartesian.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Aff_transformation_3.h>
#include <iostream>
#include <algorithm>
#include<vector>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K; 
typedef CGAL::Cartesian<double>::Vector_3 Vector_3;
typedef CGAL::Direction_3 Direction;
typedef CGAL::Linear_cell_complex_for_generalized_map<3> LCC_3;
typedef LCC_3::Dart_handle                                 Dart_handle;
typedef LCC_3::Point                                       Point;
typedef LCC_3::Vertex_attribute_handle                     VertexHandle;
typedef CGAL::Aff_transformation_3                         Transformation;
//define cell handle with a shared ptr 
//define face handle with a shared ptr
namespace HexEx{
class EachCell{

    Direction& dir();
    const Direction& dir() const;
    Parameter& parameter();
    const Parameter& parameter() const;



    VertexHandle& vertex();
    const VertexHandle& vertex() const;

    Dart_handle incidentHandle();

    Dart_handle& darthandle();  

};

class HexExtr{
  public:
    HexExtr();
    HexExtr(std::string infilename){
     //load the file into LCC
      identity(1,0,0,0,1,0,0,0,1,1);
      direction[0] = Direction(1,0,0);
      direction[1] = Direction(0,1,0);
      direction[2] = Direction(0,0,1);
      direction[3] = Direction(-1,0,0);
      direction[4] = Direction(0,-1,0);
      direction[5] = Direction(0,0,-1);
      for(int i = 0;i < 6;i++) 
        for(int j = 0;j < 6;j++) 
         for(int k = 0;k < 6;k++) 
           if(CGAL::cross_product(directions[i].vector(),directions[j].vector()) == directions[k].vector()) 
             transformations.push_back(Transformation(directions[i].dx(),directions[i].dy(),directions[i].dz(),directions[j].dx(),directions[j].dy(),directions[j].dz(),directions[k].dx(),directions[k].dy(),directions[k].dz(),1));  
    }//you have the transformations ready.



    vector<Direction> directions(6);
    Transformation identity;
    vector<Transformation> transformations;

}//namespace HexEx
#endif
