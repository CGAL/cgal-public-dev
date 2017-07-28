#ifndef HANDLES_H
#define HANDLES_H

#include"typedefs.h"
#include<vector>
#include<unordered_map>
#include<map>
//#include"hexextr.h"
class Vertex_handle{
  public:
      Vertex_handle(){}
      Vertex_handle(LCC_3 &lcc, Dart_handle &dh, int e, bool f):frame(0,0,1){
      incident_dart = dh;
      current_point = lcc.point(dh);
      next_point = lcc.point(lcc.beta(dh, 0));
      enumeration = e;
      boundary = f;
      //frame[0] = 0; frame[1] = 0; frame[2] = 1;
    }
    Dart_handle incident_dart;
    Point current_point;
    Point next_point;
    Vector_3 frame;
    int enumeration;
    bool boundary;
};

class Edge_handle{
  public:
    Edge_handle(LCC_3 &lcc, Dart_handle &dh, Vertex_handle& f, Vertex_handle& t){
      associated_dart = dh;
      from = f;
      to = t; 
    }
    Dart_handle associated_dart;
    Vertex_handle from;
    Vertex_handle to;
    //int enumeration;
    //Edge_handle next_edge;
};


class Face_handle{
  public:
    Face_handle(LCC_3 &lcc, Dart_handle& dh, int i){
      a_dart = dh;
      for(LCC_3::Dart_of_cell_range<3>:: iterator it=lcc.darts_of_cell<3>(dh).begin(), 
itend=lcc.darts_of_cell<3>(dh).end(); it!=itend; ++it){
        incident_darts.push_back(it);
      }
      enumeration = i;
    }
    std::vector<Dart_handle> incident_darts;
    Dart_handle a_dart;
    int enumeration;
    
    bool operator==(const Face_handle &other) const{ 
      return (enumeration == other.enumeration);
    }
};

class Cell_handle{
  public:
    Cell_handle(Dart_handle& dh, std::vector<Point>& p, std::vector<Point>& q, Aff_transformation& at){
      a_dart = dh;
      points = p;
      parameters = q;
      parametrization_matrix = at;
    }
    Dart_handle a_dart;
    //int enumeration;
    std::vector<Point> points, parameters;
    Aff_transformation parametrization_matrix;
};

bool comp(Vertex_handle i,Vertex_handle j){return !(i.boundary == false && j.boundary == true);}

/*
class Cell_handle{

};
*/
#endif
