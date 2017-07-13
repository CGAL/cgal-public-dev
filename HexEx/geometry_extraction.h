#ifndef GEO_EXTR_H
#define GEO_EXTR_H
#include"typedefs.h"
#include<vector>
#include<limits>
#include"func.h"
//#include"hexextr.h"
#include"handles.h"

void vertex_extraction(LCC_3& input_tet_mesh, LCC_3& output_mesh){
  //std::vector<Dart_handle> darts;
  for(LCC_3::One_dart_per_cell_range<3>::iterator it = (input_tet_mesh.one_dart_per_cell<3>()).begin(), itend = (input_tet_mesh.one_dart_per_cell<3>()).end(); it != itend; it++){
    double minx = std::numeric_limits<double>::max(), miny = std::numeric_limits<double>::max(), minz = std::numeric_limits<double>::max(); 
    double maxx = std::numeric_limits<double>::min(), maxy = std::numeric_limits<double>::min(), maxz = std::numeric_limits<double>::min();
    std::vector<Point> points;
    for(LCC_3::One_dart_per_incident_cell_range<0,3>::iterator it1 = input_tet_mesh.one_dart_per_incident_cell<0,3>(it).begin(), it1end = input_tet_mesh.one_dart_per_incident_cell<0,3>(it).end(); it1 != it1end;it1++){
      Point p = input_tet_mesh.point(it1);
      points.push_back(p); //make sure you push the parameterized points.
      Vector_3 param;// = extract_parameters(p);
      minx = (param[0]<minx)?param[0]:minx;
      miny = (param[1]<miny)?param[1]:miny;
      minz = (param[2]<minz)?param[2]:minz;
      maxx = (param[0]>maxx)?param[0]:maxx;
      maxy = (param[1]>maxy)?param[1]:maxy;
      maxz = (param[2]>maxz)?param[2]:maxz;
    }

    minx = std::round(minx); miny = std::round(miny); minz = std::round(minz);
    maxx = std::round(maxx); maxy = std::round(maxy); maxz = std::round(maxz);
    Tetrahedron_3 tet(points[0], points[1], points[2], points[3]);
   
    for(double i = minx; i<= maxx; i++){
      for(double j = miny; j<=maxy; j++){
        for(double k = minz; k<=maxz; k++){
          Point_3 p(i,j,k);
          if(tet.has_on_bounded_side(p) || tet.has_on_boundary(p)){
           // (h.hvertices).push_back(p);
            output_mesh.create_vertex_attribute(p);
            //Dart_handle dh = (h.output_mesh).create_dart(p); //for the original points, not parametrized. 
            //darts.push_back(dh); 
          }
        }
      }
    }
  }
} 

#endif
