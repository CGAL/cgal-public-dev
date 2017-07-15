#ifndef GEO_EXTR_H
#define GEO_EXTR_H
#include"typedefs.h"
#include<vector>
#include<limits>
#include<algorithm>
#include"func.h"
#include"handles.h"

bool xfn(Point_3 i, Point_3 j){ return (i[0]<j[0]); }

bool yfn(Point_3 i, Point_3 j){ return (i[1]<j[1]); }

bool zfn(Point_3 i, Point_3 j){ return (i[2]<j[2]); }

void vertex_extraction(LCC_3& input_tet_mesh, LCC_3& output_mesh, std::vector<Aff_transformation> parametrization_matrices){
  //std::vector<Dart_handle> darts;
  for(LCC_3::One_dart_per_cell_range<3>::iterator it = (input_tet_mesh.one_dart_per_cell<3>()).begin(), itend = (input_tet_mesh.one_dart_per_cell<3>()).end(); it != itend; it++){
    double minx = std::numeric_limits<double>::max(), miny = std::numeric_limits<double>::max(), minz = std::numeric_limits<double>::max(); 
    double maxx = std::numeric_limits<double>::min(), maxy = std::numeric_limits<double>::min(), maxz = std::numeric_limits<double>::min();
    std::vector<Point_3> params;
    Dart_handle temp = it;
    Point_3 param = (input_tet_mesh.info(temp)).parameters; params.push_back(param);
    temp = input_tet_mesh.alpha(it,0);
    param = (input_tet_mesh.info(temp)).parameters; params.push_back(param);
    temp = input_tet_mesh.alpha(it, 0, 1, 0);
    param = (input_tet_mesh.info(temp)).parameters; params.push_back(param);
    temp = input_tet_mesh.alpha(it, 2, 0, 1, 0);
    param = (input_tet_mesh.info(temp)).parameters; params.push_back(param);
    std::cout<<params[0]<<" "<<params[1]<<" "<<params[2]<<" "<<params[3]<<std::endl;
    minx = (*std::min_element(params.begin(), params.end(), xfn))[0];
    maxx = (*std::max_element(params.begin(), params.end(), xfn))[0];
    miny = (*std::min_element(params.begin(), params.end(), yfn))[1];
    maxy = (*std::max_element(params.begin(), params.end(), yfn))[1];
    minz = (*std::min_element(params.begin(), params.end(), zfn))[2];
    maxz = (*std::max_element(params.begin(), params.end(), zfn))[2];
    std::cout<<"minx: "<<minx<<" maxx: "<<maxx<<std::endl;
    Aff_transformation at_inv = (parametrization_matrices[(input_tet_mesh.info(it)).cell_no]).inverse();

    Tetrahedron_3 tet(params[0], params[1], params[2], params[3]);
   
    for(int i = minx; i<= maxx; i++){
      for(int j = miny; j<=maxy; j++){
        for(int k = minz; k<=maxz; k++){
          Point_3 p(i,j,k);
          if(tet.has_on_bounded_side(p) || tet.has_on_boundary(p)){
            output_mesh.create_vertex_attribute(p.transform(at_inv));

          }
        }
      }
    }
  }
} 

#endif
