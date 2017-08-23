#include "hexextr.h"

Aff_transformation HexExtr::extract_transition_function(Dart_handle dh){ //the function returns the tranformation g applied to one tet to get to the vertices of the adjacent tet.
  Aff_transformation id(1,0,0,0,1,0,0,0,1,1);
    if(input_tet_mesh.is_free(dh, 3)){//if the tet is on the boundary, it doesn't have a neighbouring tet. Returns identity.
      return id;
    }
    else{// Has an adjacent tet sharing a common face    

      Dart_const_handle dh1 = dh; 

//dart handle corresponding to the same edge but in the adjacent tet
      Dart_const_handle dh2 = input_tet_mesh.beta<3>(dh1);


      std::vector<Point> face1, face2; 
      int i=0;
     
//adding the vertices of face1 (face of the first tet)
      face1.push_back(input_tet_mesh.point(dh1));
      face1.push_back(input_tet_mesh.point(input_tet_mesh.beta(dh1,1)));
      face1.push_back(input_tet_mesh.point(input_tet_mesh.beta(dh1,1,1)));

//adding the vertices of face2 (face of the second tet)
      face2.push_back(input_tet_mesh.point(dh2));
      face2.push_back(input_tet_mesh.point(input_tet_mesh.beta(dh2,1)));
      face2.push_back(input_tet_mesh.point(input_tet_mesh.beta(dh2,1,1)));

      if(face1[0] == face2[0] && face1[1] == face2[1] && face1[2] == face2[2]){
// transition function is identity.
        return id;
      }
      Vector_3 c1 = face1[1] - face1[0];
      Vector_3 d1 = face1[2] - face1[1];
      Vector_3 c2 = face2[1] - face2[0];
      Vector_3 d2 = face2[2] - face2[1];
      auto min_dist = std::numeric_limits<double>::max();

//the min transtion given by G[i]
      int min_trans_index = 0;

//to find the transition function from cubical chiral symmetry group G
      for(auto i = 0; i < 24; ++i){
        Vector_3 transf_c1 = G[i].transform(c1);

//Frobenius norm:
        auto dist1 = CGAL::scalar_product((c2 - transf_c1),(c2 - transf_c1)); 
        Vector_3 transf_d1 = G[i].transform(d1);

//Frobenius norm:
        auto dist2 = CGAL::scalar_product((d2 - transf_d1),(d2 - transf_d1));
        if(dist1 + dist2 < min_dist){ 
          min_dist = dist1+dist2;
          min_trans_index = i;
        }
      }
      Point new_point = G[min_trans_index].transform(face1[0]);

//rounding to integer translation
      Vector_3 t(std::round((face2[0])[0] - new_point[0]), std::round((face2[0])[1] - new_point[1]), std::round((face2[0])[2] - new_point[2]));     

//Adding translation to the transformation matrix.
      Aff_transformation final_transform_for_dh1(G[min_trans_index].m(0,0), G[min_trans_index].m(0,1), G[min_trans_index].m(0,2), t[0], G[min_trans_index].m(1,0), G[min_trans_index].m(1,1), G[min_trans_index].m(1,2), t[1], G[min_trans_index].m(2,0), G[min_trans_index].m(2,1), G[min_trans_index].m(2,2), t[2], 1);
     
       return final_transform_for_dh1;

  }
}


Aff_transformation HexExtr::get_parametrization_matrix(Point p, Point q, Point r, Point s, Point u, Point v, Point w, Point x){/**
* This function returns the parametrization matrix. When this matrix is applied to four 3D points forming a tet, we get the corresponding coordinates of the tet in parametric space.
* p, q, r, s are the coordinates of the original points of the tet
* u, v, w, x  are the coordinates of the points of the tet in parametric space
*/
  Aff_transformation at1(q[0]-p[0], r[0]-p[0],s[0]-p[0], p[0], q[1]-p[1], r[1]-p[1],s[1]-p[1], p[1], q[2]-p[2], r[2]-p[2],s[2]-p[2], p[2], 1);
  Aff_transformation at2(v[0]-u[0], w[0]-u[0], x[0]-u[0], u[0], v[1]-u[1], w[1]-u[1], x[1]-u[1], u[1], v[2]-u[2], w[2]-u[2], x[2]-u[2], u[2], 1);
  return (at2*(at1.inverse()));
}


void HexExtr::set_chiral_symmetry_aff_transformations(){//six directions in 3D to obtain the 24 affine transformations in chiral cubical symmetry group:
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
          G.push_back(Aff_transformation(directions[i].dx(), directions[j].dx(), directions[k].dx(),
directions[i].dy(), directions[j].dy(), directions[k].dy(), 
directions[i].dz(), directions[j].dz(), directions[k].dz(), 1)); //chiral cubical symmetry group

}

void HexExtr::preprocess(){

//It sets G to contains 24 affine transformations under chiral cubical symmetry.
  set_chiral_symmetry_aff_transformations();

//initialize dart_info to some values
  int cell = set_dart_info();

// a vector of size of the number of tets in the input mesh, to contain transition matrices of each tet to adjacent tets
  std::vector<std::vector<Aff_transformation>> g1(cell);
  g = g1;

/** initializing each element of g to be another vector of transition matrices.
* g[i][j] will correspond to the transition fucntion(matrix) applied to tet i to transform to tet j, where i, j are cell_no of corresponding tets.
*/
  for(int j = 0; j<cell; j++){
    std::vector<Aff_transformation> temp(cell);
    g[j] = temp;
  }

// A vector of affine transformations to contain parametrization matrices of each tet indexed by cell_no
  std::vector<Aff_transformation> pm(cell);  
  parametrization_matrices = pm;

//iterate through each tet:    
  for(LCC_3::One_dart_per_cell_range<3>:: iterator it = input_tet_mesh.one_dart_per_cell<3>().begin(), itend = input_tet_mesh.one_dart_per_cell<3>().end(); it != itend; it++){
    std::vector<Point> points, parameters;

//iterate through one dart corresponding to vertices of the given tet, so as to extract the coordinates in initial and parametrized space, and stored in the vectors 'points' and 'parameters'.
    for(LCC_3::One_dart_per_incident_cell_range<0,3>::iterator it2 = input_tet_mesh.one_dart_per_incident_cell<0,3>(it).begin(), it2end = input_tet_mesh.one_dart_per_incident_cell<0,3>(it).end(); it2 != it2end; it2++){
      points.push_back(input_tet_mesh.point(it2)); parameters.push_back((input_tet_mesh.info(it2)).parameters);
    }

//finding the parametrization matrix of the tet
    Aff_transformation at = get_parametrization_matrix(points[0], points[1], points[2], points[3], parameters[0], parameters[1], parameters[2], parameters[3]); 

       
// assigning the vector parametrization_matrices with the affine transformation indexed using cell_no
    parametrization_matrices[(input_tet_mesh.info(it)).cell_no] = at;

        
  }


  for(LCC_3::One_dart_per_cell_range<2>::iterator it = input_tet_mesh.one_dart_per_cell<2>().begin(), 
itend = input_tet_mesh.one_dart_per_cell<2>().end(); it != itend; it++){//iterating over all faces, we extract trasition functions for each of the two tets adjacent on a face.
    Aff_transformation at = extract_transition_function(it);
    g[(input_tet_mesh.info(it)).cell_no][(input_tet_mesh.info(input_tet_mesh.beta(it, 3))).cell_no] = at;
    g[(input_tet_mesh.info(input_tet_mesh.beta(it, 3))).cell_no][(input_tet_mesh.info(it)).cell_no] = at.inverse();
  }

}
