#ifndef HEXEXTR_H
#define HEXEXTR_H
#include"typedefs.h"
//#include"handles.h"
#include"triangulation.h"

//#include"triangulation_to_LCC.h"
#include<unordered_map>
#include<map>
#include<CGAL/Linear_cell_complex_constructors.h>
#include"functions.h"
#include"hexahedron_extraction.h"
#include"connection_extraction.h"
#include"sanitization.h"
//#include"dart_extraction.h"
#include<vector>
#include"frame_field.h"


class HexExtr{
  public:
  
// Currently,meshing commands are executed when the constructor is called. This will be changes later to incorporate the fucntions as methos of this class.
    HexExtr(std::string infilename): identity(1,0,0,0,1,0,0,0,1,1){
//input_tet_mesh to lcc
     load_off_to_LCC(infilename, input_tet_mesh);
   
if(DEBUG)std::cout<<"beginning"<<std::endl;
  std::ifstream in;
  input_tet_mesh.clear();
  in.open("triangulation");
  if (in.is_open())
  {
    in>>input_tet_mesh; //doesn't work: segfault - figure out why
    in.close();
  }
 if(DEBUG)std::cout<<"loaded to lcc"<<std::endl;
//input_tet_mesh.display_characteristics(std::cout); std::cout<<std::endl;

//a parametrization function to test the extraction:
      dummy_parametrize(input_tet_mesh); 
if(DEBUG)std::cout<<"parametrized"<<std::endl;
//We would like to directly input parametrized meshes: so I am trying to first export LCC with parameters stored in dart_info using << in the lines 44 - 49, then try taking such a file as an input to our method in lines 51-57. If this works, parametrized_LCC file can be directly used to test.

/*std::ofstream of;
of.open("parametrized_LCC");
if(of.is_open()){
  of<<input_tet_mesh;
  of.close();
}
if(DEBUG)std::cout<<"Checkpoint: parametrization exported to file"<<std::endl;
input_tet_mesh.clear();
  in.open("parametrized_LCC");
  if (in.is_open())
  {
    in>>input_tet_mesh; //doesn't work: segfault - figure out why
    in.close();
  }
if(DEBUG)std::cout<<"Checkpoint: parametrization imported from file"<<std::endl;
input_tet_mesh.display_characteristics(std::cout); std::cout<<std::endl;

*/



//six directions in 3D to obtain the 24 affine transformations in chiral cubical symmetry group:
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


      int cell = 0;

//iterate through each tet in the input mesh:
      for(LCC_3::One_dart_per_cell_range<3>::iterator it = input_tet_mesh.one_dart_per_cell<3>().begin(), itend = input_tet_mesh.one_dart_per_cell<3>().end(); it != itend; it++){ 

//iterate through all the darts of each input mesh, to assign a cell_no unique to darts of each tet. Two darts belonging to different tets will have different cell_no
        for(LCC_3::Dart_of_cell_range<3>::iterator it1 = input_tet_mesh.darts_of_cell<3>(it).begin(), it1end = input_tet_mesh.darts_of_cell<3>(it).end(); it1 != it1end; it1++){ 
          (input_tet_mesh.info(it1)).cell_no = cell;
        }
        cell++;
      }

// a vector of size of the number of tets in the input mesh, to contain transition matrices of each tet to adjacent tets
      std::vector<std::vector<Aff_transformation>> g(cell);


/** initializing each element of g to be another vector of transition matrices.
* g[i][j] will correspond to the transition fucntion(matrix) applied to tet i to transform to tet j, where i, j are cell_no of corresponding tets.
*/
      for(int j = 0; j<cell; j++){
        std::vector<Aff_transformation> temp(cell);
        g[j] = temp;
      }

// A vector of affine transformations to contain parametrization matrices of each tet indexed by cell_no
      std::vector<Aff_transformation> parametrization_matrices(cell);  

//iterate through each tet:    
      for(LCC_3::One_dart_per_cell_range<3>:: iterator it = input_tet_mesh.one_dart_per_cell<3>().begin(), itend = input_tet_mesh.one_dart_per_cell<3>().end(); it != itend; it++){
        std::vector<Point> points, parameters;

//iterate through one dart corresponding to vertices of the given tet, so as to extract the coordinates in initial and parametrized space, and stored in the vectors 'points' and 'parameters'.
        for(LCC_3::One_dart_per_incident_cell_range<0,3>::iterator it2 = input_tet_mesh.one_dart_per_incident_cell<0,3>(it).begin(), it2end = input_tet_mesh.one_dart_per_incident_cell<0,3>(it).end(); it2 != it2end; it2++){
          points.push_back(input_tet_mesh.point(it2)); parameters.push_back((input_tet_mesh.info(it2)).parameters);
        }

//finding the parametrization matrix of the tet
        Aff_transformation at = get_parametrization_matrix(points[0], points[1], points[2], points[3], parameters[0], parameters[1], parameters[2], parameters[3]); 

        //std::cout<<points[0]<<" "<<points[0]<<" "<<points[1]<<" "<<points[2]<<" "<<parameters[0]<<" "<<parameters[1]<<" "<<parameters[2]<<std::endl; //for debug purposes
        //print_aff_transformation(at); //for debugging

// assigning the vector parametrization_matrices with the affine transformation indexed using cell_no
        parametrization_matrices[(input_tet_mesh.info(it)).cell_no] = at;

//creating a cell handle. TODO: needed? Try to do away with cell_handle, face_handle, vertex_handle and edge_handle.
       // Cell_handle ch(it, points, parameters, at);
       // cells.push_back(ch);
        
      }


      for(LCC_3::One_dart_per_cell_range<2>::iterator it = input_tet_mesh.one_dart_per_cell<2>().begin(), 
itend = input_tet_mesh.one_dart_per_cell<2>().end(); it != itend; it++){
        // Face_handle fh(input_tet_mesh, it, i); i++;
         //dart_in_face.emplace(it, fh);        
        Aff_transformation at = extract_transition_function(it, input_tet_mesh, G);
        g[(input_tet_mesh.info(it)).cell_no][(input_tet_mesh.info(input_tet_mesh.beta(it, 3))).cell_no] = at;
        g[(input_tet_mesh.info(input_tet_mesh.beta(it, 3))).cell_no][(input_tet_mesh.info(it)).cell_no] = at.inverse();
         //std::cout<<i<<std::endl;
       // print_aff_transformation(at);
        // faces_with_transitions.emplace(fh, at);
         //faces.push_back(fh);
		
      }


     /* int v=0;
      for(LCC_3::Dart_range::iterator it = lcc.darts().begin(), itend = lcc.darts().end(); it != itend; it++){
        lcc.info(it) = v; v++;
      }*/
/*
      int nv = 0;
      for(LCC_3::One_dart_per_cell_range<0>::iterator it = input_tet_mesh.one_dart_per_cell<0>().begin(), itend = input_tet_mesh.one_dart_per_cell<0>().end(); it!=itend; it++){
        Vertex_handle vh(input_tet_mesh, it, nv, input_tet_mesh.is_free(it, 3));
        vertices.push_back(vh);
      }

      
      for(LCC_3::One_dart_per_cell_range<1>::iterator it = input_tet_mesh.one_dart_per_cell<1>().begin(), itend = input_tet_mesh.one_dart_per_cell<1>().end(); it != itend; it++){
        Vertex_handle from, to;
        for(std::vector<Vertex_handle>::iterator i = vertices.begin(), iend = vertices.end(); i != iend; i++){
          if(input_tet_mesh.point((*i).incident_dart) == input_tet_mesh.point(it)) from = (*i); //there must be a better way to implement this.
          if(input_tet_mesh.point((*i).incident_dart) == input_tet_mesh.point(input_tet_mesh.beta(it, 0))) to = (*i);
        }
        Edge_handle eh(input_tet_mesh, it, from, to); 
        edges.push_back(eh); 
      }     

      optimise_frame_field(input_tet_mesh, vertices, edges, 1); 

*/

//We have a parametrization matrix associated with each tet to transform a given set of points to parametric space. Since we have rounded off the translational terms in the extract_transition_function step

//Sanitization
if(DEBUG)std::cout<<"before sanitize"<<std::endl;
    sanitize(input_tet_mesh, g);
if(DEBUG)std::cout<<"after sanitize"<<std::endl;



//Extract vertices
//Extract darts
//making hexahedrons in the output mesh incorporting vertex extraction and dart extraction in a single step:
     if(DEBUG)std::cout<<"before extracting hexes"<<std::endl;
     extract_hexes(input_tet_mesh, output_mesh, parametrization_matrices, hex_handles, output_points);
     if(DEBUG)std::cout<<"after extracting hexes"<<std::endl;

     //extract_connections(output_mesh, hex_handles, output_points);
     output_mesh.sew3_same_facets();
     std::ofstream of;
     of.open("final_output.off");
     CGAL::write_off(output_mesh, of); 
     of.close();
     if(output_mesh.is_valid()) std::cout<<"Valid!"<<std::endl;
     else std::cout<<"invalid!"<<std::endl;
     std::cout<<"*****FINAL OUTPUT MESH*****"<<std::endl;
     output_mesh.display_characteristics(std::cout); std::cout<<std::endl;

    }
   // std::unordered_map<Face_handle, Aff_transformation> faces_with_transitions; //Take this as input and make dart_handle face_handle map using this
   // std::map<Dart_handle, Face_handle> dart_in_face;
    //Dart_handle ***connections;
    std::unordered_map<Point_3, Dart_handle> hex_handles;
    std::unordered_map<Point_3, Point_3> output_points;
   // std::vector<Face_handle> faces;
    //std::vector<Edge_handle> edges;
    //std::vector<Vertex_handle> vertices;
    //std::vector<Cell_handle> cells;
    //std::vector<Point> hvertices;
    std::vector<Direction> directions;
    Aff_transformation identity;//(1,0,0,0,1,0,0,0,1,1);
    LCC_3 input_tet_mesh, output_mesh;//, parametrized_mesh, 
    std::vector<Aff_transformation> G; //chiral cubical symmetry group
   
};

#endif
