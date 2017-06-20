#ifndef HEXEXTR_H
#define HEXEXTR_H
#include"typedefs.h"
#include"handles.h"
#include"triangulation_to_LCC.h"
#include<unordered_map>
#include<map>
#include"func.h"
#include<vector>

//class Face_handle;

namespace std{
  int dart_count = 0;
  template<>
  struct hash<Face_handle>{
    std::size_t operator()(const Face_handle& fh) const{
      return (fh.enumeration)%1000;
    } 
  };

  template<>
  struct hash<Dart_handle>{
    std::size_t operator()(const Dart_handle& dh) const{
      dart_count++;
      return (dart_count)%1000;
    } 
  };
}

class HexExtr{
  public:
    //HexExtr();
    HexExtr(std::string infilename): identity(1,0,0,0,1,0,0,0,1,1){
      load_off_to_LCC(infilename, input_tet_mesh); //tetmesh to lcc
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
              G.push_back(Aff_transformation(directions[i].dx(), directions[j].dx(),  directions[k].dx(),directions[i].dy(), directions[j].dy(), directions[k].dy(), directions[i].dz(), directions[j].dz(), directions[k].dz(), 1));
      int i = 1;
      for(LCC_3::One_dart_per_cell_range<2>::iterator it = input_tet_mesh.one_dart_per_cell<2>().begin(), 
itend = input_tet_mesh.one_dart_per_cell<2>().end(); it != itend; it++){
         Face_handle fh( this->input_tet_mesh, it, i); i++;
         dart_in_face.emplace(it, fh);        
         Aff_transformation at = extract_transition_function(it, input_tet_mesh, G);
         //std::cout<<i<<std::endl;
        // print_aff_transformation(at);
         faces_with_transitions.emplace(fh, at);
         faces.push_back(fh);
		
      }

//Sanitization 
    
    }
    std::unordered_map<Face_handle, Aff_transformation> faces_with_transitions; //TAKE this as input and make dart_hanld eafce_handle map fromthis info
    std::map<Dart_handle, Face_handle> dart_in_face;
    std::vector<Face_handle> faces;
    std::vector<Direction> directions;
    Aff_transformation identity;//(1,0,0,0,1,0,0,0,1,1);
    LCC_3 input_tet_mesh;
    std::vector<Aff_transformation> G; //chiral cubical symmetry group
   
};

#endif
