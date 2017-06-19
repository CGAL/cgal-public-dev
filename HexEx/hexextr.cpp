#include"hexextr.h"
#include"handles.h"
HexExtr::HexExtr(std::string infilename): identity(1,0,0,0,1,0,0,0,1,1){
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
      for(LCC_3::One_dart_per_cell_range<3>::iterator it = input_tet_mesh.one_dart_per_cell<3>().begin(), 
itend = input_tet_mesh.one_dart_per_cell<3>().end(); it != itend; it++){
         Face_handle fh( this->input_tet_mesh, it, i); i++;
         Aff_transformation at = extract_transition_function(it, input_tet_mesh, G);
         std::cout<<i<<std::endl;
         print_aff_transformation(at);
         faces_with_transitions.emplace(fh, at);
		
      }

//Sanitization 
    
    }
