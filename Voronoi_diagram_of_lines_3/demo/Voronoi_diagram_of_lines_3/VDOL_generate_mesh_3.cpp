#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1 

// ----------------------------------------------------
// includes for VOL_3
#include <CGAL/VDOL_3/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/macros.h>

#include <CGAL/Arithmetic_kernel.h>

// ----------------------------------------------------
// include from Mesh_3
#include "debug.h" // some flags for Mesh_3 
#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_3/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

// -----------------------------------------------------
// includes from /demo/VDOL_3/include 
#include <CGAL/Mesh_criteria_3_with_balls.h>
#include <CGAL/Voronoi_diagram_of_lines_3.h> 
#include <CGAL/VDOL_to_labeled_function_wrapper_3.h>

#include <CGAL/VDOL_3/io.h> 

typedef CGAL::Arithmetic_kernel::Rational               Rational;
typedef CGAL::Cartesian<Rational>                       EK_3;
typedef CGAL::Voronoi_diagram_of_lines_3<EK_3>          VDOL_3;
typedef VDOL_3::Line_3                                  ELine_3;
typedef VDOL_3::Point_3                                 EPoint_3;

struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
// struct K: public CGAL::Exact_predicates_exact_constructions_kernel {};


// Domain: DB
typedef CGAL::VDOL_to_labeled_function_wrapper_3<VDOL_3, K> VDOL_labeled_function_3;
typedef CGAL::Mesh_3::Labeled_mesh_domain_3<VDOL_labeled_function_3, K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
//typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef CGAL::Mesh_criteria_3_with_balls<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;


int main(int argc, char **argv)
{  
  std::cout << "THIS IS A NEW COMPILATION " << std::endl;

// ---------------------------------------- 
// the voronoi cell of one line 
  
  //CGAL::set_pretty_mode(std::cout);
  //CGAL::set_pretty_mode(std::cerr);
  
  // Get the name of the input file from the command line, or use the default
  // last.dat file if no command-line parameters are given.
  const char * filename = (argc >= 2) ? argv[1] : "last.dat"; 
  if(argc > 1){
    CGAL::VDOL_3::copy_file_into(filename,"last.dat");
  }
  std::cout << "Run test with file: " << filename << std::endl;

  std::vector<ELine_3> lines;
  CGAL::VDOL_3::read_lines_3(EK_3(),filename,std::back_inserter(lines));
  if(lines.size()==0) 
    std::cerr << "Failed to read file " << filename << std::endl;
  
  std::cout << "Number of lines (incl base line): "<< lines.size() << std::endl;
  
 
  // boost::shared_ptr<VDOL_3> vdol_3;
  // vdol_3 =  boost::shared_ptr<VDOL_3>(new VDOL_3(lines.begin(),lines.end()));
  
  // Voronoi diagram of lines 3 
  // VDOL_3 vdol_3(lines.begin(),lines.end());

  // labeled function 
  VDOL_labeled_function_3 vdol_labeled_function_3(lines.begin(),lines.end());
  // VDOL_labeled_function_3 vdol_labeled_function_3(vdol_3);
   
  // mesh domain 
  Mesh_domain domain(vdol_labeled_function_3,vdol_labeled_function_3.bounding_sphere());

 
//-------------------------------------------------------------------------------------------------------
// Input mesh parameters /
  double facet_angular_bound;
  double facet_size;
  double facet_distance_approximation;
  double tet_radius_edge_ratio = 5.;
  double tet_size;
  
  if( argc <= 2 )
    {
      std::cout << "READING DEFAULT mesh_parameters" << std::endl;
      std::string mesh_parameters_filename="last.mp";
      std::ifstream input_parameters(mesh_parameters_filename.c_str());
      input_parameters>> facet_angular_bound 
                      >> facet_size
                      >> facet_distance_approximation
                      >> tet_radius_edge_ratio
                      >> tet_size ;
      input_parameters.close();
      std::cout << "parameters image file name :" << mesh_parameters_filename.c_str() << "\n";
    }
  else if( argc == 3)
    {
      std::string mesh_parameters_filename=argv[2];
      std::cout << "READING "<< mesh_parameters_filename.c_str() << std::endl;
      CGAL::VDOL_3::copy_file_into(mesh_parameters_filename.c_str(), "last.mp");
      std::ifstream input_parameters(mesh_parameters_filename.c_str());
      input_parameters>> facet_angular_bound 
                      >> facet_size
                      >> facet_distance_approximation
                      >> tet_radius_edge_ratio
                      >> tet_size ;
      input_parameters.close();
      std::cout << "parameters image file name :" << mesh_parameters_filename.c_str() << "\n";
    }
  else if( argc >3) 
    {
      facet_angular_bound = atof(argv[2]);
      facet_size = atof(argv[3]);
      facet_distance_approximation = atof(argv[4]);
      tet_radius_edge_ratio = atof(argv[5]);
      tet_size = atof(argv[6]);
      std::ofstream os("last.mp");
      os << facet_angular_bound <<" " << facet_size <<" " << facet_distance_approximation<<" " <<tet_radius_edge_ratio <<" " << tet_size << std::endl;
      
    }
  else 
    {
      std::cout <<"Enter the angular bound for surface facets\n"
        "example : 30means 30 degree \n ";
       std::cin >> facet_angular_bound;
		
      std::cout <<"Enter the size criterion for facets\n"
        "(0 means ignoring the size criterion): ";
      std::cin >> facet_size;
      assert(std::cin);
      std::cout <<"Enter the \"center-center distance\" criterion for facets\n"
        "(0 means ignoring the distance size criterion): ";
      std::cin >> facet_distance_approximation;
      assert(std::cin);
		
      std::cout <<"Enter the radius-edge ratio for tetrahedra\n"
        "(0 means ignoring the uniform size criterion): ";
      std::cin >>tet_radius_edge_ratio;
      assert(std::cin);					
	
      std::cout <<"Enter the size criterion for tetrahedra\n"
        "(0 means ignoring the size criterion): ";
      std::cin >> tet_size;
      assert(std::cin);
    }

  {	
    std::cout << "facet_angular_bound          " << facet_angular_bound << std::endl; 
    std::cout << "facet_size                   " << facet_size << std::endl; 
    std::cout << "facet_distance_approximation " << facet_distance_approximation  << std::endl; 
    std::cout << "tet_radius_edge_ratio        " << tet_radius_edge_ratio << std::endl; 
    std::cout << "tet_size                     " << tet_size  << std::endl;
  }

  //exit(1);
//--------------------------------------------------------------------------------------------

// Mesh criteria
  Facet_criteria facet_criteria(facet_angular_bound, facet_size, facet_distance_approximation); // angle, size, approximation
  Cell_criteria cell_criteria(tet_radius_edge_ratio, tet_size); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);


  C3t3 c3t3; 
  C3t3::Triangulation& tr = c3t3.triangulation();
  std::cout<<"IS VALID "  << tr.is_valid() << std::endl;
//------------------------------------------------------------------------
//  random_shuffle(labeled_function_3.wpoints().begin(),labeled_function_3.wpoints().end());
  
  // TODO make this input parameter 
  double sradius_bounding_sphere = 110;  
  double sradius_clipping_sphere = 100;
  double sdistance_generation     = 0.2; 
  double sdistance_final          = 0.41;
  double point_weight            = 0.4;   
  
  CGAL::VDOL_3::Approximation_info approximation_info;
  vdol_labeled_function_3.approximation_info().sradius_bounding_sphere = sradius_bounding_sphere;
  vdol_labeled_function_3.approximation_info().sradius_clipping_sphere = sradius_clipping_sphere;
  vdol_labeled_function_3.approximation_info().sdistance_generation = sdistance_generation; 
  vdol_labeled_function_3.approximation_info().sdistance_final = sdistance_final; 
  vdol_labeled_function_3.approximation_info().point_weight = point_weight; 
  
  // save points to file   
  std::vector<VDOL_labeled_function_3::Point_3> points; 
  vdol_labeled_function_3.generate_points(std::back_inserter(points));

  {
    std::ofstream  os ("wpoints.dat");
    os << points.size() << std::endl;
    for(int i=0; i<points.size(); i++){
      os << points[i] << " " <<  point_weight << std::endl;
    }
  }
  
  typedef C3t3::Triangulation::Point Weighted_point;        
  for(int i=0; i<points.size(); i++){
    Weighted_point wp(points[i],point_weight);	
    C3t3::Vertex_handle  v = tr.insert(wp);
    c3t3.set_dimension(v,1);
  }
 

  std::cout<<"IS VALID "  << tr.is_valid() << std::endl;
  std::cout<<"nb initial vertices = "<<tr.number_of_vertices()<<std::endl;
  std::cout<<"nb initial facet = "<<tr.number_of_facets()<<std::endl;
  std::cout<<"nb initial cells = "<<tr.number_of_cells()<<std::endl;

//------------------------------------------------------------------------
// Meshing
  CGAL::refine_mesh_3(c3t3,domain, criteria, CGAL::parameters::no_exude(), CGAL::parameters::no_perturb());

//-------------------------------------------------------------------------------------------
//Outputs to .MESH
  {	
    std::stringstream s;
    std::string input_filename("result_with_patches.mesh");
    s << input_filename;
    // <<"_"<<facet_angular_bound<<"_"
    // <<facet_size<<"_"<<facet_distance_approximation<<"_"<<tet_radius_edge_ratio<<"_"<<tet_size<< ".mesh";
	
    std::cout << "Saving to file " << s.str()<<std::endl;
    std::ofstream file_medit(s.str().c_str());
    bool show_patches = true; 
    c3t3.output_to_medit(file_medit, true, show_patches);
  }
//   {	
//     std::stringstream s;
//     std::string input_filename("result_no_patches.mesh");
//     s << input_filename;
//     // <<"_"<<facet_angular_bound<<"_"
//     // <<facet_size<<"_"<<facet_distance_approximation<<"_"<<tet_radius_edge_ratio<<"_"<<tet_size<< ".mesh";
	
//     std::cout << "Saving to file " << s.str()<<std::endl;
//     std::ofstream file_medit(s.str().c_str());
//     bool show_patches = false; 
//     c3t3.output_to_medit(file_medit, true, show_patches);
//   }
//   {	
//     std::stringstream s;
//     std::string input_filename("result");
//     s << input_filename
//       <<"_"<<facet_angular_bound<<"_"
//       <<facet_size<<"_"<<facet_distance_approximation<<"_"<<tet_radius_edge_ratio<<"_"<<tet_size<< ".mesh";
  
//     std::cout << "Saving to file " << s.str()<<std::endl;
//     std::ofstream file_medit(s.str().c_str());
//     c3t3.output_to_medit(file_medit);
//   }
//----------------------------------------------------------------------------------------
//Outputs to .OFF => ne marche plus encore avec c3t3!
// {
// 	std::stringstream s;
// 	s << input_filename << "_output.off";
// 	std::ofstream file_off(s.str().c_str());
// 	CGAL::output_surface_facets_to_off(file_off, c2t3, 0); // 0 used to force the orientation of the facets
// 	std::cout << "Saving to file " << s.str()<<std::endl;
// 	//CGAL::output_oriented_surface_facets_to_off(file_off, tr);
// }
//----------------------------------------------------------------------------------------
// Outputs to .TRIAN or .OBJ or .OFF => ne marche plus !!!
// {	
// 	cout << "Saving surfaces... \n"; 
// 	Label lab ;
// 	for(unsigned int j=0; j<I.vector_of_used_labels().size();j++)
// 	{
// 		lab = I.vector_of_used_labels()[j];
// 		if(lab!=0)
// 		{
// 			std::stringstream s;
// 			s << input_filename/*<<"_PREV"*/<<"_"<<facet_angular_bound<<"_"
// 			<<facet_size<<"_"<<facet_distance_approximation<<"_"<<tet_size<< ".trian";
// 
// 			std::cout << "Saving to file " << s.str()<<std::endl;
// 			
// 			save_TRIAN(c2t3.triangulation(),s.str().c_str(),lab);
// 			//save_OBJ(tr,s.str().c_str(),lab);
// 			//save_OFF(tr,s.str().c_str(),lab);
// 		}
// 	}
// }
//----------------------------------------------------------------------------------------

  return 0;
}
