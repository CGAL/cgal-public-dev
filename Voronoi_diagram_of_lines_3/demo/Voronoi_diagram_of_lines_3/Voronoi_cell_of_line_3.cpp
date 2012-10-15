#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1 

// ----------------------------------------------------
// includes for VOL_3
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/macros.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_curve_kernel_2_generator.h> // generates a default kernel 
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h> // traits for Arr_2

#include <CGAL/VDOL_3/Single_voronoi_cell_envelope_traits_3.h>

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
#include <CGAL/Voronoi_cell_of_line_3.h>




// -------------------------------------------------------------
// -----------------------

// --------------------------------------
// Typedefs 


// get default arithmetic 
typedef CGAL::Arithmetic_kernel AK;
typedef AK::Integer Integer;
typedef AK::Rational Rational;

// define linear kernel 
#if 1
typedef CGAL::Cartesian<Rational> EK_3;
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
#else
typedef CGAL::Cartesian<Rational> EK_3;
typedef EK_3 K;
// struct K: public CGAL::Exact_predicates_exact_constructions_kernel {};
#endif

typedef EK_3::Line_3                                    ELine_3;
typedef EK_3::Point_3                                   EPoint_3;

// define algebraic kerenel 
typedef CGAL::Algebraic_kernel_d_1_generator<Integer>     GEN_AK_1;
typedef GEN_AK_1::Default_algebraic_kernel_1            AK_1;
typedef CGAL::Algebraic_curve_kernel_2<AK_1>            AK_2;

// define curved kernel 
typedef CGAL::Curved_kernel_via_analysis_2< AK_2 >      CK_2;

// define traits for lower envelop diagram 
typedef CGAL::VDOL_3::
Single_voronoi_cell_envelope_traits_3<EK_3, CK_2>       SVCET_3;


// Domain: DB
typedef CGAL::Voronoi_cell_of_line_3<SVCET_3,K> VCL_3;
typedef CGAL::Mesh_3::Labeled_mesh_domain_3<VCL_3, K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
//typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef CGAL::Mesh_criteria_3_with_balls<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;





void copy_file_into(const char *infilename, const char *outfilename){
  std::ifstream     in_file (infilename);
  std::ofstream     out_file (outfilename);
  std::cout << "save " << infilename << " to " << outfilename << "... " << std::flush;
  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << infilename << "!" << std::endl;
    return; 
  }
  if (! out_file.is_open()) {
    std::cerr << "Failed to open " << outfilename << "!" << std::endl;
    return; 
  }
  out_file << in_file.rdbuf();
  std::cout << "  done "<< std::endl;
}

template <typename OutputIterator>
OutputIterator read_lines(const char *filename, OutputIterator oi){

  // Open the input file.
  std::ifstream     in_file (filename);
  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << filename << "!" << std::endl;
    return oi;
  }
  
  // Read the lines from the file
  // The input file format should be (all coordinate values are integers):
  // <n>                                       // number of lines.
  // <a1_x> <a1_y> <a1_z> <b1_x> <b1_y> <b1_z> // line #1.
  // <a2_x> <a2_y> <a2_z> <b2_x> <b2_y> <b2_z> // line #2.
  //   :      :       :      :
  // <an_x> <an_y> <an_z> <bn_x> <bn_y> <bn_z> // line #n.
  
  // read number of lines 
  unsigned int n; in_file >> n;

  for (int k = 0; k < n; ++k) {
    int a,b,c,d,e,f;
    in_file >> a >> b >> c>> d >> e >> f;
    typedef Rational FT; 
    *oi++ = ELine_3(EPoint_3(FT(a),FT(b),FT(c)),EPoint_3(FT(d),FT(e),FT(f)));
  }
  in_file.close();
}



int main(int argc, char **argv)
{
// ---------------------------------------- 
// the voronoi cell of one line 
  


  //CGAL::set_pretty_mode(std::cout);
  //CGAL::set_pretty_mode(std::cerr);
  
  // Get the name of the input file from the command line, or use the default
  // last.dat file if no command-line parameters are given.
  const char * filename = (argc >= 2) ? argv[1] : "last.dat"; 
  if(argc > 1){
    copy_file_into(filename,"last.dat");
  }
  std::cout << "Run test with file: " << filename << std::endl;

  std::vector<SVCET_3::Surface_3> lines;
  read_lines(filename,std::back_inserter(lines));
  if(lines.size()==0) 
    std::cerr << "Failed to read file " << filename << std::endl;
  
  std::cout << "Number of lines (incl base line): "<< lines.size() << std::endl;
  
 
 
  typedef CGAL::VDOL_3::Event_at_discontinuity_exception  Event_at_discontinuity_exception;
  boost::shared_ptr<VCL_3> vcl_3;
  
  int seed = -1; 
  bool done = false; 
  while(!done){
    std::cout << "try ..." <<std::flush; 
    try{
      vcl_3 =  boost::shared_ptr<VCL_3>(new VCL_3(seed,lines.begin(),lines.end()));
      done = true ;
    }catch(Event_at_discontinuity_exception& e){
      std::cout << "catch:" << std::endl 
                << e.what() << std::endl;
      seed++; 
      if(seed > 10) {
        std::cout << "infinite try/catch unexpected special case ? " 
                  << std::endl;
        assert(false);
      }
    } 
  }
    



//-------------------------------------------------------------------------------------------------------
// Input mesh parameters 
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
      copy_file_into(mesh_parameters_filename.c_str(), "last.mp");
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
// Domain
  Mesh_domain domain(*vcl_3,vcl_3->bounding_sphere());

// Mesh criteria
  Facet_criteria facet_criteria(facet_angular_bound, facet_size, facet_distance_approximation); // angle, size, approximation
  Cell_criteria cell_criteria(tet_radius_edge_ratio, tet_size); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);


  C3t3 c3t3; 
  C3t3::Triangulation& tr = c3t3.triangulation();
  std::cout<<"IS VALID "  << tr.is_valid() << std::endl;
//------------------------------------------------------------------------
//  random_shuffle(vcl_3->wpoints().begin(),vcl_3->wpoints().end());

  // save points to file   
  {
    std::ofstream  os ("wpoints.dat");
    os << vcl_3->wpoints().size() << std::endl;
    for(int i=0; i<vcl_3->wpoints().size(); i=i+1){
      os << vcl_3->wpoints()[i].first << " " <<  vcl_3->wpoints()[i].second << std::endl;
    }
  }
  
  
  typedef C3t3::Triangulation::Point Weighted_point;        
  for(int i=0; i<vcl_3->wpoints().size(); i=i+1){
    Weighted_point wp(vcl_3->wpoints()[i].first, vcl_3->wpoints()[i].second);	
//    std::cout<<"POINT    "  << wp << std::endl;
//     std::cout<<"IS VALID "  << tr.is_valid() << std::endl;
    C3t3::Vertex_handle  v = tr.insert(wp);
    c3t3.set_dimension(v,1);
  }
 

  std::cout<<"IS VALID "  << tr.is_valid() << std::endl;
  std::cout<<"nb initial vertices = "<<tr.number_of_vertices()<<std::endl;
  std::cout<<"nb initial facet = "<<tr.number_of_facets()<<std::endl;
  std::cout<<"nb initial cells = "<<tr.number_of_cells()<<std::endl;

//------------------------------------------------------------------------
// Meshing
  CGAL::refine_mesh_3(c3t3, domain, criteria, CGAL::parameters::no_exude(), CGAL::parameters::no_perturb());

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
