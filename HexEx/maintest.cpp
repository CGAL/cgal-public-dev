#define BOOST_PARAMETER_MAX_ARITY 12
#include <CGAL/Combinatorial_map.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Vector_3.h>
#include <unordered_map>
#include <map>
#include <cmath>
#include <chrono>
#include <vector>
#include <string>
//#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Point_3.h>
#include <cassert>
#include <CGAL/Bbox_3.h>
//#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <CGAL/Triangulation_3_to_lcc.h>
#include<CGAL/Linear_cell_complex_constructors.h>

//#define USE_MESH 1 // comment to use triangulation_3 insteat of Mesh_3
//#define DEBUG 1

#ifdef USE_MESH

// Mesh_3
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#else
// Triangulation_3
#include <CGAL/Delaunay_triangulation_3.h>

#endif

typedef struct dartinfo{

//all the darts corresponding to the same tet will have the same cell_no. This is to make it easier to assign a parametrization matrix, which can be determined from any dart of the tet.
  double cell_no;

//True the dart belongs to a flipped tet.
  bool flipped;

//a variable to store if the particular dart is singular.
  bool singular;

//Number of darts on the same edge.
  int singular_edges;

//parameters is a point corresponding to the point of gven dart in input mesh in the parametrized space 
  CGAL::Point_3<CGAL::Exact_predicates_inexact_constructions_kernel> parameters;
}dart_info;


struct myItem
{
  template < class CMap >
  struct Dart_wrapper
  {
    typedef dart_info Dart_info; 

    typedef CGAL::Cell_attribute_with_point< CMap, int, CGAL::Tag_true>   Vertex_attribute;
    
    typedef CGAL::cpp11::tuple<Vertex_attribute> Attributes;
  };
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel                    K; 
typedef CGAL::Linear_cell_complex_traits<3, K>	                               Traits;
//typedef CGAL::Linear_cell_complex_for_generalized_map<3, 3, Traits, myItem>    LCC_3;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3, 3, Traits, myItem>    LCC_3;

typedef LCC_3::Dart_handle                                                     Dart_handle;
typedef LCC_3::Dart_const_handle                                               Dart_const_handle;
typedef CGAL::Point_3<K>   Point_3;
typedef LCC_3::Point                                       Point;

#ifdef USE_MESH
// for Mesh
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
using namespace CGAL::parameters;

struct Cell_in_complex
{
  Cell_in_complex(const C3t3& ac3t3) : c3t3(ac3t3)
  {}

  bool operator() (C3t3::Cell_handle c)
  { return c3t3.is_in_complex(c); }
  
protected:
  C3t3 c3t3;
};

#else 

typedef CGAL::Delaunay_triangulation_3<LCC_3::Traits> Triangulation;

#endif


void dummy_parametrize(LCC_3& lcc){ /**dummy parametrization function: There is no guarantee that this could give good results, but for the sake of tests this might just work. 
* Every dart has an info as dart_info structure defined in typedefs.h. This contains information about the tet it belongs to (cell_no), and the parametric coordinates of the point. 
* This function updates the parametric coordinates. Any other relevant function can be used by changing the line marked with a comment below:*/
//if(DEBUG) std::cout<<"Inside dummy_parameterize"<<std::endl;
  for(LCC_3::Dart_range::iterator it = (lcc.darts()).begin(), itend = (lcc.darts()).end(); it != itend; it++){
    //if(DEBUG) std::cout<<"Before creating info"<<std::endl;
    dart_info temp;
    //if(DEBUG) std::cout<<"After creating info"<<std::endl;
    temp.cell_no = 0;
  //  if(DEBUG) std::cout<<"Accessed info attribute"<<std::endl;
    Point_3 point(round(10*(lcc.point(it))[0]), round(10*(lcc.point(it))[1]), round(10*(lcc.point(it))[2]));// This line can be changed for trying different parametrizations.
//    if(DEBUG) std::cout<<"Created a point"<<std::endl;
    temp.parameters = point;
    //  if(DEBUG) std::cout<<"Accessed parameters"<<std::endl;
    temp.singular = false;
    //if(DEBUG) std::cout<<"Accessed singular"<<std::endl;
    temp.singular_edges = 0;
//    if(DEBUG) std::cout<<"Before accessing info"<<std::endl;
    lcc.info(it) = temp;
//std::cout<<lcc.point(it)<<"    "<<(lcc.info(it)).parameters<<std::endl;
  //  if(DEBUG) std::cout<<"After accessing info"<<std::endl;
  }
}




bool load_off_to_LCC(std::string filename, LCC_3& lcc)
{
  const char* fname = filename.c_str();
#ifdef USE_MESH
  // For Mesh_3
  double fa = 25, fs = 0.15, fd = 0.008, crer = 3;
  Polyhedron polyhedron;
  std::ifstream input(fname);
  input >> polyhedron;
  if(input.fail()){
    std::cerr << "Error: Cannot read file " <<fname<< std::endl;
    return false;
  }
  input.close();
  Mesh_domain domain(polyhedron);

  Mesh_criteria criteria(fa, fs, fd, crer);
  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

 // Set tetrahedron size (keep cell_radius_edge_ratio), ignore facets
  Mesh_criteria new_criteria(crer, 1.6); //originally cell_size = 0.3

  // Mesh refinement
  CGAL::refine_mesh_3(c3t3, domain, new_criteria);

  // To convert to lcc
  Cell_in_complex cic(c3t3);  
  lcc.clear();
  C3t3::Triangulation &atr= c3t3.triangulation();
  CGAL::import_from_triangulation_3(lcc, atr, cic);
#else
  std::ifstream iFile(filename.c_str());
  if (!iFile)
  {
    std::cout << "Problem reading file " << filename << std::endl;
    return EXIT_FAILURE;
  }
  std::istream_iterator<Point> begin(iFile), end;
  Triangulation T;
  T.insert(begin, end);
  CGAL_assertion(T.is_valid(false));
  CGAL::import_from_triangulation_3(lcc, T);
#endif

  lcc.display_characteristics(std::cout);
  std::cout<<"; Number of 0-attributes: "<<lcc.attributes<0>().size()<<std::endl;
  
  return lcc.is_valid();
}

int main(int argc, char** argv)
{ std::string str;
  if (argc==1)
  {
    std::cout<<"Enter filename: (eg: data/elephant.off)"<<std::endl;
    std::cin>>str; 
  }
  else
  {
    str=argv[1];
  }

  LCC_3 lcc;
  bool res = load_off_to_LCC(str, lcc); //tetmesh to lcc
  std::cout<<"Result of load_off_to_LCC: "<<(res?"OK":"FAILED")<<std::endl;
  
  /*Dart_handle d1 = lcc.make_tetrahedron(Point(-1, 0, 0), Point(0, 2, 0), 
                                        Point(1, 0, 0), Point(1, 1, 2));
  Dart_handle d2 = lcc.make_tetrahedron(Point(0, 2, -1),
                                        Point(-1, 0, -1),
                                        Point(1, 0, -1),
                                        Point(1, 1, -3));
              
                          lcc.sew<3>(d1, d2); */
  dummy_parametrize(lcc);
  std::ofstream ofile;
  ofile.open("triangulation");
  ofile<<lcc; //works
  ofile.close();

  std::ifstream in;
  lcc.clear();
  in.open("triangulation");
  if (in.is_open())
  {
    in>>lcc; //doesn't work: segfault - figure out why
    in.close();
  }
  std::ofstream outfile;
  outfile.open("test.off");
  CGAL::write_off(lcc, outfile);
  outfile.close();
  lcc.display_characteristics(std::cout);
  
  return EXIT_SUCCESS;
}
