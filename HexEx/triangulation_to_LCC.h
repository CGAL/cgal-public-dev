#pragma once
#define BOOST_PARAMETER_MAX_ARITY 12
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include<CGAL/Triangulation_3_to_lcc.h>
// IO
#include <CGAL/IO/Polyhedron_iostream.h>
#include<cstdlib>
#include<fstream>
#include<CGAL/Triangulation_3.h>
#include<CGAL/Linear_cell_complex_constructors.h> 
#include<new>
//#define BOOST_PROPERTY_TREE_RAPIDXML_STATIC_POOL_SIZE 512

/* If you want to use a viewer, you can use qglviewer. */
#ifdef CGAL_LCC_USE_QT
#include "linear_cell_complex_3_viewer_qt.h"
#endif

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

//LCC
//typedef CGAL::Linear_cell_complex_for_generalized_map<3,3> LCC_3;

// To avoid verbose function and named parameters call
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


bool load_off_to_LCC(std::string filename, LCC_3& lcc)
{
  const char* fname = filename.c_str();
  double fa = 30, fs = 0.15, fd = 0.008, crer = 3; 

  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input(fname);
  input >> polyhedron;
  if(input.fail()){
    std::cerr << "Error: Cannot read file " <<fname<< std::endl;
    return false;
  }
  input.close();
  Mesh_domain domain(polyhedron);
  
  

  Mesh_criteria criteria(facet_angle=fa, facet_size=fs, facet_distance=fd,
                         cell_radius_edge_ratio=crer);
  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

 // Set tetrahedron size (keep cell_radius_edge_ratio), ignore facets
  Mesh_criteria new_criteria(crer, cell_size=0.3); //originally cell_size = 0.3

  // Mesh refinement
  CGAL::refine_mesh_3(c3t3, domain, new_criteria);


  // To convert to lcc
  Cell_in_complex cic(c3t3);  
  lcc.clear();
  C3t3::Triangulation &atr= c3t3.triangulation();
  CGAL::import_from_triangulation_3(lcc, atr, cic);
  std::ofstream ofile;
  ofile.open("tri.off");
  CGAL::write_off(lcc, ofile);
  ofile.close();
  //std::ofstream ofile;
/*  ofile.open("triangulation");
lcc.display_characteristics(std::cout);
  ofile<<lcc; //works
 //lcc.display_characteristics(std::cout);
  ofile.close();*/
  /*std::ifstream in;
  lcc.clear();
  in.open("triangulation");
 //lcc.display_characteristics(std::cout);
  /*if (in.is_open())
  {    
    in>>lcc;
    in.close();
  }*/
  lcc.display_characteristics(std::cout);
  std::cout<<std::endl;


#ifdef CGAL_LCC_USE_VIEWER
  display_lcc(lcc);
#endif // CGAL_LCC_USE_VIEWER

  return true;
}
