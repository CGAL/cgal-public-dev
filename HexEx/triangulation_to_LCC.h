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
typedef CGAL::Linear_cell_complex_for_generalized_map<3,3> LCC_3;

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


bool load_off_to_LCC(std::string filename, LCC_3& lcc, LCC_3& lcc2)
{  
  std::ifstream ifile1, ifile2;
  ifile1.open(("mesh/"+filename).c_str());
  ifile2.open(("mesh/parametrized"+filename).c_str());
  CGAL::load_off(lcc, ifile1);
  CGAL::load_off(lcc2, ifile2); 
  ifile1.close();
  ifile2.close();

#ifdef CGAL_LCC_USE_VIEWER
  display_lcc(lcc);
#endif // CGAL_LCC_USE_VIEWER

  return true;
}
