//#include <CGAL/Implicit_mesh_domain_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_3/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>


#include <CGAL/Mesh_criteria_3_with_balls.h>

template <typename K>
struct Function{
  typedef int return_type; 
  int operator()(typename K::Point_3 p, bool b = false) const {
    if(p.x()*p.x()+p.y()*p.y()+p.z()*p.z() > 1) return 0;
    if(p.x()>0) return 1;
    else return 2;
  }
};

// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef K::FT FT;
typedef K::Point_3 Point;
typedef CGAL::Mesh_3::Labeled_mesh_domain_3<Function<K>,K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3_with_balls<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;


int main()
{
  // Domain (Warning: Sphere_3 constructor uses square radius !)
  Mesh_domain domain(Function<K>(), K::Sphere_3(CGAL::ORIGIN, 2.));

  // Criteria
  Facet_criteria facet_criteria(30, 0.5, 0); // angle, size, approximation
  Cell_criteria cell_criteria(3, 0); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Mesh generation
  C3t3 c3t3; 
  C3t3::Triangulation& tr = c3t3.triangulation();
  typedef C3t3::Triangulation::Point Weighted_point;     
  typedef K::Point_3 Point_3;
  tr.insert(Weighted_point(Point_3(0.001,-0.01,0.01),0));
  tr.insert(Weighted_point(Point_3(1,0,0),0));
  tr.insert(Weighted_point(Point_3(0,1,0),0));
  tr.insert(Weighted_point(Point_3(0,0,1),0));     
  tr.insert(Weighted_point(Point_3(-1,0,0),0));
  tr.insert(Weighted_point(Point_3(0,-1,0),0));
  tr.insert(Weighted_point(Point_3(0,0,-1),0));
  CGAL::refine_mesh_3(c3t3, domain, criteria, CGAL::parameters::no_exude(), CGAL::parameters::no_perturb());

  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file,true,true);

  return 0;
}

