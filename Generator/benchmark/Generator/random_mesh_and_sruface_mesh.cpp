#include <iostream>
#include <CGAL/Timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_3/implicit_to_labeled_function_wrapper.h>
#include <CGAL/Mesh_3/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

using namespace std;
using namespace CGAL::parameters;

template <int Sq_radius>
double sphere_function (double x, double y, double z)
{
  double x2=x*x, y2=y*y, z2=z*z;
  return (x2+y2+z2)/Sq_radius - 1;
}

template <typename FT, typename Point>
class FT_to_point_function_wrapper : public std::unary_function<Point, FT>
{
  typedef FT (*Implicit_function)(FT, FT, FT);
  Implicit_function function;
public:
  FT_to_point_function_wrapper(Implicit_function f) : function(f) {}
  FT operator()(Point p) const { return function(p.x(), p.y(), p.z()); }
};

void benchmark_mesh_3()
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef FT_to_point_function_wrapper<K::FT, K::Point_3> Function;
	typedef CGAL::Mesh_3::Implicit_vector_to_labeled_function_wrapper<Function, K>
								Function_wrapper;
	typedef Function_wrapper::Function_vector Function_vector;
	typedef CGAL::Mesh_3::Labeled_mesh_domain_3<Function_wrapper, K> Mesh_domain;
	typedef K::Point_3 Point_3;
	typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
	typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
	typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
	typedef Mesh_criteria::Facet_criteria    Facet_criteria;
	typedef Mesh_criteria::Cell_criteria     Cell_criteria;
	typedef Tr::Geom_traits GT;
	typedef GT::Tetrahedron_3 Tetrahedron3;

	CGAL::Timer t;
	t.start();
	Function f2(&sphere_function<5>);

	Function_vector v;
	v.push_back(&f2);

	Mesh_domain domain(v, K::Sphere_3(CGAL::ORIGIN, 5.*5.), 1e-6);
	Facet_criteria facet_criteria(30, 0.2, 0.02);
	Cell_criteria cell_criteria(2., 0.4);
	Mesh_criteria criteria(facet_criteria, cell_criteria);
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());
	CGAL::perturb_mesh_3(c3t3, domain, time_limit = 10);
	CGAL::exude_mesh_3(c3t3,12);
	CGAL::Random rand;

	t.stop();
	std::cout << "Time elapsed for building the mesh_3: " << t.time() << std::endl;
	t.reset();
	
	t.start();
	int nr = 100000;
	std::vector<Point_3> points;
	points.reserve(nr);
	CGAL::Random_points_in_mesh_3<Point_3, C3t3> g(c3t3);
	CGAL::cpp11::copy_n( g, nr, std::back_inserter(points));
	t.stop();
	std::cout << "Time elapsed for generating the points in mesh_3: " << t.time() << std::endl;
	t.reset();

	points.clear();
	v.clear();
}

typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Triangle_3 Triangle_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
typedef typename C2t3::Vertex_handle Vertex;
FT sphere_function2 (Point_3 p) {
  const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  return x2+y2+z2-1;
}

void benchmark_surface_mesh_3() {
	CGAL::Timer t;
	t.start();
	Tr tr;
	C2t3 c2t3 (tr);

	Surface_3 surface(sphere_function2,
	      	    Sphere_3(CGAL::ORIGIN, 2.));

	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30., 0.1, 0.1);
	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
	CGAL::Random rand;
	
	t.stop();
	std::cout << "Time elapsed for building the mesh: " << t.time() << std::endl;
	t.reset();

	t.start();
	int nr = 100000;
	std::vector<Point_3> points;
	points.reserve(nr);
	CGAL::Random_points_on_surface_mesh_3<Point_3, C2t3> g(c2t3);
	CGAL::cpp11::copy_n( g, nr, std::back_inserter(points));
	t.stop();
	std::cout << "Time elapsed for generating the points in surface_mesh_3: " << t.time() << std::endl;
	t.reset();
	points.clear();
}

int main() {
	std::cout << "Benchmark for mesh_3:" <<std::endl;
	benchmark_mesh_3();
	std::cout << "Benchmark for surface_mesh_3:" <<std::endl;
	benchmark_surface_mesh_3();
	return 0;
}
