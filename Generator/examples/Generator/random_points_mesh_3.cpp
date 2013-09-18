#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_3/implicit_to_labeled_function_wrapper.h>
#include <CGAL/Mesh_3/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/point_generators_3.h>

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


void points_in_mesh_3(int N)
{
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
	
	std::vector<Point_3> points;
	points.reserve(N);
	CGAL::Random_points_in_mesh_3<Point_3, C3t3> g(c3t3);
	CGAL::cpp11::copy_n( g, N, std::back_inserter(points));
	std::cout << "First of the "<<N<<" generated points is: "<<std::endl;
	std::cout <<points[0].x()<<" "<<points[0].y()<<" "<<points[0].z()<<std::endl;
	points.clear();
	v.clear();
}

int main() {
	int N;
	std::cout << "Number of points to be generated inside mesh:"<<std::endl;
	std::cin >> N;
	points_in_mesh_3(N);
	return 0;
}
