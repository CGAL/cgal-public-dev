// #define CGAL_PROFILE
// #define CGAL_NO_ASSERTIONS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <boost/foreach.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_2<K>                            T;

int main(int argc, char *argv[]) {
	T t;
	
	t.insert(T::Point(0,0));
	t.insert(T::Point(1,0));
	t.insert(T::Point(0,1));
	t.insert(T::Point(1,2));

	std::cout << "Iterator range over the points in the triangulation:" << std::endl;
	BOOST_FOREACH(T::Point &p, t.points()) {
		std::cout << p.x() << ", " << p.y() << std::endl;
	}
	std::cout << "Using a non-dereferencing iterator for Finite_vertices_iterator:" << std::endl;
	T::Finite_vertices_range vhs = t.finite_vertices();
	BOOST_FOREACH(T::Finite_vertices_iterator &vh, vhs) {
		std::cout << vh->point().x() << ", " << vh->point().y() << std::endl;
	}
	std::cout << "Implicit cast to a Vertex_handle:" << std::endl;
	BOOST_FOREACH(T::Vertex_handle vh, t.finite_vertices()) {
		std::cout << vh->point().x() << ", " << vh->point().y() << std::endl;
	}

	return 0;
}
