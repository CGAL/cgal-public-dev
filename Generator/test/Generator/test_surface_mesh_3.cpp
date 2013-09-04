#include <iostream>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <CGAL/double.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cassert>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

#include <CGAL/internal/Finite_support_distribution.h>
#include <CGAL/internal/Weighted_random_generator.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

using namespace std;

//TODO:typedef CGAL::Exact_predicates_inexact_constructions_kernel 		K;
// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3			 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr>				 C2t3;

typedef Tr::Geom_traits							 GT;
typedef GT::Sphere_3							 Sphere_3;
typedef GT::Triangle_3 							 Triangle_3;
typedef GT::Point_3							 Point;
typedef GT::FT								 FT;

typedef FT (*Function)(Point);
typedef CGAL::Implicit_surface_3<GT, Function>				 Surface_3;
//TODO:typedef K::Plane_3 							Plane_3;
typedef GT::Plane_3 							 Plane_3;

typedef CGAL::FasterMemoryExpensiveTag					 FastPolicy;
typedef CGAL::SlowerMemoryEfficientTag					 SlowPolicy;
	
typedef CGAL::Random_points_in_triangle_3<Point>			 PointGen;
typedef CGAL::internal::Weighted_random_generator<PointGen>		 GeneratorWithWeight;

typedef typename C2t3::Vertex_handle					 Vertex;

FT sphere_function (Point p) {
  const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  return x2+y2+z2-1;
}

const double EPS = 1e-30;

bool inside_or_close_to_surface_mesh_3(Point pt, const Triangle_3 *tri, int size) {
	for (int i = 0; i < size; i++) {
		Plane_3 plane = tri[i].supporting_plane();
		FT dist = squared_distance(plane,pt);
		if(dist<EPS) { 
			return true;
		}
		Triangle_3 OAB = Triangle_3(pt,tri[i][0],tri[i][1]);
		Triangle_3 OAC = Triangle_3(pt,tri[i][0],tri[i][2]);
		Triangle_3 OBC = Triangle_3(pt,tri[i][1],tri[i][2]);
		FT OAB_area = sqrt(OAB.squared_area());
		FT OAC_area = sqrt(OAC.squared_area());
		FT OBC_area = sqrt(OBC.squared_area());
		FT tri_area = sqrt(tri[i].squared_area());
		if(fabs(OAB_area+OAC_area+OBC_area-tri_area)<1e-15) {
			return true;
		}
	}
	return false;
}

template <typename Triangle_3>
class WeightFunctor_triangle_3 {
	public:
		double operator() (Triangle_3 &t) {
			return sqrt(t.squared_area());
		}
};

int main() {
	Tr tr;            // 3D-Delaunay triangulation
	C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

	// defining the surface
	Surface_3 surface(sphere_function,             // pointer to function
	      	    Sphere_3(CGAL::ORIGIN, 2.)); // bounding sphere
	// Note that "2." above is the *squared* radius of the bounding sphere!

	// defining meshing criteria
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
	      					     0.1,  // radius bound
	      					     0.1); // distance bound
	// meshing surface
	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

	CGAL::Random rand;

	for (int number_points = 1; number_points < 10; number_points++) {
		std::vector<Point> points;
		points.reserve(number_points);
		CGAL::Random_points_in_surface_mesh_3<Point, C2t3> g(c2t3);
		Triangle_3 *aux = new Triangle_3[c2t3.number_of_facets()];

		CGAL::cpp11::copy_n(g, number_points, std::back_inserter(points));

		for (int i = 0; i < number_points; i++) {
			typename C2t3::Facet_iterator iter = c2t3.facets_begin();
			int count = 0;
			for (; iter != c2t3.facets_end(); ++iter) {
				Vertex v[4];
				int k = 0;
				for(int j = 0; j < 4; j++)
				{
					if(j == iter->second) continue;
					v[k] = iter->first->vertex(j); // vertices of the facet
					k++;
				}
				aux[count] = Triangle_3(v[0]->point(), v[1]->point(),
						v[2]->point());
				count++;
			}
			assert(inside_or_close_to_surface_mesh_3(points[i],aux,count));
		}
	}

	const int MIN_POINTS = 50;
	const int MAX_POINTS = 100;
	const int number_points = rand.get_int(MIN_POINTS, MAX_POINTS);
	std::vector<Point> points;
	points.reserve(number_points);
	CGAL::Random_points_in_surface_mesh_3<Point, C2t3> g(c2t3);
	Triangle_3 *aux = new Triangle_3[c2t3.number_of_facets()];

	CGAL::cpp11::copy_n(g, number_points, std::back_inserter(points));

	for (int i = 0; i < number_points; i++) {
		typename C2t3::Facet_iterator iter = c2t3.facets_begin();
		int count = 0;
		for (; iter != c2t3.facets_end(); ++iter) {
			Vertex v[4];
			int k = 0;
			for(int j = 0; j < 4; j++)
			{
				if(j == iter->second) continue;
				v[k] = iter->first->vertex(j); // vertices of the facet
				k++;
			}
			aux[count] = Triangle_3(v[0]->point(), v[1]->point(),
					v[2]->point());
			count++;
		}
		assert(inside_or_close_to_surface_mesh_3(points[i],aux,count));
	}

// Testing the copy-constructor
	points.clear();
	points.reserve(number_points);
	CGAL::Random_points_in_surface_mesh_3<Point, C2t3> g1(g);
	delete[] aux;
	aux = new Triangle_3[c2t3.number_of_facets()];

	CGAL::cpp11::copy_n(g1, number_points, std::back_inserter(points));

	for (int i = 0; i < number_points; i++) {
		typename C2t3::Facet_iterator iter = c2t3.facets_begin();
		int count = 0;
		for (; iter != c2t3.facets_end(); ++iter) {
			Vertex v[4];
			int k = 0;
			for(int j = 0; j < 4; j++)
			{
				if(j == iter->second) continue;
				v[k] = iter->first->vertex(j); // vertices of the facet
				k++;
			}
			aux[count] = Triangle_3(v[0]->point(), v[1]->point(),
					v[2]->point());
			count++;
		}
		assert(inside_or_close_to_surface_mesh_3(points[i],aux,count));
	}

// Testing the constructor that has FSD as argument
	points.clear();
	points.reserve(number_points);

	int Nr_facets = c2t3.number_of_facets();
	WeightFunctor_triangle_3<Triangle_3> weightElem;
	GeneratorWithWeight* containing_structure;
	containing_structure = new GeneratorWithWeight[Nr_facets];
	typename C2t3::Facet_iterator iter = c2t3.facets_begin();
	int i = 0;
	for (; iter != c2t3.facets_end(); ++iter) {
		Vertex v[4];
		int k = 0;
		for(int j = 0; j < 4; j++)
		{
			if(j == iter->second) continue;
			v[k] = iter->first->vertex(j); // vertices of the facet
			k++;
		}
		Triangle_3 aux(v[0]->point(), v[1]->point(),
				v[2]->point());
		double weight = weightElem(aux);
		PointGen randGen(aux);
		containing_structure[i] = GeneratorWithWeight (randGen, weight);
		i++;
	}
	int N = 1;
//	int N = 1<<10;
	CGAL::internal::Finite_support_distribution<GeneratorWithWeight> fsd =
		CGAL::internal::Finite_support_distribution<GeneratorWithWeight>
		(containing_structure, i, N);

	CGAL::Random_points_in_surface_mesh_3<Point, C2t3> g2(fsd);
	delete[] aux;
	aux = new Triangle_3[c2t3.number_of_facets()];

	CGAL::cpp11::copy_n(g2, number_points, std::back_inserter(points));

	for (int i = 0; i < number_points; i++) {
		typename C2t3::Facet_iterator iter = c2t3.facets_begin();
		int count = 0;
		for (; iter != c2t3.facets_end(); ++iter) {
			Vertex v[4];
			int k = 0;
			for(int j = 0; j < 4; j++)
			{
				if(j == iter->second) continue;
				v[k] = iter->first->vertex(j); // vertices of the facet
				k++;
			}
			aux[count] = Triangle_3(v[0]->point(), v[1]->point(),
					v[2]->point());
			count++;
		}
		assert(inside_or_close_to_surface_mesh_3(points[i],aux,count));
	}

	points.clear();
	return 0;
}

