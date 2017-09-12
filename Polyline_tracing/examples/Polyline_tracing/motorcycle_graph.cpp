#define CGAL_MOTORCYCLE_GRAPH_VERBOSE

//#define CGAL_MOTORCYCLE_GRAPH_USE_MANUAL_POINTS
//#define CGAL_MOTORCYCLE_GRAPH_RANDOM_POINTS_ON_SEGMENT
//#define CGAL_MOTORCYCLE_GRAPH_RANDOM_POINTS_IN_SQUARE
#define CGAL_MOTORCYCLE_GRAPH_RANDOM_POINTS_IN_TRIANGLE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyline_tracing/Motorcycle.h>
#include <CGAL/Polyline_tracing/Motorcycle_graph.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h> // for 'copy_n_unique()'

#include <fstream>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;

typedef K::Point_2                                               Point_2;
typedef K::Vector_2                                              Vector_2;

typedef CGAL::Surface_mesh<Point_2>                              PolygonMesh;
typedef CGAL::Polyline_tracing::Motorcycle<K, PolygonMesh>       Motorcycle;
typedef CGAL::Polyline_tracing::Motorcycle_graph<K, PolygonMesh> Motorcycle_graph;

namespace CP = CGAL::parameters;

int main()
{
  std::cout.precision(17);
  std::cerr.precision(17);

  PolygonMesh pm;
  std::ifstream in("data/two_triangles.off");
  in >> pm;
  std::cout << pm.number_of_vertices() << " vertices" << std::endl;
  std::cout << pm.number_of_edges() << " edges" << std::endl;
  std::cout << pm.number_of_faces() << " faces" << std::endl;
  CGAL_precondition(pm.is_valid());

  std::ofstream out("polygon_mesh.off");
  out << pm;

  bool is_loop_infinite = true;
  while(is_loop_infinite)
  {
    is_loop_infinite = false;

    std::vector<Motorcycle> motorcycles;

#ifdef CGAL_MOTORCYCLE_GRAPH_USE_MANUAL_POINTS
    // add some motorcycles
    motorcycles.push_back(Motorcycle(CP::speed = 1.,
                                     CP::source = Point_2(0, 0),
                                     CP::destination = Point_2(1, 0),
                                     CP::initial_time = 0.));
    motorcycles.push_back(Motorcycle(CP::source = Point_2(-1, 0),
                                     CP::direction = Vector_2(0, 1),
                                     CP::speed = 1.));
#else // random stuff below
    const int size = 1; // number of random points
    motorcycles.reserve(size);

 #ifdef CGAL_MOTORCYCLE_GRAPH_RANDOM_POINTS_IN_SQUARE
    typedef CGAL::Random_points_in_square_2<Point_2>            Generator;
    const int side = 1;
  #if CGAL_MOTORCYCLE_GRAPH_USE_FIXED_SEEDS
    CGAL::Random rand_s(6), rand_d(7);
    std::cerr << "Seeds = " << rand_s.get_seed() << " " << rand_d.get_seed() << std::endl;
    Generator gen_s(side, rand_s), gen_d(side, rand_d);
  #else
    Generator gen_s(side), gen_d(side);
  #endif

    for(int i=0; i<size; ++i)
    {
      motorcycles.push_back(Motorcycle(CP::source = *gen_s++,
                                       CP::destination = *gen_d++));
      motorcycles.push_back(Motorcycle(CP::source = *gen_s++,
                                       CP::direction = Vector_2(CGAL::ORIGIN, *gen_d++)));
    }
 #elif defined(CGAL_MOTORCYCLE_GRAPH_RANDOM_POINTS_IN_TRIANGLE)
    typedef CGAL::Random_points_in_triangle_2<Point_2>            Generator;
    const Point_2 t0(0,0), t1(1,0), t2(0, 1);
    Generator gen_s(t0, t1, t2), gen_d(t0, t1, t2);
    for(int i=0; i<size; ++i)
    {
//      motorcycles.push_back(Motorcycle(CP::source = *gen_s++,
//                                       CP::destination = *gen_d++));
      motorcycles.push_back(Motorcycle(CP::source = *gen_s++,
                                       CP::direction = Vector_2(Point_2(1./3.,1./3.), *gen_d++)));
    }
 #elif defined(CGAL_MOTORCYCLE_GRAPH_RANDOM_POINTS_ON_SEGMENT)
    typedef CGAL::Random_points_on_segment_2<Point_2>            Generator;
    const Point_2 s0(0,0), t0(4,0), s1(CGAL_PI, -0.1), t1(-CGAL::sqrt(3.), std::cos(0.1));
    Generator gen_s0(s0, t0), gen_d0(s0, t0);
    Generator gen_s1(s1, t1), gen_d1(s1, t1);

    for(int i=0; i<size; ++i)
    {
      motorcycles.push_back(Motorcycle(CP::source = *gen_s0++,
                                       CP::destination = *gen_d1++));
      motorcycles.push_back(Motorcycle(CP::source = *gen_s1++,
                                       CP::destination = *gen_d1++));
    }
 #endif

#endif // CGAL_MOTORCYCLE_GRAPH_USE_MANUAL_POINTS

    // trace the graph
    Motorcycle_graph motorcycle_graph(pm);
    motorcycle_graph.trace_graph(motorcycles.begin(), motorcycles.end());
    CGAL_postcondition(motorcycle_graph.is_valid());

    // output
    motorcycle_graph.output_motorcycles_sources_and_destinations();
    motorcycle_graph.output_all_dictionary_points();
    for(std::size_t i=0; i<motorcycles.size(); ++i)
    {
  //    motorcycle_graph.motorcycle(i).output_intended_track();
      motorcycle_graph.motorcycle(i).output_track();
    }
  }

  return 0;
}
