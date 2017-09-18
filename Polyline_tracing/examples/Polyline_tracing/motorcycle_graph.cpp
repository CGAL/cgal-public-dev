#define CGAL_MOTORCYCLE_GRAPH_VERBOSE
#define CGAL_MOTORCYCLE_GRAPH_OUTPUT

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

typedef K::FT                                                    FT;
typedef K::Point_2                                               Point_2;
typedef K::Vector_2                                              Vector_2;
typedef K::Triangle_2                                            Triangle_2;

typedef CGAL::Surface_mesh<Point_2>                              PolygonMesh;
typedef CGAL::Polyline_tracing::Motorcycle<K, PolygonMesh>       Motorcycle;
typedef CGAL::Polyline_tracing::Motorcycle_graph<K, PolygonMesh> Motorcycle_graph;

namespace CP = CGAL::parameters;

void motorcycle_club_1(std::vector<Motorcycle>& motorcycles)
{
  // add some motorcycles
  motorcycles.push_back(Motorcycle(CP::speed = 1.,
                                   CP::source = Point_2(0.1, 0.1),
                                   CP::direction = Vector_2(1, 0),
                                   CP::initial_time = 0.));
  motorcycles.push_back(Motorcycle(CP::source = Point_2(0.9, 0.9),
                                   CP::direction = Vector_2(0., -1),
                                   CP::speed = 1.));
}

void random_motorcycles_in_triangle(std::vector<Motorcycle>& motorcycles,
                                    const Triangle_2& triangle = Triangle_2(Point_2(0, 0), Point_2(1, 0), Point_2(0, 1)))
{
  typedef CGAL::Random_points_in_triangle_2<Point_2>            Generator;

  const int size = 10; // number of random points
  motorcycles.reserve(size);

  FT third = 1./3.;
  Point_2 bar = CGAL::barycenter(triangle[0], third, triangle[1], third, triangle[2], third);

  Generator gen_s(triangle), gen_d(triangle);
  for(int i=0; i<size; ++i)
  {
    motorcycles.push_back(Motorcycle(CP::source = *gen_s++,
                                     CP::destination = *gen_d++));
    motorcycles.push_back(Motorcycle(CP::source = *gen_s++,
                                     CP::direction = Vector_2(bar, *gen_d++)));
  }
}

void random_motorcycles_on_segment(std::vector<Motorcycle>& motorcycles)
{
  const int size = 1; // number of random points
  motorcycles.reserve(size);

  typedef CGAL::Random_points_on_segment_2<Point_2>            Generator;
  const Point_2 s0(0,0), t0(4,0), s1(CGAL_PI, -0.1), t1(-CGAL::sqrt(3.), std::cos(0.1));
  Generator gen_s0(s0, t0), gen_d0(s0, t0);
  Generator gen_s1(s1, t1), gen_d1(s1, t1);

  for(int i=0; i<size; ++i)
  {
    motorcycles.push_back(Motorcycle(CP::source = *gen_s0++,
                                     CP::destination = *gen_d0++));
    motorcycles.push_back(Motorcycle(CP::source = *gen_s1++,
                                     CP::destination = *gen_d1++));
  }
}

void random_motorcycles_in_square(std::vector<Motorcycle>& motorcycles)
{
  const int size = 1; // number of random points
  motorcycles.reserve(size);

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
}

void random_motorcycles_on_face(std::vector<Motorcycle>& motorcycles,
                                const PolygonMesh& mesh,
                                const boost::graph_traits<PolygonMesh>::face_descriptor fd)
{
  boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type vpm = get(CGAL::vertex_point, mesh);

  Triangle_2 tr(get(vpm, source(halfedge(fd, mesh), mesh)),
                get(vpm, target(halfedge(fd, mesh), mesh)),
                get(vpm, target(next(halfedge(fd, mesh), mesh), mesh)));
  return random_motorcycles_in_triangle(motorcycles, tr);
}

int main()
{
  std::cout.precision(17);
  std::cerr.precision(17);

  PolygonMesh pm;
  std::ifstream in("data/eight_triangles.off");
  in >> pm;
  CGAL_precondition(pm.is_valid());

//  std::ofstream out("polygon_mesh.off");
//  out << pm;

  bool is_loop_infinite = true;
  while(is_loop_infinite)
  {
    is_loop_infinite = false;

    std::vector<Motorcycle> motorcycles;
//    motorcycle_club_1(motorcycles);
//    random_motorcycles_on_segment(motorcycles);
//    random_motorcycles_in_triangle(motorcycles);
//    random_motorcycles_in_square(motorcycles);
    random_motorcycles_on_face(motorcycles, pm, *(faces(pm).begin()));
    random_motorcycles_on_face(motorcycles, pm, *(++(++(++(++(++faces(pm).begin()))))));

    Motorcycle_graph motorcycle_graph(pm);
    motorcycle_graph.trace_graph(motorcycles.begin(), motorcycles.end());
    CGAL_postcondition(motorcycle_graph.is_valid());

#ifdef CGAL_MOTORCYCLE_GRAPH_OUTPUT
    motorcycle_graph.output_all_dictionary_points();
    for(std::size_t i=0; i<motorcycles.size(); ++i)
    {
//    motorcycle_graph.motorcycle(i).output_intended_track();
      motorcycle_graph.motorcycle(i).output_track();
    }
  }
#endif

  return 0;
}
