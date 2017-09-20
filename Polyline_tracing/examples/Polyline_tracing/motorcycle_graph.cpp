#define CGAL_MOTORCYCLE_GRAPH_VERBOSE
#define CGAL_MOTORCYCLE_GRAPH_OUTPUT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyline_tracing/Motorcycle.h>
#include <CGAL/Polyline_tracing/Motorcycle_graph.h>

#include <CGAL/Origin.h>
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
  // This is a casual motorcycle club that just likes to use all of its options

  motorcycles.push_back(Motorcycle(CP::speed = 1.,
                                   CP::source = Point_2(0.1, 0.1),
                                   CP::direction = Vector_2(1, 0),
                                   CP::initial_time = 0.));
  motorcycles.push_back(Motorcycle(CP::source = Point_2(0.9, 0.9),
                                   CP::direction = Vector_2(0., -1),
                                   CP::speed = 1.));
}

void motorcycle_club_2(std::vector<Motorcycle>& motorcycles)
{
  // This is a motorcycle club with nasty positions and nasty intersections

  // The next two should ram into each other (same supporting line, opposite directions)
  motorcycles.push_back(Motorcycle(CP::source = Point_2(CGAL_PI/15., CGAL_PI/31.),
                                   CP::destination = Point_2(0.5, 0.2)));
  motorcycles.push_back(Motorcycle(CP::source = Point_2(1. - CGAL_PI/15., 0.2 + CGAL_PI/31.),
                                   CP::destination = Point_2(0.5, 0.2)));

  // This motorcycle should crash at the source of motorcycle #1
  motorcycles.push_back(Motorcycle(CP::source = Point_2(CGAL_PI/30., CGAL_PI/62.),
                                   CP::destination = Point_2(CGAL_PI/7.5, CGAL_PI/15.5)));

  // The next motorcycle starts at the same point as motorcycle #2, but in another direction
  motorcycles.push_back(Motorcycle(CP::source = Point_2(1. - CGAL_PI/15., 0.2 + CGAL_PI/31.),
                                   CP::direction = Vector_2(0.5, -0.2)));

  // The following do NOT have the same supporting lines, but impact each other at the same time
  motorcycles.push_back(Motorcycle(CP::source = Point_2(CGAL::sqrt(3.)/5., CGAL::sqrt(5.)/5.),
                                   CP::direction = Vector_2(1./3., 1./3.)));
  motorcycles.push_back(Motorcycle(CP::source = Point_2(1. - CGAL::sqrt(3.)/5., CGAL::sqrt(5.)/5.),
                                   CP::direction = Vector_2(-1./3., 1./3.)));

  // Intersects the same point as the last two, but at a later time
  motorcycles.push_back(Motorcycle(CP::source = Point_2(0.5, 0.99),
                                   CP::direction = Vector_2(0., -1.)));

  // The next motorcycles are collinear and the second rams into the first's source
  motorcycles.push_back(Motorcycle(CP::source = Point_2(0.9, 0.3),
                                   CP::destination = Point_2(0.95, 0.6)));
  motorcycles.push_back(Motorcycle(CP::source = Point_2(0.9+0.1/3., 0.5),
                                   CP::direction = Vector_2(1., 6.)));

  // The following motorcycles move in the same direction and from the same point
  // but for numerical reasons they don't see it...
  motorcycles.push_back(Motorcycle(CP::source = Point_2(0.6, 0.02),
                                   CP::destination = Point_2(0.6 + 1./4., 0.02 - 1./100.)));
  motorcycles.push_back(Motorcycle(CP::source = Point_2(0.6, 0.02),
                                   CP::direction = Vector_2(1., -0.04)));

  // The following motorcycles move in the same direction and from the same point,
  // but for numerical reasons they don't see it...
  motorcycles.push_back(Motorcycle(CP::source = Point_2(0.1, 0.02),
                                   CP::destination = Point_2(0.1 + 1./3., 0.02 - 1./97.)));
  motorcycles.push_back(Motorcycle(CP::source = Point_2(0.1, 0.02),
                                   CP::direction = Vector_2(1., -3./97.)));
}

void motorcycle_club_3(std::vector<Motorcycle>& motorcycles)
{
  // This motorcycle club is all about starting from weird locations

  FT eps = std::numeric_limits<FT>::epsilon();

  motorcycles.push_back(Motorcycle(CP::source = Point_2(0., 0.),
                                   CP::direction = Vector_2(1., 0.5)));

  motorcycles.push_back(Motorcycle(CP::source = Point_2(1., 1./3.),
                                   CP::direction = Vector_2(0., 1.)));

  motorcycles.push_back(Motorcycle(CP::source = Point_2(1., 1./4.),
                                   CP::direction = Vector_2(1., 1.)));

  motorcycles.push_back(Motorcycle(CP::source = Point_2(1., 1./5.),
                                   CP::direction = Vector_2(eps, 1.)));

  motorcycles.push_back(Motorcycle(CP::source = Point_2(1., 1./6.),
                                   CP::direction = Vector_2(-10 * eps, 1.)));

  motorcycles.push_back(Motorcycle(CP::source = Point_2(1., 0.),
                                   CP::direction = Vector_2(-1, 1.)));
}

void random_motorcycles_in_triangle(std::vector<Motorcycle>& motorcycles,
                                    const Triangle_2& triangle,
                                    CGAL::Random& rnd)
{
  const int size = 15; // number of random points
  motorcycles.reserve(size);

  FT third = 1./3.;
  Point_2 bar = CGAL::barycenter(triangle[0], third, triangle[1], third, triangle[2], third);

  CGAL::Random_points_in_triangle_2<Point_2> gen(triangle, rnd);

  std::cout << "seed: " << rnd.get_seed() << std::endl;
  std::cout << "triangle: " << triangle << std::endl;

  for(int i=0; i<size; ++i)
  {
    const Point_2& s1 = *gen++; const Point_2& d1 = *gen++;
    const Point_2& s2 = *gen++; const Vector_2 di2(bar, *gen++);
    motorcycles.push_back(Motorcycle(CP::source = s1, CP::destination = d1));
    motorcycles.push_back(Motorcycle(CP::source = s2, CP::direction = di2));
  }
}

void random_motorcycles_on_segment(std::vector<Motorcycle>& motorcycles,
                                   CGAL::Random& rnd)
{
  const int size = 1; // number of random points
  motorcycles.reserve(size);

  const Point_2 s0(0,0), t0(4,0), s1(CGAL_PI, -0.1), t1(-CGAL::sqrt(3.), std::cos(0.1));
  CGAL::Random_points_on_segment_2<Point_2> gen_s0(s0, t0, rnd), gen_s1(s1, t1, rnd);

  for(int i=0; i<size; ++i)
  {
    const Point_2& s1 = *gen_s0++; const Point_2& d1 = *gen_s0++;
    const Point_2& s2 = *gen_s1++; const Point_2& d2 = *gen_s1++;
    motorcycles.push_back(Motorcycle(CP::source = s1, CP::destination = d1));
    motorcycles.push_back(Motorcycle(CP::source = s2, CP::destination = d2));
  }
}

void random_motorcycles_in_square(std::vector<Motorcycle>& motorcycles,
                                  const FT square_side,
                                  CGAL::Random& rnd)
{
  const int size = 1; // number of random points
  motorcycles.reserve(size);

  CGAL::Random_points_in_square_2<Point_2> gen(square_side, rnd);

  for(int i=0; i<size; ++i)
  {
    const Point_2& s1 = *gen++; const Point_2& d1 = *gen++;
    const Point_2& s2 = *gen++; const Vector_2 di2(CGAL::ORIGIN, *gen++);
    motorcycles.push_back(Motorcycle(CP::source = s1, CP::destination = d1));
    motorcycles.push_back(Motorcycle(CP::source = s2, CP::direction = di2));
  }
}

void random_motorcycles_on_face(std::vector<Motorcycle>& motorcycles,
                                const PolygonMesh& mesh,
                                const boost::graph_traits<PolygonMesh>::face_descriptor fd,
                                CGAL::Random& rnd)
{
  boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type vpm = get(CGAL::vertex_point, mesh);

  Triangle_2 tr(get(vpm, source(halfedge(fd, mesh), mesh)),
                get(vpm, target(halfedge(fd, mesh), mesh)),
                get(vpm, target(next(halfedge(fd, mesh), mesh), mesh)));
  return random_motorcycles_in_triangle(motorcycles, tr, rnd);
}

int main()
{
  std::cout.precision(17);
  std::cerr.precision(17);

#if 1//def CGAL_MOTORCYCLE_GRAPH_USE_FIXED_SEEDS
  CGAL::Random rnd(1505826388);
#else
  CGAL::Random rnd(CGAL::get_default_random());
#endif
  std::ofstream seed_out("seeds.txt");
  seed_out << rnd.get_seed() << std::endl;

  PolygonMesh pm;
  std::ifstream in("data/eight_triangles.off");
  in >> pm;
  CGAL_precondition(pm.is_valid());

//  std::ofstream out("polygon_mesh.off");
//  out << pm;

  bool is_loop_infinite = true;
  while(is_loop_infinite)
  {
    is_loop_infinite = true;

    std::vector<Motorcycle> motorcycles;
//    motorcycle_club_1(motorcycles);
//    motorcycle_club_2(motorcycles);
    motorcycle_club_3(motorcycles);

//    random_motorcycles_on_segment(motorcycles, rnd);
//    random_motorcycles_in_triangle(motorcycles, rnd);
//    random_motorcycles_in_square(motorcycles, rnd);

    random_motorcycles_on_face(motorcycles, pm, *(faces(pm).begin()), rnd);
    random_motorcycles_on_face(motorcycles, pm, *(++(++(++(++(++faces(pm).begin()))))), rnd);
    random_motorcycles_on_face(motorcycles, pm, *(++(++faces(pm).begin())), rnd);

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
