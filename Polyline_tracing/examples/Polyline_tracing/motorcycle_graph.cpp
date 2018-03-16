#define CGAL_CHECK_EXPENSIVE

#define CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
#define CGAL_MOTORCYCLE_GRAPH_VERBOSE
#define CGAL_MOTORCYCLE_GRAPH_OUTPUT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyline_tracing/Motorcycle_graph.h>
#include <CGAL/Polyline_tracing/Motorcycle_graph_traits_2.h>
#include <CGAL/Polyline_tracing/Point_set_tracer.h>
#include <CGAL/Polyline_tracing/Uniform_direction_tracer_visitor.h>

#include <CGAL/Origin.h>
#include <CGAL/point_generators_2.h>

#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>

#include <fstream>
#include <iostream>
#include <vector>
#include <utility>

namespace CP = CGAL::parameters;
namespace PL = CGAL::Polyline_tracing;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef CGAL::Surface_mesh<K::Point_2>                           PolygonMesh;
typedef PL::Motorcycle_graph_traits_2<K, PolygonMesh>            MGT;

typedef MGT::FT                                                  FT;
typedef MGT::Point_d                                             Point_2;
typedef MGT::Vector_d                                            Vector_2;
typedef MGT::Triangle_d                                          Triangle_2;
typedef MGT::Face_location                                       Face_location;

typedef PL::Motorcycle_graph<MGT>                                Motorcycle_graph;

typedef Motorcycle_graph::Motorcycle                             Motorcycle;
typedef boost::shared_ptr<Motorcycle>                            Motorcycle_ptr;
typedef std::vector<Motorcycle_ptr>                              Motorcycle_container;

typedef PL::Uniform_direction_tracer_visitor<MGT>                Uniform_tracer;
typedef PL::Motorcycle<MGT, Uniform_tracer>                      Motorcycle_U;
typedef PL::Point_set_tracer<MGT>                                Point_set_tracer;
typedef PL::Motorcycle<MGT, Point_set_tracer>                    Motorcycle_PS;

typedef boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;
typedef boost::graph_traits<PolygonMesh>::face_descriptor        face_descriptor;

void motorcycle_club_1(Motorcycle_container& motorcycles, const PolygonMesh& mesh)
{
  // This is a casual motorcycle club that likes to use all of its available options.
  // should be used used with 'eight_triangles.off'

  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::speed = 1.,
                                                        CP::source = Point_2(0.1, 0.1),
                                                        CP::direction = Vector_2(1., 0.),
                                                        CP::initial_time = 0.)));

  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = Point_2(0.9, 0.9),
                                                         CP::direction = Vector_2(0., -1.),
                                                         CP::speed = 1.)));

  face_descriptor fd = *(faces(mesh).begin());
  Face_location loc = std::make_pair(fd, CGAL::make_array(0.4, 0.4, 0.2));
  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = loc,
                                                         CP::direction = Vector_2(1., 1.))));

  Uniform_tracer uft;
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.5, 0.2),
                                                        CP::destination = Point_2(0.95, 0.7),
                                                        CP::tracer = uft)));

  // need a tree or something ! @todo
//  std::vector<Face_location> destinations;
//  destinations.push_back(PMP::locate(Point_2(0.3, 0.6), mesh));
//  Point_set_tracer pst(destinations);
//  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_PS(CP::source = Point_2(0.4, 0.6),
//                                                         CP::tracer = pst)));
}

void motorcycle_club_2(Motorcycle_container& motorcycles)
{
  // This is a motorcycle club with nasty positions and nasty intersections.
  // should be used used with 'eight_triangles.off'

  // The next two should ram into each other (same supporting line, opposite directions)
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(CGAL_PI/15., CGAL_PI/31.),
                                                        CP::destination = Point_2(0.5, 0.2))));
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(1. - CGAL_PI/15., 0.4 - CGAL_PI/31.),
                                                        CP::destination = Point_2(CGAL_PI/15., CGAL_PI/31.))));

  // This motorcycle should crash at the source of motorcycle #1
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(CGAL_PI/30., CGAL_PI/62.),
                                                        CP::destination = Point_2(CGAL_PI/7.5, CGAL_PI/15.5))));

  // The next motorcycle starts at the same point as motorcycle #2, but in another direction
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(1. - CGAL_PI/15., 0.4 - CGAL_PI/31.),
                                                        CP::direction = Vector_2(0.5, -0.2))));

  // The following do NOT have the same supporting lines, but impact each other at the same time
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(CGAL::sqrt(3.)/5., CGAL::sqrt(5.)/5.),
                                                        CP::direction = Vector_2(1./3., 1./3.))));
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(1. - CGAL::sqrt(3.)/5., CGAL::sqrt(5.)/5.),
                                                        CP::direction = Vector_2(-1./3., 1./3.))));

  // Intersects the same point as the last two, but at a later time
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.5, 0.99),
                                                        CP::direction = Vector_2(0., -1.))));

  // The next motorcycles are collinear and the second rams into the first's source
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.9, 0.3),
                                                        CP::destination = Point_2(0.95, 0.6))));
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.9+0.1/3., 0.5),
                                                        CP::direction = Vector_2(1., 6.))));

  // The following motorcycles move in the same direction and from the same point
  // but for numerical reasons they don't see it...
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.6, 0.02),
                                                        CP::destination = Point_2(0.6 + 1./4., 0.02 - 1./100.))));
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.6, 0.02),
                                                        CP::direction = Vector_2(1., -0.04))));

  // The following motorcycles move in the same direction and from the same point,
  // but for numerical reasons they don't see it...
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.1, 0.02),
                                                        CP::destination = Point_2(0.1 + 1./3., 0.02 - 1./97.))));
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.1, 0.02),
                                                        CP::direction = Vector_2(1., -3./97.))));

  // The following motorcycles intersect at an edge
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.6, 0.4),
                                                        CP::direction = Vector_2(-1., 1.))));
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.5, 0.6),
                                                        CP::direction = Vector_2(0., -1.))));
}

void motorcycle_club_3(Motorcycle_container& motorcycles)
{
  // This motorcycle club is all about starting from weird locations.
  // should be used used with 'eight_triangles.off'

  FT eps = std::numeric_limits<FT>::epsilon();

  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0., 0.),
                                                        CP::direction = Vector_2(1., 0.5))));

  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(1., 1./3.),
                                                        CP::direction = Vector_2(0., 1.))));

  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(1., 1./4.),
                                                        CP::direction = Vector_2(1., 1.))));

  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(1., 1./5.),
                                                        CP::direction = Vector_2(eps, 1.))));

  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(1., 1./6.),
                                                        CP::direction = Vector_2(-10 * eps, 1.))));

  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(1., 0.),
                                                        CP::direction = Vector_2(-1, 1.))));
}

void motorcycle_club_4(Motorcycle_container& motorcycles)
{
  // Some configuration that is nastier than it looks
  // should be used used with 'triangle.off'

  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0., 0.1),
                                                        CP::destination = Point_2(0.4, 4.95))));

  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0., 4.95),
                                                        CP::destination = Point_2(0.5, 4.95))));

  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.4, 0.08),
                                                        CP::direction = Vector_2(0, 1.))));

  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.25, 0.2),
                                                        CP::destination = Point_2(0.0, 0.2))));

  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.5, 0.15),
                                                        CP::destination = Point_2(0.3, 0.15))));

//  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(1., 0.),
//                                                        CP::destination = Point_2(-1, 1.))));
}

void motorcycle_club_5(Motorcycle_container& motorcycles)
{
  // This motorcycle club is all about walking the edge(s).
  // should be used used with 'eight_triangles.off'

  face_descriptor fd0 = PolygonMesh::Face_index(0);
  face_descriptor fd1 = PolygonMesh::Face_index(1);

  // #0 Motorcycle walking an edge
  Face_location source_loc = std::make_pair(fd0, CGAL::make_array(0.6, 0.4, 0.));
  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = source_loc,
                                                         CP::direction = Vector_2(1., 1.))));

  // #1 motorcycle walking the same edge as #0 with #1.target = #0.source
  source_loc = std::make_pair(fd0, CGAL::make_array(0.5, 0.5, 0.));
  Face_location destination_loc = std::make_pair(fd0, CGAL::make_array(0.6, 0.4, 0.));
  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = source_loc,
                                                         CP::destination = destination_loc)));

  // #2 starting from inside a face and intersecting on an edge a track on another face
  source_loc = std::make_pair(fd1, CGAL::make_array(0.3, 0.3, 0.4));
  destination_loc = std::make_pair(fd1, CGAL::make_array(0., 0.45, 0.55));
  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = source_loc,
                                                         CP::destination = destination_loc)));

  // #3 starting from inside a face and intersecting on an edge a source on another face
  source_loc = std::make_pair(fd1, CGAL::make_array(0.3, 0.3, 0.4));
  destination_loc = std::make_pair(fd1, CGAL::make_array(0., 0.5, 0.5));
  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = source_loc,
                                                         CP::destination = destination_loc)));

  // #4 starting from inside a face and intersecting on an edge a destination on another face
  destination_loc = std::make_pair(fd1, CGAL::make_array(0., 0.4, 0.6));
  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = source_loc,
                                                         CP::destination = destination_loc)));

  // #5 starting from inside a face and intersecting on an vertex a track on another face
  destination_loc = std::make_pair(fd1, CGAL::make_array(0., 0., 1.));
  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = source_loc,
                                                         CP::destination = destination_loc)));

  // #6 walking an edge and intersecting motorcycle #0 at the vertex [0,0]
  face_descriptor fd7 = PolygonMesh::Face_index(7);
  source_loc = std::make_pair(fd7, CGAL::make_array(0., 0.6, 0.4));
  destination_loc = std::make_pair(fd7, CGAL::make_array(0., 1., 0.));
  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = source_loc,
                                                         CP::destination = destination_loc)));

  // #7 and #8 walk the same edge in the same direction but from different faces
  face_descriptor fd4 = PolygonMesh::Face_index(4);
  source_loc = std::make_pair(fd4, CGAL::make_array(1., 0., 0.));
  destination_loc = std::make_pair(fd4, CGAL::make_array(0.9, 0.1, 0.));
  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = source_loc,
                                                         CP::destination = destination_loc)));

  // #8 and #7 walk the same edge in the same direction but from different faces
  face_descriptor fd5 = PolygonMesh::Face_index(5);
  source_loc = std::make_pair(fd5, CGAL::make_array(0., 0.6, 0.4));
  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = source_loc,
                                                         CP::direction = Vector_2(1., 1.))));

  // #9 starts from the corner [1,-1] and aimes at [0,0]
  face_descriptor fd2 = PolygonMesh::Face_index(2);
  source_loc = std::make_pair(fd2, CGAL::make_array(0., 0., 1.));
  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = source_loc,
                                                         CP::direction = Vector_2(-1., 1.))));

  // #10 and #11 walk in the same direction, on opposite sides of an edge
  source_loc = std::make_pair(fd0, CGAL::make_array(0.1, 0.9, 0.));
  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = source_loc,
                                                         CP::direction = Vector_2(1., 1.))));

  // #11 and #10 walk in the same direction, on opposite sides of an edge
  source_loc = std::make_pair(fd1, CGAL::make_array(0., 0.75, 0.25));
  destination_loc = std::make_pair(fd1, CGAL::make_array(0., 0.8, 0.2));
  motorcycles.push_back(Motorcycle_ptr(new  Motorcycle_U(CP::source = source_loc,
                                                         CP::destination = destination_loc)));
}

void motorcycle_club_6(Motorcycle_container& motorcycles)
{
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.2, 0.18),
                                                        CP::destination = Point_2(0.3, 0.18))));
  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_2(0.6, 0.18),
                                                        CP::destination = Point_2(0.2, 0.18))));
}

void random_motorcycles_in_triangle(Motorcycle_container& motorcycles,
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
    motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = s1, CP::destination = d1)));
    motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = s2, CP::direction = di2)));
  }
}

void random_motorcycles_on_segment(Motorcycle_container& motorcycles,
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
    motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = s1, CP::destination = d1)));
    motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = s2, CP::destination = d2)));
  }
}

void random_motorcycles_in_square(Motorcycle_container& motorcycles,
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
    motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = s1, CP::destination = d1)));
    motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = s2, CP::direction = di2)));
  }
}

void random_motorcycles_on_face(Motorcycle_container& motorcycles,
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
  std::cout.precision(18);
  std::cerr.precision(18);

#if 1//def CGAL_MOTORCYCLE_GRAPH_USE_FIXED_SEEDS
  CGAL::Random rnd(1516813504);
#else
  CGAL::Random rnd(CGAL::get_default_random());
#endif
  std::ofstream seed_out("results_2/seed.txt");

  // read input mesh
  PolygonMesh pm;
  std::ifstream in("data/eight_triangles.off");
  in >> pm;
  CGAL_precondition(pm.number_of_faces() && pm.is_valid());

  bool is_loop_infinite = true;
  while(is_loop_infinite)
  {
    is_loop_infinite = true;

    seed_out << rnd.get_seed() << std::endl;

    Motorcycle_container motorcycles;

    motorcycle_club_1(motorcycles, pm);
    motorcycle_club_2(motorcycles);
    motorcycle_club_3(motorcycles);
//    motorcycle_club_4(motorcycles);
    motorcycle_club_5(motorcycles);
    motorcycle_club_6(motorcycles);

    random_motorcycles_on_face(motorcycles, pm, *(faces(pm).begin()), rnd);
    random_motorcycles_on_face(motorcycles, pm, *(++(++(++(++(++faces(pm).begin()))))), rnd);
    random_motorcycles_on_face(motorcycles, pm, *(++(++faces(pm).begin())), rnd);

    Motorcycle_graph motorcycle_graph(pm);
    motorcycle_graph.trace_graph(motorcycles.begin(), motorcycles.end());

#ifdef CGAL_MOTORCYCLE_GRAPH_OUTPUT
    motorcycle_graph.output_all_dictionary_points();
    for(std::size_t i=0; i<motorcycles.size(); ++i)
    {
//    motorcycle_graph.motorcycle(i).output_intended_track();
      motorcycle_graph.motorcycle(i).output_track();
    }
#endif

    CGAL_postcondition(motorcycle_graph.is_valid());

    rnd = CGAL::Random();
  }

  return 0;
}
