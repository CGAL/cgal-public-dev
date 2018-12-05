#define CGAL_CHECK_EXPENSIVE

#define CGAL_MOTORCYCLE_GRAPH_VERBOSE
#define CGAL_MOTORCYCLE_GRAPH_OUTPUT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/polyline_tracing.h>

#include <CGAL/Origin.h>
#include <CGAL/point_generators_2.h>

#include <fstream>
#include <iostream>
#include <vector>

namespace CP = CGAL::parameters;
namespace PL = CGAL::Polyline_tracing;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef CGAL::Surface_mesh<K::Point_2>                           Triangle_mesh;
typedef PL::Motorcycle_graph_traits_2<K, Triangle_mesh>          MGT;

typedef MGT::FT                                                  FT;
typedef MGT::Point_d                                             Point_2;
typedef MGT::Vector_d                                            Vector_2;
typedef MGT::Triangle_d                                          Triangle_2;
typedef MGT::Face_location                                       Face_location;


typedef PL::Motorcycle<MGT>                                      Motorcycle;
typedef PL::Motorcycle_with_info<MGT, int>                       Motorcycle_with_info;
typedef PL::Motorcycle_graph<MGT, Motorcycle>                    Motorcycle_graph;

typedef PL::Uniform_direction_tracer_visitor<Motorcycle_graph>   Uniform_tracer;
typedef PL::Point_set_tracer<Motorcycle_graph>                   Point_set_tracer;

typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;
typedef boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;

// A temporary motorcycle club out to break some stuff
void motorcycle_club_0(Motorcycle_graph& motorcycle_graph)
{
  motorcycle_graph.add_motorcycle(Point_2(0.1, 0.1) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0., 0.))
                                     .speed(1.));

  motorcycle_graph.add_motorcycle(Point_2(0.4, 0.4),
                                  Uniform_tracer(Vector_2(-1., -1.)),
                                  CP::speed(1.));

  motorcycle_graph.add_motorcycle(Point_2(-0.5, -0.),
                                  Uniform_tracer(Vector_2(1., 0.)),
                                  CP::speed(1.));

  motorcycle_graph.add_motorcycle(Point_2(0.9, -0.9),
                                  Uniform_tracer(Vector_2(-0.9, 0.9)),
                                  CP::speed(1.));

  motorcycle_graph.add_motorcycle(Point_2(-0.034911, 0.9918),
                                  Uniform_tracer(),
                                  CP::destination(Point_2(-0.02, 0.98)));
}

// This is a casual motorcycle club that likes to use all of its available options.
// should be used used with 'eight_triangles.off'
void motorcycle_club_1(Motorcycle_graph& motorcycle_graph)
{
  const Triangle_mesh& mesh = motorcycle_graph.mesh();

  motorcycle_graph.add_motorcycle(Point_2(0.1, 0.1),
                                  Uniform_tracer(Vector_2(1., 0.)),
                                  CP::speed(1.)
                                     .initial_time(0.));

  motorcycle_graph.add_motorcycle(Point_2(0.9, 0.9),
                                  Uniform_tracer(Vector_2(0., -1.)));

  face_descriptor fd = *(faces(mesh).begin());
  Face_location loc = std::make_pair(fd, CGAL::make_array(0.4, 0.4, 0.2));
  motorcycle_graph.add_motorcycle(loc, Uniform_tracer(Vector_2(1., 1.)));

  Uniform_tracer uft;
  motorcycle_graph.add_motorcycle(Point_2(0.5, 0.2),
                                  uft,
                                  CP::destination(Point_2(0.95, 0.7)));

  // need a tree or something ! @todo
//  std::vector<Face_location> destinations;
//  destinations.push_back(PMP::locate(Point_2(0.3, 0.6), mesh));
//  Point_set_tracer pst(destinations);
//  motorcycle_graph.add_motorcycle(Point_2(0.4, 0.6), pst));
}

// This is a motorcycle club with nasty positions and nasty intersections.
// should be used used with 'eight_triangles.off'
void motorcycle_club_2(Motorcycle_graph& motorcycle_graph)
{
  // The next two should ram into each other (same supporting line, opposite directions)
  motorcycle_graph.add_motorcycle(Point_2(CGAL_PI/15., CGAL_PI/31.),
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.5, 0.2)));

  motorcycle_graph.add_motorcycle(Point_2(1. - CGAL_PI/15., 0.4 - CGAL_PI/31.),
                                  Uniform_tracer(),
                                  CP::destination(Point_2(CGAL_PI/15., CGAL_PI/31.)));

  // This motorcycle should crash at the source of motorcycle #1
  motorcycle_graph.add_motorcycle(Point_2(CGAL_PI/30., CGAL_PI/62.),
                                  Uniform_tracer(),
                                  CP::destination(Point_2(CGAL_PI/7.5, CGAL_PI/15.5)));

  // The next motorcycle starts at the same point as motorcycle #2, but in another direction
  motorcycle_graph.add_motorcycle(Point_2(1. - CGAL_PI/15., 0.4 - CGAL_PI/31.),
                                  Uniform_tracer(Vector_2(0.5, -0.2)));

  // The following do NOT have the same supporting lines, but impact each other at the same time
  motorcycle_graph.add_motorcycle(Point_2(CGAL::sqrt(3.)/5., CGAL::sqrt(5.)/5.),
                                  Uniform_tracer(Vector_2(1./3., 1./3.)));

  motorcycle_graph.add_motorcycle(Point_2(1. - CGAL::sqrt(3.)/5., CGAL::sqrt(5.)/5.),
                                  Uniform_tracer(Vector_2(-1./3., 1./3.)));

  // Intersects the same point as the last two, but at a later time
  motorcycle_graph.add_motorcycle(Point_2(0.5, 0.99),
                                  Uniform_tracer(Vector_2(0., -1.)));

  // The next motorcycles are collinear and the second rams into the first's source
  motorcycle_graph.add_motorcycle(Point_2(0.9, 0.3),
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.95, 0.6)));
  motorcycle_graph.add_motorcycle(Point_2(0.9+0.1/3., 0.5),
                                  Uniform_tracer(Vector_2(1., 6.)));

  // The following motorcycles move in the same direction and from the same point
  // but for numerical reasons they don't see it...
  motorcycle_graph.add_motorcycle(Point_2(0.6, 0.02),
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.6 + 1./4., 0.02 - 1./100.)));
  motorcycle_graph.add_motorcycle(Point_2(0.6, 0.02),
                                  Uniform_tracer(Vector_2(1., -0.04)));

  // The following motorcycles move in the same direction and from the same point,
  // but for numerical reasons they don't see it...
  motorcycle_graph.add_motorcycle(Point_2(0.1, 0.02),
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.1 + 1./3., 0.02 - 1./97.)));
  motorcycle_graph.add_motorcycle(Point_2(0.1, 0.02),
                                  Uniform_tracer(Vector_2(1., -3./97.)));

  // The following motorcycles intersect at an edge
  motorcycle_graph.add_motorcycle(Point_2(0.6, 0.4),
                                  Uniform_tracer(Vector_2(-1., 1.)));
  motorcycle_graph.add_motorcycle(Point_2(0.5, 0.6),
                                  Uniform_tracer(Vector_2(0., -1.)));
}

// This motorcycle club is all about starting from weird locations.
// should be used used with 'eight_triangles.off'
void motorcycle_club_3(Motorcycle_graph& motorcycle_graph)
{
  FT eps = std::numeric_limits<FT>::epsilon();

  motorcycle_graph.add_motorcycle(Point_2(0., 0.),
                                  Uniform_tracer(Vector_2(1., 0.5)));

  motorcycle_graph.add_motorcycle(Point_2(1., 1./3.),
                                  Uniform_tracer(Vector_2(0., 1.)));

  motorcycle_graph.add_motorcycle(Point_2(1., 1./4.),
                                  Uniform_tracer(Vector_2(1., 1.)));

  motorcycle_graph.add_motorcycle(Point_2(1., 1./5.),
                                  Uniform_tracer(Vector_2(eps, 1.)));

  motorcycle_graph.add_motorcycle(Point_2(1., 1./6.),
                                  Uniform_tracer(Vector_2(-10 * eps, 1.)));

  motorcycle_graph.add_motorcycle(Point_2(1., 0.),
                                  Uniform_tracer(Vector_2(-1, 1.)));
}

// Some configuration that is nastier than it looks
// should be used used with 'triangle.off'
void motorcycle_club_4(Motorcycle_graph& motorcycle_graph)
{
  motorcycle_graph.add_motorcycle(Point_2(0., 0.1),
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.4, 4.95)));

  motorcycle_graph.add_motorcycle(Point_2(0., 4.95),
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.5, 4.95)));

  motorcycle_graph.add_motorcycle(Point_2(0.4, 0.08),
                                  Uniform_tracer(Vector_2(0, 1.)));

  motorcycle_graph.add_motorcycle(Point_2(0.25, 0.2),
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.0, 0.2)));

  motorcycle_graph.add_motorcycle(Point_2(0.5, 0.15),
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.3, 0.15)));

//  motorcycle_graph.add_motorcycle(Point_2(1., 0.),
//                                  Uniform_tracer(),
//                                  CP::destination(Point_2(-1, 1.)));
}

// This motorcycle club is all about walking the edge(s).
// should be used used with 'eight_triangles.off'
void motorcycle_club_5(Motorcycle_graph& motorcycle_graph)
{
  face_descriptor fd0 = Triangle_mesh::Face_index(0);
  face_descriptor fd1 = Triangle_mesh::Face_index(1);

  // #0 Motorcycle walking an edge
  Face_location source_loc = std::make_pair(fd0, CGAL::make_array(0.6, 0.4, 0.));
  motorcycle_graph.add_motorcycle(source_loc,
                                  Uniform_tracer(Vector_2(1., 1.)));

  // #1 motorcycle walking the same edge as #0 with #1.target = #0.source
  source_loc = std::make_pair(fd0, CGAL::make_array(0.5, 0.5, 0.));
  Face_location destination_loc = std::make_pair(fd0, CGAL::make_array(0.6, 0.4, 0.));
  motorcycle_graph.add_motorcycle(source_loc,
                                  Uniform_tracer(),
                                  CP::destination(destination_loc));

  // #2 starting from inside a face and intersecting on an edge a track on another face
  source_loc = std::make_pair(fd1, CGAL::make_array(0.3, 0.3, 0.4));
  destination_loc = std::make_pair(fd1, CGAL::make_array(0., 0.45, 0.55));
  motorcycle_graph.add_motorcycle(source_loc,
                                  Uniform_tracer(),
                                  CP::destination(destination_loc));

  // #3 starting from inside a face and intersecting on an edge a source on another face
  source_loc = std::make_pair(fd1, CGAL::make_array(0.3, 0.3, 0.4));
  destination_loc = std::make_pair(fd1, CGAL::make_array(0., 0.5, 0.5));
  motorcycle_graph.add_motorcycle(source_loc,
                                  Uniform_tracer(),
                                  CP::destination(destination_loc));

  // #4 starting from inside a face and intersecting on an edge a destination on another face
  destination_loc = std::make_pair(fd1, CGAL::make_array(0., 0.4, 0.6));
  motorcycle_graph.add_motorcycle(source_loc,
                                  Uniform_tracer(),
                                  CP::destination(destination_loc));

  // #5 starting from inside a face and intersecting on an vertex a track on another face
  destination_loc = std::make_pair(fd1, CGAL::make_array(0., 0., 1.));
  motorcycle_graph.add_motorcycle(source_loc,
                                  Uniform_tracer(),
                                  CP::destination(destination_loc));

  // #6 walking an edge and intersecting motorcycle #0 at the vertex [0,0]
  face_descriptor fd7 = Triangle_mesh::Face_index(7);
  source_loc = std::make_pair(fd7, CGAL::make_array(0., 0.6, 0.4));
  destination_loc = std::make_pair(fd7, CGAL::make_array(0., 1., 0.));
  motorcycle_graph.add_motorcycle(source_loc,
                                  Uniform_tracer(),
                                  CP::destination(destination_loc));

  // #7 and #8 walk the same edge in the same direction but from different faces
  face_descriptor fd4 = Triangle_mesh::Face_index(4);
  source_loc = std::make_pair(fd4, CGAL::make_array(1., 0., 0.));
  destination_loc = std::make_pair(fd4, CGAL::make_array(0.9, 0.1, 0.));
  motorcycle_graph.add_motorcycle(source_loc,
                                  Uniform_tracer(),
                                  CP::destination(destination_loc));

  // #8 and #7 walk the same edge in the same direction but from different faces
  face_descriptor fd5 = Triangle_mesh::Face_index(5);
  source_loc = std::make_pair(fd5, CGAL::make_array(0., 0.6, 0.4));
  motorcycle_graph.add_motorcycle(source_loc,
                                  Uniform_tracer(),
                                  Uniform_tracer(Vector_2(1., 1.)));

  // #9 starts from the corner [1,-1] and aimes at [0,0]
  face_descriptor fd2 = Triangle_mesh::Face_index(2);
  source_loc = std::make_pair(fd2, CGAL::make_array(0., 0., 1.));
  motorcycle_graph.add_motorcycle(source_loc,
                                  Uniform_tracer(),
                                  Uniform_tracer(Vector_2(-1., 1.)));

  // #10 and #11 walk in the same direction, on opposite sides of an edge
  source_loc = std::make_pair(fd0, CGAL::make_array(0.1, 0.9, 0.));
  motorcycle_graph.add_motorcycle(source_loc,
                                  Uniform_tracer(),
                                  Uniform_tracer(Vector_2(1., 1.)));

  // #11 and #10 walk in the same direction, on opposite sides of an edge
  source_loc = std::make_pair(fd1, CGAL::make_array(0., 0.75, 0.25));
  destination_loc = std::make_pair(fd1, CGAL::make_array(0., 0.8, 0.2));
  motorcycle_graph.add_motorcycle(source_loc,
                                  Uniform_tracer(),
                                  CP::destination(destination_loc));
}

void motorcycle_club_6(Motorcycle_graph& motorcycle_graph)
{
  motorcycle_graph.add_motorcycle(Point_2(0.2, 0.18),
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.3, 0.18)));
  motorcycle_graph.add_motorcycle(Point_2(0.6, 0.18),
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.2, 0.18)));
}

void random_motorcycles_in_triangle(Motorcycle_graph& motorcycle_graph,
                                    const Triangle_2& triangle,
                                    CGAL::Random& rnd)
{
  const int size = 1; // number of random points

  FT third = 1./3.;
  Point_2 bar = CGAL::barycenter(triangle[0], third, triangle[1], third, triangle[2], third);

  CGAL::Random_points_in_triangle_2<Point_2> gen(triangle, rnd);

  std::cout << "seed: " << rnd.get_seed() << std::endl;
  std::cout << "triangle: " << triangle << std::endl;

  for(int i=0; i<size; ++i)
  {
    const Point_2& s1 = *gen++; const Point_2& d1 = *gen++;
    const Point_2& s2 = *gen++; const Vector_2 di2(bar, *gen++);
    motorcycle_graph.add_motorcycle(s1, Uniform_tracer(), CP::destination(d1));
    motorcycle_graph.add_motorcycle(s2, Uniform_tracer(di2));
  }
}

void random_motorcycles_on_segment(Motorcycle_graph& motorcycle_graph,
                                   CGAL::Random& rnd)
{
  const int size = 50; // number of random points

  const Point_2 s0(0,0), t0(4,0), s1(CGAL_PI, -0.1), t1(-CGAL::sqrt(3.), std::cos(0.1));
  CGAL::Random_points_on_segment_2<Point_2> gen_s0(s0, t0, rnd), gen_s1(s1, t1, rnd);

  for(int i=0; i<size; ++i)
  {
    const Point_2& s1 = *gen_s0++; const Point_2& d1 = *gen_s0++;
    const Point_2& s2 = *gen_s1++; const Point_2& d2 = *gen_s1++;
    motorcycle_graph.add_motorcycle(s1, Uniform_tracer(), CP::destination(d1));
    motorcycle_graph.add_motorcycle(s2, Uniform_tracer(), CP::destination(d2));
  }
}

void random_motorcycles_in_square(Motorcycle_graph& motorcycle_graph,
                                  const FT square_side,
                                  CGAL::Random& rnd)
{
  const int size = 50; // number of random points

  CGAL::Random_points_in_square_2<Point_2> gen(square_side, rnd);

  for(int i=0; i<size; ++i)
  {
    const Point_2& s1 = *gen++;
    const Point_2& d1 = *gen++;
    const Point_2& s2 = *gen++;
    const Vector_2 di2(CGAL::ORIGIN, *gen++);
    motorcycle_graph.add_motorcycle(s1, Uniform_tracer(), CP::destination(d1));
    motorcycle_graph.add_motorcycle(s2, Uniform_tracer(di2));
  }
}

void random_motorcycles_on_face(Motorcycle_graph& motorcycle_graph,
                                const boost::graph_traits<Triangle_mesh>::face_descriptor fd,
                                CGAL::Random& rnd)
{
  const Triangle_mesh& mesh = motorcycle_graph.mesh();

  boost::property_map<Triangle_mesh, CGAL::vertex_point_t>::const_type vpm = get(CGAL::vertex_point, mesh);

  Triangle_2 tr(get(vpm, source(halfedge(fd, mesh), mesh)),
                get(vpm, target(halfedge(fd, mesh), mesh)),
                get(vpm, target(next(halfedge(fd, mesh), mesh), mesh)));
  return random_motorcycles_in_triangle(motorcycle_graph, tr, rnd);
}

void clear_data_files()
{
  std::ofstream oof, odf;

  oof.open("results_2/motorcycles_origins.xyz", std::ofstream::out | std::ofstream::trunc);
  odf.open("results_2/motorcycles_destinations.xyz", std::ofstream::out | std::ofstream::trunc);
}

int main()
{
  clear_data_files();

  std::cout << std::fixed;
  std::cout.precision(17);
  std::cerr << std::fixed;
  std::cerr.precision(17);

#if 1//def CGAL_MOTORCYCLE_GRAPH_USE_FIXED_SEEDS
  CGAL::Random rnd(1522335154); // THIS IS A BUGGY CASE @FIXME (clubs 1&3)
  //CGAL::Random rnd(1518004014); // THIS IS A BUGGY CASE @FIXME
#else
  CGAL::Random rnd(CGAL::get_default_random());
#endif
  std::ofstream seed_out("results_2/seed.txt");

  // read input mesh
  Triangle_mesh tm;
  std::ifstream in("data/eight_triangles.off");
  if(!in.good())
  {
    std::cerr << "Error: Failed to read input mesh" << std::endl;
    return EXIT_FAILURE;
  }

  in >> tm;
  CGAL_precondition(tm.number_of_faces() && tm.is_valid());

  bool is_loop_infinite = true; // loop till it breaks !
  while(is_loop_infinite)
  {
    is_loop_infinite = false;

    seed_out << rnd.get_seed() << std::endl;

    Motorcycle_graph motorcycle_graph(CP::input_mesh(&tm));

    motorcycle_club_0(motorcycle_graph);
//    motorcycle_club_1(motorcycle_graph);
    motorcycle_club_2(motorcycle_graph);
//    motorcycle_club_3(motorcycle_graph);
//    motorcycle_club_4(motorcycle_graph);
//    motorcycle_club_5(motorcycle_graph);
//    motorcycle_club_6(motorcycle_graph);

    random_motorcycles_on_face(motorcycle_graph, *(faces(tm).begin()), rnd);
    random_motorcycles_on_face(motorcycle_graph, *(++(++(++(++(++faces(tm).begin()))))), rnd);
    random_motorcycles_on_face(motorcycle_graph, *(++(++faces(tm).begin())), rnd);

    motorcycle_graph.construct_motorcycle_graph();

    rnd = CGAL::Random();
  }

  return EXIT_SUCCESS;
}
