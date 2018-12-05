// Testing intersections of multiple motorcycles with different speeds

#define CGAL_CHECK_EXPENSIVE
#define CGAL_MOTORCYCLE_GRAPH_VERBOSE
#define CGAL_MOTORCYCLE_GRAPH_OUTPUT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/polyline_tracing.h>

#include <CGAL/number_utils.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

namespace CP = CGAL::parameters;
namespace PL = CGAL::Polyline_tracing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;

typedef K::FT                                                    FT;
typedef K::Point_2                                               Point_2;
typedef K::Vector_2                                              Vector_2;

typedef CGAL::Surface_mesh<Point_2>                              Triangle_mesh;

typedef PL::Motorcycle_graph_traits_2<K, Triangle_mesh>          MGT;
typedef PL::Motorcycle<MGT>                                      Motorcycle;
typedef PL::Motorcycle_graph<MGT, Motorcycle>                    Motorcycle_graph;

typedef PL::Uniform_direction_tracer_visitor<Motorcycle_graph>   Uniform_tracer;
typedef PL::Point_set_tracer<Motorcycle_graph>                   PS_tracer;

typedef MGT::Face_location                                       Face_location;

typedef boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;

void motorcycle_club(Motorcycle_graph& motorcycle_graph)
{
  // Motorcycle #0
  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.1, 0.2))
                                     .speed(0.1));
  assert(motorcycle_graph.motorcycles().back().speed() == 0.1);

  // Motorcycle #1
  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0., 1.))
                                     .speed(1.));

  // Motorcycle #2
  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0., -1.))
                                     .speed(2.));

  // Motorcycle #3
  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(Vector_2(1., -1.)),
                                  CP::speed(0.5*CGAL::sqrt(2.)));

  // Motorcycle #4 & #5
  motorcycle_graph.add_motorcycle(Point_2(0., 1.) /*origin*/,
                                  Uniform_tracer(Vector_2(-0.5, -0.5)),
                                  CP::speed(CGAL::sqrt(2.)));
  motorcycle_graph.add_motorcycle(Point_2(-1., -1.) /*origin*/,
                                  Uniform_tracer(Vector_2(0., 1.))); // speed = 1
  assert(motorcycle_graph.motorcycles().back().speed() == 1); // test the default value

  // Motorcycle #6 & #7
  face_descriptor fd0 = Triangle_mesh::Face_index(0);
  motorcycle_graph.add_motorcycle(Face_location(fd0, CGAL::make_array(0.8, 0.2, 0.)) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Face_location(fd0, CGAL::make_array(0.2, 0.8, 0.)))
                                     .speed(2.));
  face_descriptor fd1 = Triangle_mesh::Face_index(1);
  motorcycle_graph.add_motorcycle(Face_location(fd1, CGAL::make_array(0., 0.8, 0.2)) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Face_location(fd1, CGAL::make_array(0., 0.2, 0.8)))
                                     .speed(1.));

  // Motorcycle #8 & #9
  motorcycle_graph.add_motorcycle(Point_2(0.9, 0.75) /*origin*/,
                                  Uniform_tracer(Vector_2(0.1, -1.)),
                                  CP::speed(3.));
  motorcycle_graph.add_motorcycle(Point_2(0.5, -0.1) /*origin*/,
                                  Uniform_tracer(Vector_2(1., 0.))); // speed = 1

  // Motorcycle #-2 and #-1, fail (intentionally) on an assertion
//  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
//                                  Uniform_tracer(),
//                                  CP::destination(Point_2(0.1, 0.2))
//                                     .speed(-1.));
//  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
//                                  Uniform_tracer(),
//                                  CP::destination(Point_2(0.1, 0.2))
//                                     .speed(0.));
}

void is_valid_graph(const Motorcycle_graph& motorcycle_graph)
{
  // Motorcycle #1 should reach its final destination at t=1
  const Motorcycle& mc1 = motorcycle_graph.motorcycle(1);
  assert(mc1.track().back().time_at_target() == 1.);

  // Motorcycle #2 should reach its final destination at t=0.5
  const Motorcycle& mc2 = motorcycle_graph.motorcycle(2);
  assert(mc2.track().back().time_at_target() == 0.5);

  // Motorcycle #3 should reach its final destination at t=2
  const Motorcycle& mc3 = motorcycle_graph.motorcycle(3);
  assert(mc3.track().back().time_at_target() == 2.);

  // Motorcycle #4 & #5 should crash into each other at (-1., 0.)
  const Motorcycle& mc4 = motorcycle_graph.motorcycle(4);
  const Motorcycle& mc5 = motorcycle_graph.motorcycle(5);
  assert(mc4.track().back().target()->point() == Point_2(-1., 0.));
  assert(mc4.track().back().target()->point() == mc5.track().back().target()->point());
  assert(mc4.track().back().time_at_target() == mc5.track().back().time_at_target());

  // Motorcycle #6 & #7 should crash into each other at 4/6th of the edge
  const Motorcycle& mc6 = motorcycle_graph.motorcycle(6);
  const Motorcycle& mc7 = motorcycle_graph.motorcycle(7);
  assert(mc6.track().back().target()->point() == mc7.track().back().target()->point());
  assert(CGAL::abs(mc6.track().back().target()->barycentric_coordinate(0) - 0.4) < std::numeric_limits<FT>::epsilon());
  assert(CGAL::abs(mc6.track().back().target()->barycentric_coordinate(1) - 0.6) < std::numeric_limits<FT>::epsilon());

  // Motorcycle #8 cuts the track of motorcycle #9
  const Motorcycle& mc8 = motorcycle_graph.motorcycle(8);
  const Motorcycle& mc9 = motorcycle_graph.motorcycle(9);
  assert(mc8.current_position() == mc8.destination());
  assert(mc9.current_position() != mc9.destination());
}

void clear_data_files()
{
  std::ofstream oof, odf;

  // @fixme results_2 not acceptable
  oof.open("results_2/motorcycles_origins.xyz", std::ofstream::out | std::ofstream::trunc);
  odf.open("results_2/motorcycles_destinations.xyz", std::ofstream::out | std::ofstream::trunc);
}

int main()
{
  std::cout.precision(17);
  std::cout << std::fixed;
  std::cerr.precision(17);
  std::cout << std::fixed;

  clear_data_files();

  // read input mesh
  Triangle_mesh tm;
  std::ifstream in("data/eight_triangles.off");
  in >> tm;
  CGAL_precondition(tm.is_valid());

  Motorcycle_graph motorcycle_graph(CP::input_mesh(&tm));
  motorcycle_club(motorcycle_graph);

  motorcycle_graph.construct_motorcycle_graph();

  assert(motorcycle_graph.is_valid());
  assert(CGAL::Polyline_tracing::internal::is_valid_hds(motorcycle_graph.graph()));
  is_valid_graph(motorcycle_graph);

  return EXIT_SUCCESS;
}
