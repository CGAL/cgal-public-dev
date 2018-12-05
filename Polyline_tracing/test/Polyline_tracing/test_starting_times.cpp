// Testing motorcycles with different starting times

#define CGAL_CHECK_EXPENSIVE
#define CGAL_MOTORCYCLE_GRAPH_VERBOSE
#define CGAL_MOTORCYCLE_GRAPH_OUTPUT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/polyline_tracing.h>

#include <fstream>
#include <iostream>
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
  motorcycle_graph.add_motorcycle(Point_2(0., 1.) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(-1., 0.))
                                     .speed(0.5 * CGAL::sqrt(2.))
                                     .initial_time(-2.));

  // Motorcycle #1
  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(Vector_2(-0.1, 0.)),
                                  CP::initial_time(-1.));
  assert(motorcycle_graph.motorcycles().back().current_time() == -1.);

  // Motorcycle #2
  motorcycle_graph.add_motorcycle(Point_2(-1., 1.) /*origin*/,
                                  Uniform_tracer(Vector_2(0., -0.1)),
                                  CP::initial_time(-1.));

  // Motorcycle #3
  motorcycle_graph.add_motorcycle(Point_2(-0.8, 0.1) /*origin*/,
                                  Uniform_tracer(Vector_2(-0.1, 0.1)));
  assert(motorcycle_graph.motorcycles().back().current_time() == 0.);

  // Motorcycle #4, #5, #6 and #7
  face_descriptor fd0 = Triangle_mesh::Face_index(0);
  face_descriptor fd1 = Triangle_mesh::Face_index(1);

  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(Vector_2(-1., -1.)),
                                  CP::initial_time(1.)
                                     .speed(CGAL::sqrt(2.)));
  motorcycle_graph.add_motorcycle(Point_2(-0.2, -0.2) /*origin*/,
                                  Uniform_tracer(Vector_2(-1., 1.)),
                                  CP::initial_time(2.));
  motorcycle_graph.add_motorcycle(Face_location(fd1, CGAL::make_array(0., 0.5, 0.5)) /*origin*/,
                                  Uniform_tracer(Vector_2(0., -1.)),
                                  CP::initial_time(1.5));
  motorcycle_graph.add_motorcycle(Face_location(fd0, CGAL::make_array(0.2, 0.8, 0.)) /*origin*/,
                                  Uniform_tracer(Vector_2(1., -1.)),
                                  CP::initial_time(1.));

  // Motorcycles #8, #9, and #10
  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(Vector_2(0., -1.)),
                                  CP::initial_time(-5.));
  motorcycle_graph.add_motorcycle(Point_2(0., -0.25) /*origin*/,
                                  Uniform_tracer(Vector_2(0., -1.)),
                                  CP::initial_time(0.));
  motorcycle_graph.add_motorcycle(Point_2(0., -0.75) /*origin*/,
                                  Uniform_tracer(Vector_2(0., 1.)));
  assert(motorcycle_graph.motorcycles().back().current_time() == 0.); // tests the default value

  // Motorcycles #11, #12, #13
  motorcycle_graph.add_motorcycle(Point_2(0.25, -0.125) /*origin*/,
                                  Uniform_tracer(Vector_2(1., 0.)),
                                  CP::initial_time(-1.));
  motorcycle_graph.add_motorcycle(Point_2(0.3, -0.125) /*origin*/,
                                  Uniform_tracer(Vector_2(1., 0.)),
                                  CP::initial_time(1.));
  motorcycle_graph.add_motorcycle(Point_2(0.5, -0.125) /*origin*/,
                                  Uniform_tracer(Vector_2(-1., 0.)),
                                  CP::initial_time(1.));

  // Motorcycle #14
  motorcycle_graph.add_motorcycle(Point_2(-0.1, 0.1) /*origin*/,
                                  Uniform_tracer(Vector_2(1., -1.)),
                                  CP::initial_time(-5.1)
                                     .speed(CGAL::sqrt(2.)));

  // Motorcycle #15
  face_descriptor fd2 = Triangle_mesh::Face_index(2);
  motorcycle_graph.add_motorcycle(Face_location(fd2, CGAL::make_array(0.5, 0.5, 0.)) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.9,-0.9))
                                     .initial_time(-5.)
                                     .speed(CGAL::sqrt(2.)));

  // Motorcycle #16
  motorcycle_graph.add_motorcycle(Point_2(0.95, -0.95) /*origin*/,
                                  Uniform_tracer(Vector_2(-1., 1.)),
                                  CP::initial_time(-4.55)
                                      .speed(CGAL::sqrt(2.)));

  // Motorcycle #17, #18, #19, #20, #21
  motorcycle_graph.add_motorcycle(Point_2(0.15, 0.65) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.65, 0.65))
                                      .initial_time(2.));
  motorcycle_graph.motorcycles().back().is_destination_final() = true;

  motorcycle_graph.add_motorcycle(Point_2(0.15, 0.65) /*origin*/,
                                  Uniform_tracer(Vector_2(1., 0.)),
                                  CP::initial_time(10.));
  motorcycle_graph.add_motorcycle(Point_2(0.15, 0.65) /*origin*/,
                                  Uniform_tracer(Vector_2(-1., 0.)),
                                  CP::initial_time(10.));
  motorcycle_graph.add_motorcycle(Point_2(0.65, 0.65) /*origin*/,
                                  Uniform_tracer(Vector_2(-1., 0.)),
                                  CP::initial_time(10.));
  motorcycle_graph.add_motorcycle(Point_2(0.65, 0.65) /*origin*/,
                                  Uniform_tracer(Vector_2(1., 0.)),
                                  CP::initial_time(10.));
}

void is_valid_graph(const Motorcycle_graph& motorcycle_graph)
{
  // Motorcycle #0 reaches its destination at t=0
  const Motorcycle& mc0 = motorcycle_graph.motorcycle(0);
  assert(mc0.track().back().time_at_target() == 0.);
  assert(mc0.track().back().target()->point() == Point_2(-1., 0.));

  // Motorcycle #1 reaches its destination (same as #0) at t=0 and both crashes
  const Motorcycle& mc1 = motorcycle_graph.motorcycle(1);
  assert(mc1.track().back().time_at_target() == 0.);
  assert(mc0.track().back().target()->point() == mc1.track().back().target()->point());

  // Motorcycle #2 crashes in the same point as #0 and #1 but it was not its final destination
  const Motorcycle& mc2 = motorcycle_graph.motorcycle(2);
  assert(mc2.track().back().time_at_target() == 0.);
  assert(mc2.track().back().target()->point() == mc1.track().back().target()->point());

  // Motorcycle #3 crashes into Motorcycle #0, because it starts later
  const Motorcycle& mc3 = motorcycle_graph.motorcycle(3);
  assert(CGAL::squared_distance(mc3.track().back().target()->point(), Point_2(-0.85, 0.15)) <
           std::numeric_limits<FT>::epsilon());

  // Motorcycles #5, #6, and #7 start respectively after, at the same time, and before
  // motorcycle #4 passes. Only the motorcycle starting before (#7) blocks #4.
  const Motorcycle& mc4 = motorcycle_graph.motorcycle(4);
  const Motorcycle& mc5 = motorcycle_graph.motorcycle(5);
  const Motorcycle& mc6 = motorcycle_graph.motorcycle(6);
  const Motorcycle& mc7 = motorcycle_graph.motorcycle(7);
  assert(mc4.current_position() == mc7.track().front().source());
  assert(mc5.current_position() != mc5.track().front().source());
  assert(mc6.current_position() != mc6.track().front().source());
  assert(mc6.track().back().target()->point().y() == -1.);

  // Motorcycles #9 and #10 start on an already existing track(#7's) on an edge,
  // moving in a collinear direction and should crash instantly.
  const Motorcycle& mc8 = motorcycle_graph.motorcycle(8);
  assert(mc8.current_position()->point() == Point_2(0., -1.));
  const Motorcycle& mc9 = motorcycle_graph.motorcycle(9);
  assert(mc9.track().size() == 1 && mc9.current_position() == mc9.track().front().source());
  const Motorcycle& mc10 = motorcycle_graph.motorcycle(10);
  assert(mc10.track().size() == 1 && mc10.current_position() == mc10.track().front().source());

  // Motorcycles #12 and #13 start on an already existing track (#11's), moving
  // in a collinear direction and should crash instantly.
  const Motorcycle& mc11 = motorcycle_graph.motorcycle(11);
  assert(mc11.current_position()->point() == Point_2(1., -0.125));
  const Motorcycle& mc12 = motorcycle_graph.motorcycle(12);
  assert(mc12.current_position()->point() == Point_2(0.3, -0.125));
  const Motorcycle& mc13 = motorcycle_graph.motorcycle(13);
  assert(mc13.current_position() == mc13.track().front().source());

  // Motorcycle #14 goes through the origin and should not block anyone,
  // but it will crash on motorcycle #15 that is on its track and in the same direction
  const Motorcycle& mc14 = motorcycle_graph.motorcycle(14);
  const Motorcycle& mc15 = motorcycle_graph.motorcycle(15);
  assert(mc14.current_position()->point() == mc15.track().front().source()->point());

  // Motorcycle #16 starts in the opposite direction as motorcycle #15 just as #15
  // reaches it. Motorcycle #16 should instantly crash but #15 continues.
  const Motorcycle& mc16 = motorcycle_graph.motorcycle(16);
  assert(mc15.current_position()->point() == Point_2(1., -1.));
  assert(mc16.current_position() == mc16.track().front().source());
  assert(mc16.track().size() == 1);

  // Motorcycle #17's track is a segment, with 2 motorcycles starting from
  // its extremities in collinear, both directions, much later.
  const Motorcycle& mc17 = motorcycle_graph.motorcycle(17);
  const Motorcycle& mc18 = motorcycle_graph.motorcycle(18);
  const Motorcycle& mc19 = motorcycle_graph.motorcycle(19);
  const Motorcycle& mc20 = motorcycle_graph.motorcycle(20);
  const Motorcycle& mc21 = motorcycle_graph.motorcycle(21);
  assert(mc17.track().size() == 2);
  assert(mc18.current_position() == mc18.track().front().source()); // crashes
  assert(mc19.current_position() != mc19.track().front().source()); // doesn't crash
  assert(mc20.current_position() == mc20.track().front().source()); // crashes
  assert(mc21.current_position() != mc21.track().front().source()); // doesn't crash
}

void clear_data_files()
{
  std::ofstream oof, odf;

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
  assert(CGAL::Polyline_tracing::internal::is_valid(motorcycle_graph));

  is_valid_graph(motorcycle_graph);

  return EXIT_SUCCESS;
}
