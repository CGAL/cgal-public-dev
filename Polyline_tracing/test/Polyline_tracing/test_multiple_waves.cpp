// Testing intersections of multiple motorcycles with different speeds

#define CGAL_CHECK_EXPENSIVE
#define CGAL_MOTORCYCLE_GRAPH_VERBOSE
#define CGAL_MOTORCYCLE_GRAPH_OUTPUT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/polyline_tracing.h>

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

void motorcycle_club_1(Motorcycle_graph& motorcycle_graph)
{
  // Motorcycle #0
  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.1, 0.2))
                                     .speed(0.3));

  // Motorcycle #1
  motorcycle_graph.add_motorcycle(Point_2(0.2, 0.3) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0., 1.)));

  // Motorcycle #2
  motorcycle_graph.add_motorcycle(Point_2(-0.7, 0.3) /*origin*/,
                                  Uniform_tracer(Vector_2(1., -0.1)));

  // Motorcycle #3
  motorcycle_graph.add_motorcycle(Point_2(-0.4, -1.) /*origin*/,
                                  Uniform_tracer(Vector_2(1., -0.1)));

  // Motorcycle #4
  motorcycle_graph.add_motorcycle(Point_2(0.7, -0.3) /*origin*/,
                                  Uniform_tracer(Vector_2(0.3, -0.2)));

  // Motorcycle #5
  motorcycle_graph.add_motorcycle(Point_2(-0.1, -0.3) /*origin*/,
                                  Uniform_tracer(Vector_2(0.1, -0.1)));
}

void motorcycle_club_2(Motorcycle_graph& motorcycle_graph)
{
  const FT time_after_first_wave = 1000; // @todo something nicer ("mg.latest_event_time()")
  std::cout << "time after first wave: " << time_after_first_wave << std::endl;

  // Motorcycle #6 & #7
  face_descriptor fd0 = Triangle_mesh::Face_index(0);
  motorcycle_graph.add_motorcycle(Face_location(fd0, CGAL::make_array(0.8, 0.2, 0.)) /*origin*/,
                                  Uniform_tracer(Vector_2(0.4, 0.6)),
                                  CP::speed(2.)
                                     .initial_time(time_after_first_wave));

  face_descriptor fd1 = Triangle_mesh::Face_index(1);
  motorcycle_graph.add_motorcycle(Face_location(fd1, CGAL::make_array(0., 0.1, 0.9)) /*origin*/,
                                  Uniform_tracer(Vector_2(0., 1.)),
                                  CP::initial_time(time_after_first_wave + 1.));

  // Motorcycle #8 & #9
  motorcycle_graph.add_motorcycle(Point_2(-0.3, -0.75) /*origin*/,
                                  Uniform_tracer(Vector_2(0.1, 1.)),
                                  CP::initial_time(time_after_first_wave)
                                     .speed(3.));

  motorcycle_graph.add_motorcycle(Point_2(0.5, -0.1) /*origin*/,
                                  Uniform_tracer(Vector_2(1., 0.)),
                                  CP::initial_time(time_after_first_wave));

  // Motorcycle #10
  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.1, 0.2))
                                     .initial_time(time_after_first_wave));

  // Motorcycle #11
  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(Vector_2(0., -1.)),
                                  CP::initial_time(time_after_first_wave));
}

bool is_valid_graph(Motorcycle_graph& motorcycle_graph)
{
  FT latest_time_of_first_wave = - std::numeric_limits<FT>::infinity();
  for(std::size_t id=0, nm=motorcycle_graph.number_of_motorcycles(); id<nm; ++id)
  {
    const Motorcycle& mc = motorcycle_graph.motorcycle(id);
    assert(mc.track().size() > 0);

    if(id < 6) // first wave
    {
      const FT arrival_time = mc.track().back().time_at_target();
      if(arrival_time > latest_time_of_first_wave)
        latest_time_of_first_wave = arrival_time;
    }
    else
    {
      const FT begin_time = mc.track().back().time_at_source();
      assert(begin_time >= latest_time_of_first_wave);
    }
  }

  // Motorcycle #10 should crash instantly
  const Motorcycle& mc10 = motorcycle_graph.motorcycle(10);
  assert(mc10.track().size() == 1.);

  // Motorcycle #11 should be able to leave its starting position
  const Motorcycle& mc11 = motorcycle_graph.motorcycle(11);
  assert(mc11.current_position()->point() != Point_2(0., 0.));

  return true;
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

  // add the first wave
  motorcycle_club_1(motorcycle_graph);
  motorcycle_graph.construct_motorcycle_graph();
  assert(CGAL::Polyline_tracing::internal::is_valid(motorcycle_graph));

  std::ofstream out("intermediary_mg.polylines.txt");
  motorcycle_graph.print_motorcycle_graph(out);

  // now with a second wave
  motorcycle_club_2(motorcycle_graph);
  motorcycle_graph.construct_motorcycle_graph();
  assert(CGAL::Polyline_tracing::internal::is_valid(motorcycle_graph));

  motorcycle_graph.output_all_points();

  return EXIT_SUCCESS;
}
