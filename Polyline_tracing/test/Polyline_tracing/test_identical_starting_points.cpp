// All motorcycles here are pairs of motorcycles starting from various locations.
// For each pair, the motorcycles have the same starting point, time and directions.
// They all must crash instantly.

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

// motorcycles starting from a vertex
void vertex_motorcycle_club(Motorcycle_graph& motorcycle_graph)
{
  // #0/#1 exactly the same starting point and destinations
  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.1, 0.2)));
  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.1, 0.2)));

  // #2 add another one on top!
  // @fixme does not work because #0 and #1 crash each other and then nothing prevents
  // this one from moving forward. Don't know (yet) how to handle this _very_ specific case
  // without changing the base algorithm.
//  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
//                                  Uniform_tracer(),
//                                  CP::destination(Point_2(0.1, 0.2)));

  // #3/#4 same starting point as above but with another pair of directions
  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(-0.5, 0.5)));
  motorcycle_graph.add_motorcycle(Point_2(0., 0.) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(-0.5, 0.5)));

  // #5/#6 same starting point, but with different direction types
  motorcycle_graph.add_motorcycle(Point_2(1., 1.) /*origin*/,
                                  Uniform_tracer(Vector_2(-0.5, -1.)));
  motorcycle_graph.add_motorcycle(Point_2(1., 1.) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.9, 0.8)));

  // #7/#8 same starting point (-1., 0.), but on different faces (edge walking)
  face_descriptor fd1 = Triangle_mesh::Face_index(1);
  Face_location first_loc = std::make_pair(fd1, CGAL::make_array(1., 0., 0.));

  face_descriptor fd4 = Triangle_mesh::Face_index(4);
  Face_location second_loc = std::make_pair(fd4, CGAL::make_array(0., 1., 0.));

  motorcycle_graph.add_motorcycle(first_loc /*origin*/,
                                  Uniform_tracer(Vector_2(0.5, 0.5)));
  motorcycle_graph.add_motorcycle(second_loc /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(-0.5, 0.5)));

  // #9/#10 same starting point, different tracers
  motorcycle_graph.add_motorcycle(Point_2(0., 1.) /*origin*/,
                                  Uniform_tracer(Vector_2(0., -1.)));

  PS_tracer pst;
  pst.add_destination(Face_location(fd1, CGAL::make_array(0., 0., 1.)));
  motorcycle_graph.add_motorcycle(Point_2(0., 1.) /*origin*/, pst);
}

// motorcycles starting from a halfedge
void edge_motorcycle_club(Motorcycle_graph& motorcycle_graph)
{
  // #11/#12 exactly the same starting point and destinations
  motorcycle_graph.add_motorcycle(Point_2(0.5, 0.5) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0., 1.)));
  motorcycle_graph.add_motorcycle(Point_2(0.5, 0.5) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0., 1.)));

  // #13 add another one on top!
//  motorcycle_graph.add_motorcycle(Point_2(0.5, 0.5) /*origin*/,
//                                  Uniform_tracer(),
//                                  CP::destination(Point_2(0., 1.)));

  // #14/#15 same starting point as above but with another pair of directions
  motorcycle_graph.add_motorcycle(Point_2(0.5, 0.5) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0., 1.)));
  motorcycle_graph.add_motorcycle(Point_2(0.5, 0.5) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0., 1.)));

  // #16/#17 same starting point, but with different direction types
  motorcycle_graph.add_motorcycle(Point_2(-0.5, -0.5) /*origin*/,
                                  Uniform_tracer(Vector_2(-0.5, -0.5)));
  motorcycle_graph.add_motorcycle(Point_2(-0.5, -0.5) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(-1., -1.)));

  // #18/#19 same starting point (0.3, -0.7), but on different faces (edge walking)
  face_descriptor fd2 = Triangle_mesh::Face_index(2);
  Face_location first_loc_o = std::make_pair(fd2, CGAL::make_array(0.3, 0.7, 0.));

  face_descriptor fd3 = Triangle_mesh::Face_index(3);
  Face_location second_loc_o = std::make_pair(fd3, CGAL::make_array(0., 0.7, 0.3));
  Face_location second_loc_d = std::make_pair(fd3, CGAL::make_array(0., 0.5, 0.5));

  motorcycle_graph.add_motorcycle(first_loc_o /*origin*/,
                                  Uniform_tracer(Vector_2(0.5, 0.5)));
  motorcycle_graph.add_motorcycle(second_loc_o /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(second_loc_d));

  // #20/#21 same starting point, different tracers
  motorcycle_graph.add_motorcycle(Point_2(-0.8, 0.2) /*origin*/,
                                  Uniform_tracer(Vector_2(1., 4.)));

  PS_tracer pst;
  face_descriptor fd5 = Triangle_mesh::Face_index(5);
  pst.add_destination(Face_location(fd5, CGAL::make_array(0.6, 0., 0.4)));
  motorcycle_graph.add_motorcycle(Point_2(-0.8, 0.2) /*origin*/, pst);
}

// motorcycles starting within a face
void face_motorcycle_club(Motorcycle_graph& motorcycle_graph)
{
  // #22/#23 exactly the same starting point and destinations
  motorcycle_graph.add_motorcycle(Point_2(0.1, -0.1) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.8, -1./11.)));
  motorcycle_graph.add_motorcycle(Point_2(0.1, -1./10.) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.8, -1./11.)));

  // #24 add another one on top!
//  motorcycle_graph.add_motorcycle(Point_2(0.1, -0.1) /*origin*/,
//                                  Uniform_tracer(),
//                                  CP::destination(Point_2(0.8, -1./11.)));

  // #25/#26 same starting point as above but with another pair of directions
  motorcycle_graph.add_motorcycle(Point_2(0.1, -0.1) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.5, -0.5)));
  motorcycle_graph.add_motorcycle(Point_2(0.1, -0.1) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(0.5, -0.5)));

  // #27/#28 same starting point, but with different direction types
  motorcycle_graph.add_motorcycle(Point_2(-0.3, -0.2) /*origin*/,
                                  Uniform_tracer(Vector_2(0., 1.)));
  motorcycle_graph.add_motorcycle(Point_2(-0.3, -0.2) /*origin*/,
                                  Uniform_tracer(),
                                  CP::destination(Point_2(-0.3, 0.)));

  // #29/#30 same starting point, different tracers
  motorcycle_graph.add_motorcycle(Point_2(0.5, 0.8) /*origin*/,
                                  Uniform_tracer(Vector_2(1., 0.4)));

  PS_tracer pst;
  face_descriptor fd7 = Triangle_mesh::Face_index(7);
  pst.add_destination(Face_location(fd7, CGAL::make_array(0., 0., 1.)));
  motorcycle_graph.add_motorcycle(Point_2(0.5, 0.8) /*origin*/, pst);
}

void is_valid_graph(const Motorcycle_graph& motorcycle_graph)
{
  Motorcycle_graph::MCC_cit mcc_it = motorcycle_graph.motorcycles().begin(),
                            mcc_end = motorcycle_graph.motorcycles().end();
  for(; mcc_it!=mcc_end; ++mcc_it)
  {
    const Motorcycle& mc = motorcycle_graph.motorcycle(mcc_it);
    const Motorcycle_graph::Track& track = mc.track();
    assert(track.size() == 1);
    assert(track.front().is_degenerate());
    assert(track.front().source() == mc.origin());
  }
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
  assert(tm.is_valid());

  Motorcycle_graph motorcycle_graph(CP::input_mesh(&tm));
  vertex_motorcycle_club(motorcycle_graph);
  edge_motorcycle_club(motorcycle_graph);
//  face_motorcycle_club(motorcycle_graph);

  motorcycle_graph.construct_motorcycle_graph();
  assert(CGAL::Polyline_tracing::internal::is_valid(motorcycle_graph));
  assert(num_edges(motorcycle_graph.graph()) == 0);

  is_valid_graph(motorcycle_graph);

  return EXIT_SUCCESS;
}
