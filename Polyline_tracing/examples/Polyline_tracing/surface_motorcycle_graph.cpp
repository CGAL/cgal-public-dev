#define CGAL_MOTORCYCLE_GRAPH_VERBOSE
#define CGAL_MOTORCYCLE_GRAPH_OUTPUT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyline_tracing/Motorcycle_graph.h>
#include <CGAL/Polyline_tracing/Motorcycle_graph_traits_3.h>

#include <fstream>
#include <iostream>
#include <vector>

namespace PL = CGAL::Polyline_tracing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef CGAL::Surface_mesh<K::Point_3>                           PolygonMesh;
typedef PL::Motorcycle_graph_traits_3<K, PolygonMesh>            MGT;

typedef MGT::FT                                                  FT;
typedef MGT::Point_d                                             Point_3;
typedef MGT::Vector_d                                            Vector_3;

typedef PL::Motorcycle_graph<MGT>                                Motorcycle_graph;
typedef Motorcycle_graph::Motorcycle                             Motorcycle;


namespace CP = CGAL::parameters;

void motorcycle_club_1(std::vector<Motorcycle>& motorcycles)
{
  // This is a casual motorcycle club that just likes to use all of its options

  motorcycles.push_back(Motorcycle(CP::speed = 1.,
                                   CP::source = Point_3(0.1, 0.1, 0),
                                   CP::direction = Vector_3(1, 0, 0),
                                   CP::initial_time = 0.));
  motorcycles.push_back(Motorcycle(CP::source = Point_3(0.1, 0.8, 0),
                                   CP::direction = Vector_3(1., -1, 0),
                                   CP::speed = 1.));
}

int main()
{
  std::cout.precision(18);
  std::cerr.precision(18);

  PolygonMesh pm;
  std::ifstream in("data/triangle_3.off");
  in >> pm;
  CGAL_precondition(pm.is_valid());

//  std::ofstream out("polygon_mesh.off");
//  out << pm;

  std::vector<Motorcycle> motorcycles;
  motorcycle_club_1(motorcycles);

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

  return 0;
}
