#define CGAL_CHECK_EXPENSIVE

#define CGAL_MOTORCYCLE_GRAPH_VERBOSE
#define CGAL_MOTORCYCLE_GRAPH_OUTPUT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyline_tracing/Motorcycle_graph.h>
#include <CGAL/Polyline_tracing/Motorcycle_graph_traits_3.h>
#include <CGAL/Polyline_tracing/Point_set_tracer.h>
// #include <CGAL/Polyline_tracing/Uniform_direction_tracer_visitor.h>

#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Random.h>

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

#include <fstream>
#include <iostream>
#include <vector>

namespace CP = CGAL::parameters;
namespace PL = CGAL::Polyline_tracing;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef CGAL::Surface_mesh<K::Point_3>                           PolygonMesh;
typedef PL::Motorcycle_graph_traits_3<K, PolygonMesh>            MGT;

typedef MGT::FT                                                  FT;
typedef MGT::Point_d                                             Point_3;
typedef MGT::Vector_d                                            Vector_3;
typedef MGT::Face_location                                       Face_location;

typedef boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;

//typedef PL::Uniform_direction_tracer_visitor<MGT>                Uniform_tracer;
//typedef PL::Motorcycle<MGT, Uniform_tracer>                      Motorcycle_U;
typedef PL::Point_set_tracer<MGT>                                PS_tracer;
typedef PL::Motorcycle<MGT, PS_tracer>                           PS_Motorcycle;

typedef PL::Motorcycle_graph<MGT>                                Motorcycle_graph;

typedef Motorcycle_graph::Motorcycle                             Motorcycle;
typedef boost::shared_ptr<Motorcycle>                            Motorcycle_ptr;
typedef std::vector<Motorcycle_ptr>                              Motorcycle_container;

void motorcycle_club_1(Motorcycle_container& /*motorcycles*/)
{
  // This is a casual motorcycle club that just likes to use all of its options

//  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::speed = 1.,
//                                                        CP::source = Point_3(0.1, 0.1, 0),
//                                                        CP::direction = Vector_3(1, 0, 0),
//                                                        CP::initial_time = 0.)));
//  motorcycles.push_back(Motorcycle_ptr(new Motorcycle_U(CP::source = Point_3(0.1, 0.8, 0),
//                                                        CP::direction = Vector_3(1., -1, 0),
//                                                        CP::speed = 1.)));
}

// -----------------------------------------------------------------------------

void random_motorcycle_club(Motorcycle_container& motorcycles,
                            const PolygonMesh& mesh,
                            CGAL::Random& rnd)
{
  // This club has a single motorcycle starting from the barycenter of each face
  typedef CGAL::property_map_selector<PolygonMesh, boost::vertex_point_t>::const_type VertexPointMap;
  VertexPointMap vpmap = get_const_property_map(CGAL::vertex_point, mesh);

  std::size_t max_number_of_motorcycles = 50;

  boost::graph_traits<PolygonMesh>::face_iterator fit, end;
  boost::tie(fit, end) = faces(mesh);
  for(; fit!=end; ++fit)
  {
    // Motorcycle
    const FT third = 1./3.;
    halfedge_descriptor hd = halfedge(*fit, mesh);
    const Point_3 bar = CGAL::barycenter(get(vpmap, source(hd, mesh)), third,
                                         get(vpmap, target(hd, mesh)), third,
                                         get(vpmap, target(next(hd, mesh), mesh)), third);

    // Generate some targets
    PS_tracer pst;
    std::size_t dest_n = 20;
    std::vector<Face_location> destinations;
    destinations.reserve(dest_n);
    halfedge_descriptor entry_hd = hd, exit_hd = boost::graph_traits<PolygonMesh>::null_halfedge();

    for(std::size_t i=0; i<dest_n; ++i)
    {
      exit_hd = next(entry_hd, mesh);

      if(rnd.uniform_int(0,2) == 0)
        exit_hd = next(exit_hd, mesh);

      Face_location loc = PMP::random_location_on_halfedge(exit_hd, mesh, rnd);
      destinations.push_back(loc);
      entry_hd = opposite(exit_hd, mesh);
    }
    pst.set_destinations(destinations);

    motorcycles.push_back(Motorcycle_ptr(new PS_Motorcycle(CP::source = bar,
                                                           CP::tracer = pst)));

    if(motorcycles.size() >= max_number_of_motorcycles)
      break;
  }
}

int main()
{
  std::cout.precision(18);
  std::cerr.precision(18);

#ifdef CGAL_MOTORCYCLE_GRAPH_USE_FIXED_SEEDS
  CGAL::Random rnd(1506422869);
#else
  CGAL::Random rnd(CGAL::get_default_random());
#endif
  std::ofstream seed_out("results_3/seed.txt");
  seed_out << rnd.get_seed() << std::endl;

  // read input mesh
  PolygonMesh pm;
  std::ifstream in("data/rough_bunny.off");
  in >> pm;
  CGAL_precondition(pm.is_valid());

  Motorcycle_container motorcycles;

//  motorcycle_club_1(motorcycles);
  random_motorcycle_club(motorcycles, pm, rnd);

  Motorcycle_graph motorcycle_graph(pm);
  motorcycle_graph.construct_motorcycle_graph(motorcycles.begin(), motorcycles.end());

#ifdef CGAL_MOTORCYCLE_GRAPH_OUTPUT
  motorcycle_graph.output_all_dictionary_points();
  for(std::size_t i=0; i<motorcycles.size(); ++i)
  {
//    motorcycle_graph.motorcycle(i).output_intended_track();
    motorcycle_graph.motorcycle(i).output_track();
  }
#endif

  return 0;
}
