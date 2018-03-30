#define CGAL_CHECK_EXPENSIVE

//#define CGAL_MOTORCYCLE_GRAPH_VERBOSE
//#define CGAL_MOTORCYCLE_GRAPH_OUTPUT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyline_tracing/Motorcycle_graph.h>
#include <CGAL/Polyline_tracing/Motorcycle_graph_traits_3.h>
#include <CGAL/Polyline_tracing/Point_set_tracer.h>
//#include <CGAL/Polyline_tracing/Uniform_direction_tracer_visitor.h>

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

//@todo extend uniform to surfaces
//typedef PL::Uniform_direction_tracer_visitor<MGT>                Uniform_tracer;
//typedef PL::Motorcycle<MGT, Uniform_tracer>                      U_Motorcycle;
typedef PL::Point_set_tracer<MGT>                                PS_tracer;
typedef PL::Motorcycle<MGT, PS_tracer>                           PS_Motorcycle;

typedef PL::Motorcycle_graph<MGT>                                Motorcycle_graph;

typedef Motorcycle_graph::Motorcycle                             Motorcycle;
typedef boost::shared_ptr<Motorcycle>                            Motorcycle_ptr;
typedef std::vector<Motorcycle_ptr>                              Motorcycle_container;

// a bunch of motorcycle all starting from the same mesh vertex
// useful to test the ordering of edges in the ouput graph (@todo move to 'tests')
void motorcycle_club_1(Motorcycle_container& motorcycles, const PolygonMesh& mesh)
{
  boost::graph_traits<PolygonMesh>::halfedge_descriptor hd = *(halfedges(mesh).begin());
  boost::graph_traits<PolygonMesh>::vertex_descriptor vd = target(hd, mesh);

  CGAL::Face_around_target_iterator<PolygonMesh> fit, fend;
  boost::tie(fit, fend) = CGAL::faces_around_target(hd, mesh);

  while(fit != fend)
  {
    boost::graph_traits<PolygonMesh>::face_descriptor fd = *fit++;
    CGAL::cpp11::array<FT, 3> arr = {{0.,0.,0.}};

    int v_id = CGAL::Polygon_mesh_processing::vertex_index_in_face(vd, fd, mesh);
    arr[v_id] = 1.;
    Face_location loc = std::make_pair(fd, arr);

    arr[v_id] = 0.;
    arr[(v_id+1)%3] = 1.;
    Face_location loc2 = std::make_pair(fd, arr);
    PS_tracer pst;
    pst.add_destination(loc2);

    arr[(v_id+1)%3] = 0.3;
    arr[(v_id+2)%3] = 0.7;
    Face_location loc3 = std::make_pair(fd, arr);
    PS_tracer pst2;
    pst2.add_destination(loc3);

    motorcycles.push_back(Motorcycle_ptr(new PS_Motorcycle(CP::source = loc,
                                                           CP::tracer = pst)));

    motorcycles.push_back(Motorcycle_ptr(new PS_Motorcycle(CP::source = loc,
                                                           CP::tracer = pst2)));
  }
}

// -----------------------------------------------------------------------------

void random_motorcycle_club(Motorcycle_container& motorcycles,
                            const PolygonMesh& mesh,
                            CGAL::Random& rnd)
{
  // This club has a single motorcycle starting from the barycenter of each face
  typedef CGAL::property_map_selector<PolygonMesh, boost::vertex_point_t>::const_type VertexPointMap;
  VertexPointMap vpmap = get_const_property_map(CGAL::vertex_point, mesh);

  std::size_t max_number_of_motorcycles = 3000;

  boost::graph_traits<PolygonMesh>::face_iterator first, fit, last;
  boost::tie(first, last) = faces(mesh);
  fit = first;
  --last;
  for(;;)
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

    fit = (fit == last) ? first : ++fit;
  }

  std::cout << motorcycles.size() << " motorcycles enter the game" << std::endl;
}

int main()
{
  std::cout.precision(18);
  std::cerr.precision(18);

#if 1//def CGAL_MOTORCYCLE_GRAPH_USE_FIXED_SEEDS
  CGAL::Random rnd(1522325815);
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

//  motorcycle_club_1(motorcycles, pm);
  random_motorcycle_club(motorcycles, pm, rnd);

  Motorcycle_graph motorcycle_graph(pm);
  motorcycle_graph.construct_motorcycle_graph(motorcycles.begin(), motorcycles.end());

#if 0//def CGAL_MOTORCYCLE_GRAPH_OUTPUT
  motorcycle_graph.output_all_dictionary_points();
  for(std::size_t i=0; i<motorcycles.size(); ++i)
  {
//    motorcycle_graph.motorcycle(i).output_intended_track();
    motorcycle_graph.motorcycle(i).output_track();
  }
#endif

  return 0;
}
