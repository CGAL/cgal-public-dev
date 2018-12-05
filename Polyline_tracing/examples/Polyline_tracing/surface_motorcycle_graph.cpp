#define CGAL_CHECK_EXPENSIVE

//#define CGAL_MOTORCYCLE_GRAPH_VERBOSE
//#define CGAL_MOTORCYCLE_GRAPH_OUTPUT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyline_tracing/Motorcycle_graph.h>
#include <CGAL/Polyline_tracing/Motorcycle_graph_traits_3.h>
#include <CGAL/Polyline_tracing/Point_set_tracer.h>

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
typedef CGAL::Surface_mesh<K::Point_3>                           Triangle_mesh;
typedef PL::Motorcycle_graph_traits_3<K, Triangle_mesh>          MGT;

typedef MGT::FT                                                  FT;
typedef MGT::Point_d                                             Point_3;
typedef MGT::Vector_d                                            Vector_3;
typedef MGT::Face_location                                       Face_location;

typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;
typedef boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;

typedef PL::Motorcycle<MGT>                                      Motorcycle;
typedef PL::Motorcycle_graph<MGT, Motorcycle>                    Motorcycle_graph;
typedef PL::Point_set_tracer<Motorcycle_graph>                   PS_tracer;

// a bunch of motorcycle all starting from the same mesh vertex
// useful to test the ordering of edges in the ouput graph
void motorcycle_club_1(Motorcycle_graph& motorcycle_graph)
{
  const Triangle_mesh& mesh = motorcycle_graph.mesh();

  halfedge_descriptor hd = *(halfedges(mesh).begin());
  vertex_descriptor vd = target(hd, mesh);

  CGAL::Face_around_target_iterator<Triangle_mesh> fit, fend;
  boost::tie(fit, fend) = CGAL::faces_around_target(hd, mesh);

  while(fit != fend)
  {
    face_descriptor fd = *fit++;
    CGAL::cpp11::array<FT, 3> arr = {{0.,0.,0.}};

    // location for the vertex
    int v_id = CGAL::Polygon_mesh_processing::vertex_index_in_face(vd, fd, mesh);
    arr[v_id] = 1.;
    Face_location origin = std::make_pair(fd, arr);

    // location on the opposite halfedge
    arr[v_id] = 0.;
    arr[(v_id+1)%3] = 1.;
    Face_location first_destination = std::make_pair(fd, arr);

    // another location on the opposite halfedge
    arr[(v_id+1)%3] = 0.3;
    arr[(v_id+2)%3] = 0.7;
    Face_location second_destination = std::make_pair(fd, arr);

    PS_tracer pst;
    pst.add_destination(first_destination);
    motorcycle_graph.add_motorcycle(origin, pst);

    PS_tracer pst2;
    pst2.add_destination(second_destination);
    motorcycle_graph.add_motorcycle(origin, pst2);
  }
}

void random_motorcycle_club(Motorcycle_graph& motorcycle_graph, CGAL::Random& rnd)
{
  const Triangle_mesh& mesh = motorcycle_graph.mesh();

  // This club has a single motorcycle starting from the barycenter of each face
  typedef CGAL::property_map_selector<Triangle_mesh, boost::vertex_point_t>::const_type VertexPointMap;
  VertexPointMap vpmap = get_const_property_map(CGAL::vertex_point, mesh);

  std::size_t max_number_of_motorcycles = 3000;

  boost::graph_traits<Triangle_mesh>::face_iterator first, fit, last;
  boost::tie(first, last) = faces(mesh);
  fit = first;
  --last;
  for(;;)
  {
    // origin
    const FT third = 1./3.;
    halfedge_descriptor hd = halfedge(*fit, mesh);
    const Point_3 origin = CGAL::barycenter(get(vpmap, source(hd, mesh)), third,
                                            get(vpmap, target(hd, mesh)), third,
                                            get(vpmap, target(next(hd, mesh), mesh)), third);

    // Generate some targets
    PS_tracer pst;
    std::size_t dest_n = 20;
    std::vector<Face_location> destinations;
    destinations.reserve(dest_n);
    halfedge_descriptor entry_hd = hd, exit_hd = boost::graph_traits<Triangle_mesh>::null_halfedge();

    for(std::size_t i=0; i<dest_n; ++i)
    {
      exit_hd = next(entry_hd, mesh);

      if(rnd.uniform_int(0,2) == 0)
        exit_hd = next(exit_hd, mesh);

      Face_location loc = PMP::random_location_on_halfedge(exit_hd, mesh, rnd);
      destinations.push_back(loc);
      entry_hd = opposite(exit_hd, mesh);
    }
    pst.destinations() = destinations;

    motorcycle_graph.add_motorcycle(origin, pst);

    if(motorcycle_graph.number_of_motorcycles() >= max_number_of_motorcycles)
      break;

    // loop the faces until the number of motorcycles is satisfactory
    fit = (fit == last) ? first : ++fit;
  }

  std::cout << motorcycle_graph.number_of_motorcycles() << " motorcycles enter the game" << std::endl;
}

int main()
{
  std::cout << std::fixed;
  std::cout.precision(17);
  std::cerr << std::fixed;
  std::cerr.precision(17);

#if 1//def CGAL_MOTORCYCLE_GRAPH_USE_FIXED_SEEDS
  CGAL::Random rnd(1522325815);
#else
  CGAL::Random rnd(CGAL::get_default_random());
#endif
  std::ofstream seed_out("results_3/seed.txt");
  seed_out << rnd.get_seed() << std::endl;

  // read input mesh
  Triangle_mesh pm;
  std::ifstream in("data/rough_bunny.off");
  if(!in.good())
  {
    std::cerr << "Error: Failed to read input" << std::endl;
    return EXIT_FAILURE;
  }

  in >> pm;
  CGAL_precondition(pm.is_valid());

  Motorcycle_graph motorcycle_graph(CP::input_mesh(&pm));

  // Add some motorcycles
//  motorcycle_club_1(motorcycle_graph);
  random_motorcycle_club(motorcycle_graph, rnd);

  motorcycle_graph.construct_motorcycle_graph();

  return 0;
}
