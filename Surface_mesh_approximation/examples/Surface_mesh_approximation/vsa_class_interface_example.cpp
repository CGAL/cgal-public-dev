#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/vsa_mesh_approximation_traits.h>
#include <CGAL/VSA_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef Polyhedron_3::Facet_handle Facet_handle;
typedef boost::associative_property_map<std::map<Facet_handle, FT> > FacetAreaMap;

typedef CGAL::PlaneProxy<Polyhedron_3> PlaneProxy;
typedef CGAL::L2Metric<Polyhedron_3, FacetAreaMap> L2Metric;
typedef CGAL::L2ProxyFitting<Polyhedron_3> L2ProxyFitting;
typedef CGAL::VSA_approximation<Polyhedron_3, PlaneProxy, L2Metric, L2ProxyFitting> VSA;

int main(int argc, char *argv[])
{
  if (argc < 5)
    return EXIT_FAILURE;

  // read Polyhedron_3
  Polyhedron_3 mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // facet area map
  std::map<Facet_handle, FT> facet_areas;
  for (Polyhedron_3::Facet_iterator fitr = mesh.facets_begin();
    fitr != mesh.facets_end(); ++fitr) {
    Polyhedron_3::Halfedge_handle he = fitr->halfedge();
    const Point_3 &p0 = he->opposite()->vertex()->point();
    const Point_3 &p1 = he->vertex()->point();
    const Point_3 &p2 = he->next()->vertex()->point();

    FT farea(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
    facet_areas.insert(std::pair<Facet_handle, FT>(fitr, farea));
  }
  FacetAreaMap area_pmap(facet_areas);

  const std::size_t num_proxies = std::atoi(argv[3]);
  const std::size_t num_iterations = std::atoi(argv[4]);
  std::vector<int> tris;
  std::vector<Kernel::Point_3> anchor_pos;
  int init = std::atoi(argv[2]);
  if (init < 0 || init > 2)
    return EXIT_FAILURE;

  L2Metric metric(mesh, area_pmap);
  L2ProxyFitting proxy_fitting(mesh);

  // create VSA L2 metric approximation algorithm instance
  VSA l2_approx;
  l2_approx.set_mesh(mesh);
  l2_approx.set_error_metric(metric);
  l2_approx.set_proxy_fitting(proxy_fitting);

  // initialize proxies randomly on the mesh
  l2_approx.init_proxies(num_proxies, VSA::RandomInit);
  
  // run the iteration to minimize the error
  for (std::size_t i = 0; i < num_iterations; ++i)
    l2_approx.run_one_step();

  // add proxies to the one with the maximum fitting error
  l2_approx.add_proxies(VSA::IncrementalInit, 3);
  for (std::size_t i = 0; i < 5; ++i)
    l2_approx.run_one_step();

  // teleport the proxies from local minimal
  l2_approx.teleport_proxies(3, false);

  // extract the approximation polyhedron
  Polyhedron_3 out_mesh;
  l2_approx.meshing(out_mesh);

  return EXIT_SUCCESS;
}