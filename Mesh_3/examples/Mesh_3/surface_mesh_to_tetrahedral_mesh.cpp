#include <CGAL/config.h>
#include <CGAL/Default.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/io.h>
#include <CGAL/IO/output_to_vtu.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/iterator.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Sizing_field_with_aabb_tree.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/tags.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <utility>

// Domain
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Surface_mesh = CGAL::Surface_mesh<K::Point_3>;
using Mesh_domain = CGAL::Polyhedral_mesh_domain_with_features_3<K, Surface_mesh>;


#ifdef CGAL_CONCURRENT_MESH_3
using Concurrency_tag = CGAL::Parallel_tag;
#else
using Concurrency_tag = CGAL::Sequential_tag;
#endif

// Triangulation
using Tr = CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type;

using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index>;

// Criteria
using Features_sizing_field = CGAL::Sizing_field_with_aabb_tree<K, Mesh_domain>;
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

namespace params = CGAL::parameters;

int main(int argc, char*argv[])
{
  const std::string input_file_name = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/fandisk.off");
  const std::string output_file_name = (argc > 2) ? argv[2] : "out-tetmesh.vtu";
  Surface_mesh surface_mesh;
  if(!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(input_file_name, surface_mesh)) {
    std::cerr << "Error: Cannot read file " << input_file_name << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cout << "Read surface mesh from file " << input_file_name << std::endl;
    std::cout << "  - number of vertices: " << num_vertices(surface_mesh) << std::endl;
    std::cout << "  - number of faces:    " << num_faces(surface_mesh) << std::endl;
  }

  if(!CGAL::is_triangle_mesh(surface_mesh)) {
    std::cerr << "ERROR: the input surface mesh is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  const bool is_closed = CGAL::is_closed(surface_mesh);
  if(is_closed) {
    std::cout << "  - the input surface mesh is closed." << std::endl;
  } else {
    std::cout << "  - WARNING: the input surface mesh is not closed." << std::endl;
  }

  // check if the input surface mesh is manifold
  std::size_t nb_of_non_manifold_vertices{0};
  CGAL::Counting_output_iterator count_it{&nb_of_non_manifold_vertices};
  count_it = CGAL::Polygon_mesh_processing::non_manifold_vertices(surface_mesh, count_it);
  if(nb_of_non_manifold_vertices > 0) {
    std::cout << "  - WARNING: the input surface mesh has " << nb_of_non_manifold_vertices
              << " non-manifold vertices."<< std::endl;
  }
  const bool self_intersecting = CGAL::Polygon_mesh_processing::does_self_intersect(surface_mesh);
  if(self_intersecting) {
    std::cout << "  - WARNING: the input surface mesh is self-intersecting." << std::endl;
  }

  const bool manifold = !self_intersecting && nb_of_non_manifold_vertices == 0;

  if(manifold && is_closed) {
    std::cout << "  - the input surface mesh is manifold." << std::endl;
  }

  // Create domain
  Mesh_domain domain(std::move(surface_mesh));

  // Get sharp features
  domain.detect_features();

  // Mesh criteria
  auto bbox = domain.bbox();
  const double diagonal_length =
      CGAL::sqrt(CGAL::square(bbox.x_span()) + CGAL::square(bbox.y_span()) + CGAL::square(bbox.z_span()));
  const double size_bound = diagonal_length * 0.05;
  const Features_sizing_field edge_sizing_field(0.07, domain);
  const Mesh_criteria criteria(params::edge_size(edge_sizing_field)
                                      .facet_angle(25)
                                      .facet_size(size_bound)
                                      .facet_distance(size_bound * 0.1)
                                      .cell_radius_edge_ratio(3)
                                      .cell_size(size_bound));
  const auto manifold_criteria = manifold ? (is_closed ? params::manifold()
                                                       : params::manifold_with_boundary())
                                          : (params::non_manifold());

  // Mesh generation
  std::cout << "Start meshing..." << std::endl;
  CGAL::Real_timer timer;
  timer.start();
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, manifold_criteria,
                                      params::mesh_3_options(params::nonlinear_growth_of_balls = true));
  timer.stop();

  std::cout << "Meshing completed in " << timer.time() << " seconds." << std::endl;
  std::cout << "  number of vertices:       " << c3t3.triangulation().number_of_vertices() << std::endl;
  std::cout << "  number of surface facets: " << c3t3.number_of_facets() << std::endl;
  std::cout << "  number of cells:          " << c3t3.number_of_cells() << std::endl;

  // Output
  std::ofstream ofs(output_file_name);
  ofs.precision(17);
  CGAL::IO::output_to_vtu(ofs, c3t3, CGAL::IO::ASCII);
  if(ofs.fail()) {
    std::cerr << "Error: Cannot write the surface mesh to file " << output_file_name << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cout << "Tet mesh written to file " << output_file_name << std::endl;
  }

  return EXIT_SUCCESS;
}
