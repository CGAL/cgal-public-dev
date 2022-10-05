#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Curvatures/interpolated_corrected_curvature_measures.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_Kernel;
typedef CGAL::Surface_mesh<Epic_Kernel::Point_3> Surface_Mesh;
typedef boost::graph_traits<Surface_Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Surface_Mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[])
{
  Surface_Mesh g1;
  const std::string filename = (argc>1) ?
      argv[1] :
      CGAL::data_file_path("meshes/small_bunny.obj");

  if(!CGAL::IO::read_polygon_mesh(filename, g1))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  PMP::Interpolated_corrected_curvatures_computer<Surface_Mesh> icc(g1, true, true, true);

  icc.compute_all_curvatures();

  // uncomment this to compute a curvature while specifying named parameters
  // Example: an expansion ball radius of 0.5 and a vertex normals map (does not have to depend on positions)
  
  /*Surface_Mesh::Property_map<vertex_descriptor, Epic_Kernel::Vector_3> vnm;
  boost::tie(vnm, created) = g1.add_property_map<vertex_descriptor, Epic_Kernel::Vector_3>(
      "v:vnm", Epic_Kernel::Vector_3(0, 0, 0)
  );

  assert(created);

  PMP::interpolated_corrected_mean_curvature(
      g1,
      mean_curvature_map,
      CGAL::parameters::ball_radius(0.5).vertex_normal_map(vnm)
  );*/


  for (vertex_descriptor v : vertices(g1))
  {
    const auto& PC = icc.principal_curvature_map[v];
      std::cout << v.idx() << ": HC = " << icc.mean_curvature_map[v]
                           << ", GC = " << icc.gaussian_curvature_map[v] << "\n"
                           << ", PC = [ " << PC.min_curvature << " , " << PC.max_curvature << " ]\n";
  }
}
