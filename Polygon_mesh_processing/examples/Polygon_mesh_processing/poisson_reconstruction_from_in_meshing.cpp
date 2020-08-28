//includes in-meshing
#include <core/util.h>
#include <core/output.h>
#include <core/mesh_completion.h>
#include <core/adjust_geometry.h>
#include <core/dualcontouring/connectivity.h>
#include <core/dualcontouring/dual_contouring.h>
#include <poisson_recon/Bridge.h>

//includes CGAL
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/IO/read_xyz_points.h>

namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::Simple_cartesian<double>                                        Kernel;
typedef Kernel::Point_3                                                       Point_3;
typedef Kernel::Vector_3                                                      Vector_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>                                   Surface_mesh;

struct parameters
{
  const char* filename = nullptr;
  uint32_t in_base_vertices = uint32_t(-1);
  uint32_t in_base_triangles = uint32_t(-1);
  int max_depth = 8;
  bool export_raw_poisson_surface = true;
  bool out_enable_export = true;
  bool out_enable_log_quality = true;
  bool out_enable_log_timing = true;
  bool out_enable_log_data = true;
  bool allow_recompute_normals = true;
  bool allow_same_boundary_chains = false;
  float salient_angle_deg = 84.0f;

  parameters(int argc, char** argv) :
      filename(argc > 1 ? argv[1] : "../misc/caterpillar.obj")
  {
    switch(argc)
    {
      case 3:
        max_depth = std::atoi(argv[2]);
        break;
      case 4:
        in_base_vertices = std::atoi(argv[2]);
        in_base_triangles = std::atoi(argv[3]);
        break;
    }
  }

  bool input_already_completed() const { return in_base_vertices != uint32_t(-1) && in_base_triangles != uint32_t(-1); }
  std::string output_filename() const { return filename_no_ext(filename_no_dir(filename)) + "-stitched.ply"; }
  void print() const
  {
    std::cout << "Parameters:\n"
              << "\tfile: " << filename << "\n"
              << "\t------------\n"
              << "\tdepth: " << max_depth << "\n"
              << "\tbase vertices: " << in_base_vertices << "\n"
              << "\tbase triangles: " << in_base_triangles << "\n"
              << "\texport raw poisson: " << (export_raw_poisson_surface ? "yes" : "no") << "\n"
              << "\t------------\n"
              << "\tallow chains in same boundary: " << (allow_same_boundary_chains ? "yes" : "no") << "\n"
              << "\tsalient point angle (degrees): " << salient_angle_deg << " deg\n"
              << "\t------------\n"
              << "\texport intermediate objects: " << (out_enable_export ? "yes" : "no") << "\n"
              << "\tlog mesh quality: " << (out_enable_log_quality ? "yes" : "no") << "\n"
              << "\tlog timings: " << (out_enable_log_timing ? "yes" : "no") << "\n"
              << "\tlog data: " << (out_enable_log_data ? "yes" : "no") << "\n"
              << "\t------------\n"
              << "\tskip poisson: " << (input_already_completed() ? "yes" : "no") << "\n"
              << "\toutput file: " << output_filename() << "\n";
  }
};


std::vector<oriented_point> read_ply_points(const std::string& file)
{
  std::vector<oriented_point> points;

  std::ifstream is(file);
  float x, y, z, nx, ny, nz;

  std::string line;
  while(!is.eof() && line != "end_header")
  {
    is >> line;
  }
  while(!is.eof())
  {
    is >> x >> y >> z >> nx >> ny >> nz;
    points.emplace_back(Eigen::Vector3f(x, y, z), Eigen::Vector3f(nx, ny, nz));
  }

  return points;
}

std::vector<oriented_point> read_xyz_points(const std::string& file)
{
  std::vector<oriented_point> points;

  std::ifstream is(file);
  float x, y, z, nx, ny, nz;

  while(!is.eof())
  {
    is >> x >> y >> z >> nx >> ny >> nz;
    points.emplace_back(Eigen::Vector3f(x, y, z), Eigen::Vector3f(nx, ny, nz));
  }

  return points;
}

/// Return an octree that correspond to the 'space nodes' of the reference Poisson reconstruction
template<typename T>
sorted_octree<T> space_nodes(const poisson_outputs& poisson, const T& def = T())
{
  std::vector<octree_id> snodes;
  poisson.tree.for_each([&](const octree_id& node, const poisson_outputs::node&)
                        {
                          if(node == poisson.FEM_root || poisson.FEM_root.is_ancestor_of(node))
                            snodes.emplace_back(node.relative_pos(poisson.FEM_root));
                        });
  return sorted_octree<T>(snodes, def);
}

Surface_mesh point_mesh_to_surface_mesh(const point_mesh& pm, const Eigen::Matrix4f& back_transform)
{
  Surface_mesh sm;

  std::vector<Point_3> points;
  for(auto& v : pm.vertices)
  {
    double x = from_hpos<float>(back_transform * to_hpos(v.position))[0];
    double y = from_hpos<float>(back_transform * to_hpos(v.position))[1];
    double z = from_hpos<float>(back_transform * to_hpos(v.position))[2];
    points.emplace_back(x, y, z);
  }

  std::vector<std::vector<std::size_t>> polygons;
  for(unsigned i = 0; i < pm.indices.size(); ++i)
  {
    if(i % 3 == 0)
    {
      polygons.emplace_back();
    }
    polygons.back().push_back(pm.indices[i]);
  }

  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, sm);

  return sm;
}


int main(int argc, char** argv)
{
  unsigned max_depth = std::stoi(argv[3]);

  point_mesh mesh;
  std::string test_dirr = "/home/felix/Bureau/Geo_Facto/PSR/tests-poisson/jeux_de_test/";
  std::string name = argv[1];
  mesh.vertices = read_xyz_points(test_dirr + name + ".xyz");
//  mesh.vertices = read_ply_points(test-dirr + "0-1-input-points.ply");

  Eigen::Matrix4f back_transform = mesh.transform_to_unit(1.25f);

  completed_mesh out_mesh = std::move(mesh);

  // Run reference implementation for Poisson reconstruction
  poisson_outputs poisson = run_poisson({ out_mesh.vertices, max_depth, poisson_boundary_type::dirichlet, poisson_output_tree });

  // Upscale solution for faster access
  upscale_1(poisson.tree);

  sorted_octree<uint8_t> tree = space_nodes<uint8_t>(poisson);
  point_mesh m;
  dual_contouring(tree, poisson.implicit_value(), poisson.implicit_normal()).extract(m.vertices, m.indices);

  Surface_mesh sm = point_mesh_to_surface_mesh(m, back_transform);

  std::string res_dirr = argv[2];

  std::ofstream os (res_dirr + '/' +    name + std::to_string(max_depth) + ".off");
  CGAL::write_off(os, sm);
  os.close();
}
