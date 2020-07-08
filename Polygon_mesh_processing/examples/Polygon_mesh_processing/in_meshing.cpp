#include <core/util.h>
#include <core/output.h>
#include <core/mesh_completion.h>
#include <core/adjust_geometry.h>

#include <core/dualcontouring/connectivity.h>

#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>

typedef CGAL::Simple_cartesian<double>                                        Kernel;
typedef Kernel::Point_3                                                       Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>                                   Surface_mesh;
typedef CGAL::dynamic_face_property_t<bool>                                   Face_bool_tag;
typedef typename boost::property_map<Surface_mesh, Face_bool_tag>::type   Mark_map;

//<editor-fold desc="Ajouts pour CGAL, rien à voir avec le code original">
//<editor-fold desc="pour exporter un petit bout">
std::string extract(const std::string& str, char separator, unsigned& i)
{
  std::string res;
  while(str[i] != separator && str[i] != ')')
  {
    res += str[i];
    ++i;
  }
  ++i;
  return res;
}

std::vector<Eigen::Vector3f> load_boundary()
{
  std::string boundary_file = "/home/felix/Bureau/Geo_Facto/PSR/tests-code/jeux-de-test/tests-bord-unique/bord.wkt";
  std::ifstream is(boundary_file);

  std::vector<Eigen::Vector3f> boundary;
  std::string data;
  std::getline(is, data);
  unsigned i = 0;

  while(data[i] != '(')
  {
    ++i;
  }
  i+=2;

  while(data[i] != ')')
  {
    double x, y, z;
    x = std::stod(extract(data, ' ', i));
    y = std::stod(extract(data, ' ', i));
    z = std::stod(extract(data, ',', i));
    boundary.emplace_back(x, y, z);
  }

  return boundary;
}

std::vector<unsigned> make_boundary_indices(const std::vector<Eigen::Vector3f>& boundary, const point_mesh& mesh)
{
  std::vector<unsigned> indices;
  for(auto& v : boundary)
  {
    double min = 100;
    for(unsigned i = 0; i < mesh.vertices.size(); ++i)
    {
      double a = std::abs(v[0] - mesh.vertices[i].position[0]);
      double b = std::abs(v[1] - mesh.vertices[i].position[1]);
      double c = std::abs(v[2] - mesh.vertices[i].position[2]);
      double max = std::max(a, std::max(b, c));
      min = std::min(min, max);

      if(v[0] == mesh.vertices[i].position[0]
      && v[1] == mesh.vertices[i].position[1]
      && v[2] == mesh.vertices[i].position[2])
      {
        indices.push_back(i);
        break;
      }
    }
    std::cout << min << ' ';
  }

  std::cout << std::endl;
  return indices;
}

std::vector<unsigned> make_boundary_indices_transformed (const std::vector<Eigen::Vector3f>& boundary, const output& out, const completed_mesh& mesh)
{
  std::vector<unsigned> indices;
  for(auto& v : boundary)
  {
    double min = 100;
    for(unsigned i = 0; i < mesh.vertices.size(); ++i)
    {
      double a = std::abs(v[0] - from_hpos<float>(out.back_transform * to_hpos(mesh.vertices[i].position))[0]);
      double b = std::abs(v[1] - from_hpos<float>(out.back_transform * to_hpos(mesh.vertices[i].position))[1]);
      double c = std::abs(v[2] - from_hpos<float>(out.back_transform * to_hpos(mesh.vertices[i].position))[2]);
      double max = std::max(a, std::max(b, c));
      min = std::min(min, max);

      if(std::abs(v[0] - from_hpos<float>(out.back_transform * to_hpos(mesh.vertices[i].position))[0]) < 0.0000008
         && std::abs(v[1] - from_hpos<float>(out.back_transform * to_hpos(mesh.vertices[i].position))[1]) < 0.0000008
         && std::abs(v[2] - from_hpos<float>(out.back_transform * to_hpos(mesh.vertices[i].position))[2]) < 0.0000008)
      {
        indices.push_back(i);
        break;
      }
    }
    std::cout << min << ' ';
  }

  std::cout << std::endl;
  return indices;
}

void dump_mesh(const point_mesh& mesh, const std::string& filename = "dump")
{
  std::ofstream os("/home/felix/Bureau/Geo_Facto/PSR/tests-code/dumps/" + filename);
  for(auto& v : mesh.vertices)
  {
    os << v.position[0] << ' ' << v.position[1] << ' ' << v.position[2] << ' ';
    os << v.normal[0] << ' ' << v.normal[1] << ' ' << v.normal[2] << std::endl;
  }
}

void save_interesting_part(const std::vector<unsigned>& boundary_indices, const output& out, const completed_mesh& mesh)
{
  Surface_mesh sm;

  std::vector<Point_3> points;
  for(auto& v : mesh.vertices)
  {
    double a = from_hpos<float>(out.back_transform * to_hpos(v.position))[0];
    double b = from_hpos<float>(out.back_transform * to_hpos(v.position))[1];
    double c = from_hpos<float>(out.back_transform * to_hpos(v.position))[2];
    points.emplace_back(a, b, c);
  }

  std::vector<std::vector<std::size_t> > polygons;
  for(unsigned i = 0; i < mesh.indices.size(); ++i)
  {
    if(i % 3 == 0)
    {
      polygons.emplace_back();
    }
    polygons.back().push_back(mesh.indices[i]);
  }

  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, sm);

//  Surface_mesh::Vertex_index v = *(S.vertices().begin() + boundary_indices[0]);
//  Surface_mesh::Halfedge_index h = S.halfedge(v);
//  S.hal

  Mark_map to_delete = get(Face_bool_tag(), sm);
  std::vector<Surface_mesh::Face_index> marked;

  for(unsigned i = 0; i < boundary_indices.size() - 1; ++i)
  {
    Surface_mesh::Vertex_index v1 = *(sm.vertices().begin() + boundary_indices[i]);
    Surface_mesh::Vertex_index v2 = *(sm.vertices().begin() + boundary_indices[i + 1]);
    Surface_mesh::Halfedge_index h = sm.halfedge(v1, v2);
    put(to_delete, sm.face(h), true);
    put(to_delete, sm.face(sm.opposite(sm.next(h))), true);
    marked.push_back(sm.face(sm.opposite(sm.next(h))));
  }

  while(!marked.empty())
  {
    Surface_mesh::Face_index f = marked.back();
    marked.pop_back();
    Surface_mesh::Halfedge_index h = sm.opposite(sm.halfedge(f));
    for(unsigned i = 0; i < 3; ++i, h = sm.opposite(sm.next(sm.opposite(h))))
    {
      if(get(to_delete, sm.face(h)))
      {
        continue;
      }
      marked.push_back(sm.face(h));
      put(to_delete, sm.face(h), true);
    }
  }

  for(auto& f : sm.faces())
  {
    if(get(to_delete, f))
    {
      CGAL::Euler::remove_face(halfedge(f, sm), sm);
    }
  }


  std::ofstream os("/home/felix/Bureau/Geo_Facto/PSR/tests-code/bord-unique/reconstruction.off");
  CGAL::write_off(os, sm);

}
//</editor-fold>
//<editor-fold desc="pour choisir gérer seulement un petit bout">

//</editor-fold>
//</editor-fold>

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

int main(int argc, char** argv)
{
  /// 0) Parameters + loading
  /// =============================================================================================================================================
  std::cout.setf(std::ios::unitbuf);
  const parameters cmd(argc, argv);
  cmd.print();
  std::cout << "All timings in milliseconds\n";

  output out;
  out.enable_export = cmd.out_enable_export;
  out.enable_log_data = cmd.out_enable_log_data;
  out.enable_log_quality = cmd.out_enable_log_quality;
  out.enable_log_timing = cmd.out_enable_log_timing;

  /// Load & transform points
  out.start_timing();
  point_mesh mesh = point_mesh::load(cmd.filename);

  dump_mesh(mesh, "mesh");

  out.back_transform = mesh.transform_to_unit(1.25f);
  out.stop_timing("Loading & transforming points");
  if(mesh.vertices.empty())
  {
    std::cout << "Empty mesh !\n";
    exit(0);
  }
  if(cmd.allow_recompute_normals)
  {
    size_t invalid_normals = 0;
    for(size_t v = 0; v < mesh.vertices.size(); v++)
      if(mesh.vertices[v].normal.squaredNorm() <= 0.1)
        invalid_normals++;
    if(invalid_normals > 0)
    {
      printf("%i invalid normals => recompute them\n", (int)invalid_normals);
      const std::vector<Eigen::Vector3f> normals = compute_normals(mesh.vertices, mesh.indices);
      for(size_t v = 0; v < mesh.vertices.size(); v++)
        mesh.vertices[v].normal = normals[v];
    }
  }
  if(out.enable_export)
    export_faceless_vertices(mesh, 0.05f, "0-control-points.ply", out.back_transform);
  out.log_metric("Input/Vertices", mesh.vertices.size());
  out.log_metric("Input/Triangles", mesh.indices.size() / 3);
  out.log_metric("Input/Depth", cmd.max_depth);

  /// 1) Mesh completion
  /// =============================================================================================================================================
  std::cout << "\n1. Mesh completion\n===============================================================================\n";
  completed_mesh out_mesh = std::move(mesh);

  if(!cmd.input_already_completed())
    mesh_in_filling(out_mesh, out, cmd.max_depth, cmd.export_raw_poisson_surface);
  else
  {
    std::cout << "Skipped\n";
    out_mesh.num_base_vertices = cmd.in_base_vertices;
    out_mesh.num_base_triangles = cmd.in_base_triangles;
  }

  /// 2) Large scale geometric adjustments
  /// =============================================================================================================================================
  std::cout << "\n2. Large scale geometry adjustments\n===============================================================================\n";
  out.log_metric("Total vertices", out_mesh.vertices.size());
  out.log_metric("Total triangles", out_mesh.indices.size() / 3);
  out.log_metric("Base vertices", out_mesh.num_base_vertices);
  out.log_metric("Base triangles", out_mesh.num_base_triangles);
  geometry_new(out, out_mesh, cmd.allow_same_boundary_chains, cmd.salient_angle_deg);

  /// 3) Final output
  /// =============================================================================================================================================
  out.enable_export = true;
  out.log_metric("Output/Vertices", out_mesh.vertices.size());
  out.log_metric("Output/Triangles", out_mesh.indices.size() / 3);
  out.save_object("Final object", cmd.output_filename(), out_mesh);

  //reconstruction qui travaille sur la surface entière
//  std::vector<Eigen::Vector3f> boundary = load_boundary();
//  std::vector<unsigned> boundary_indices2 = make_boundary_indices_transformed(boundary, out, out_mesh);
//  save_interesting_part(boundary_indices2, out, out_mesh);


}
