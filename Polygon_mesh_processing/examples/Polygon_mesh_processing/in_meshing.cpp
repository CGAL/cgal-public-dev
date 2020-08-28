//includes in-meshing
#include <core/util.h>
#include <core/output.h>
#include <core/mesh_completion.h>
#include <core/adjust_geometry.h>
#include <core/dualcontouring/connectivity.h>

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

namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::Simple_cartesian<double>                                        Kernel;
typedef Kernel::Point_3                                                       Point_3;
typedef Kernel::Vector_3                                                      Vector_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>                                   Surface_mesh;
typedef CGAL::dynamic_face_property_t<bool>                                   Face_bool_tag;
typedef typename boost::property_map<Surface_mesh, Face_bool_tag>::type       Mark_map;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor                  vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor                    face_descriptor;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor                halfedge_descriptor;
typedef std::vector<vertex_descriptor>                                        vertex_range;
typedef std::vector<face_descriptor>                                          face_range;
typedef std::vector<halfedge_descriptor>                                      halfedge_range;
typedef CGAL::Face_filtered_graph<Surface_mesh>                               Filtered_graph;

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

//<editor-fold desc="Version 2.0">
vertex_range make_hole_points(const Surface_mesh& sm, const std::string& holes_file)
{
  std::cout << "Acquisition des points qui bordent les trous\n";

  vertex_range hole_points;
  std::ifstream is(holes_file);
  double x, y, z;

  unsigned i = 0;
  std::string line;
  while(line != "end_header")
  {
    is >> line;
  }
  is >> x >> y >> z;
  while(!is.eof())
  {
    Point_3 p(x, y, z);
    for(auto& v : sm.vertices())
    {
      if(sm.point(v) == p)
      {
        hole_points.push_back(v);
      }
    }
    ++i;
    is >> x >> y >> z;
  }
  std::cout << i << " points lus, " << hole_points.size() << " points trouvés sur la surface, "
  << hole_points.size() - i << " ratés\n\n";

  return hole_points;
}

face_range make_first_ring(const Surface_mesh& sm, const vertex_range& hole_points, std::vector<std::vector<unsigned>>& holes_indices)
{
  std::set<Surface_mesh::Face_index> added_faces;
  face_range first_ring;

  for(const auto& v : hole_points)
  {
    //récupérer le halfedge h qui borde le trou
    Surface_mesh::Halfedge_index h;
    for(auto& h_ : sm.halfedges_around_target(sm.halfedge(v)))
    {
      if(sm.face(h_) == Surface_mesh::null_face())
      {
        h = h_;
        break;
      }
    }

    holes_indices.emplace_back();

    //se balader autour du trou
    Surface_mesh::Halfedge_index it = h;
    do{
      holes_indices.back().push_back(it.idx());

      for(const auto& f : sm.faces_around_target(it))
      {
        if(f == Surface_mesh::null_face() || added_faces.count(f) > 0)
        {
          continue;
        }
        first_ring.push_back(f);
        added_faces.insert(f);
      }
      it = sm.next(it);
    }while(it != h);
  }

  return first_ring;
}

point_mesh make_point_mesh(const Surface_mesh& sm, const face_range& rings,
    std::vector<std::vector<unsigned>>& sm_holes_indices, std::vector<std::vector<unsigned>>& pm_holes_indices)
{
  // sm_holes_indices[i] contient la liste des indices des halfedges qui bordent le ieme trou

  // map from the vertices indices of sm to the ones of mesh
  std::map<unsigned, unsigned> sm_indices_to_pm_indices;

  unsigned i = 0;
  std::vector<oriented_point> points;
  std::set<Surface_mesh::Vertex_index> treated;
  for(const auto& f : rings)
  {
    for(const auto& v : sm.vertices_around_face(sm.halfedge(f)))
    {
      if(treated.count(v) > 0)
      {
        continue;
      }

      Point_3 p = sm.point(v);
      auto n = CGAL::Polygon_mesh_processing::compute_vertex_normal(v, sm);
      Eigen::Vector3f p_(p.x(), p.y(), p.z());
      Eigen::Vector3f n_(n.x(), n.y(), n.z());

      points.emplace_back(p_, n_);
      sm_indices_to_pm_indices[v.idx()] = i;
      ++i;
      treated.insert(v);
    }
  }

  std::vector<unsigned> faces;
  for(const auto& f : rings)
  {
    for(const auto& v : sm.vertices_around_face(sm.halfedge(f)))
    {
      faces.push_back(sm_indices_to_pm_indices[v.idx()]);
    }
  }

  for(unsigned j = 0; j < sm_holes_indices.size(); ++j)
  {
    pm_holes_indices.emplace_back();
    for(unsigned k = 0; k < sm_holes_indices[j].size(); ++k)
    {
      Surface_mesh::Halfedge_index h = *(sm.halfedges().begin() + sm_holes_indices[j][k]);
      Surface_mesh::Vertex_index v = sm.source(h);
      pm_holes_indices.back().push_back(sm_indices_to_pm_indices[v.idx()]);
    }
  }

  return point_mesh(points, faces);
}

point_mesh make_point_mesh_for_in_meshing(const Surface_mesh& sm, const std::string& holes_file, const std::string& guide_file,
                                          std::vector<std::vector<unsigned>>& sm_holes_indices, std::vector<std::vector<unsigned>>& pm_holes_indices,
                                          unsigned expand_degree = 2)
{
  std::ofstream os;
  //récupérer les trous
  vertex_range hole_points = make_hole_points(sm, holes_file);

  //récupérer la première couronne
  face_range rings = make_first_ring(sm, hole_points, sm_holes_indices);
  //maintenant sm_hole_indices[i] = [id0, id2, ... idn] où id0, ... idn sont les indices dans sm des halfedges qui forment le i-eme trou, et id0 = idn

  //expand
  typedef boost::graph_traits<Surface_mesh>::face_descriptor     face_descriptor;
  auto selected = get(Face_bool_tag(), sm);
  for(auto& f : rings)
  {
    put(selected, f, true);
  }
  CGAL::expand_face_selection(rings, sm, expand_degree, selected, std::back_inserter(rings));

  Filtered_graph ffg_second_rings(sm, rings);
  Surface_mesh second_rings_sm;
  CGAL::copy_face_graph(ffg_second_rings, second_rings_sm);
  os = std::ofstream("/home/felix/Bureau/Geo_Facto/PSR/tests-code/couronnes/dump_second_rings.off");
  CGAL::write_off(os, second_rings_sm);
  os.close();


  //faire le point_mesh
  point_mesh mesh = make_point_mesh(sm, rings, sm_holes_indices,pm_holes_indices);

  //ajouter le guide
  if(!guide_file.empty())
  {
    mesh.add_guide(guide_file);
  }

  return mesh;
}

Surface_mesh make_surface_mesh_from_completed_mesh(const output& out, const completed_mesh& out_mesh)
{
  Surface_mesh sm;

  std::vector<Point_3> points;
  for(auto& v : out_mesh.vertices)
  {
    double x = from_hpos<float>(out.back_transform * to_hpos(v.position))[0];
    double y = from_hpos<float>(out.back_transform * to_hpos(v.position))[1];
    double z = from_hpos<float>(out.back_transform * to_hpos(v.position))[2];
    points.emplace_back(x, y, z);
  }

  std::vector<std::vector<std::size_t>> polygons;
  for(unsigned i = 0; i < out_mesh.indices.size(); ++i)
  {
    if(i % 3 == 0)
    {
      polygons.emplace_back();
    }
    polygons.back().push_back(out_mesh.indices[i]);
  }

  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, sm);

  return sm;
}

Surface_mesh extract_reconstruction(const output& out, const completed_mesh& out_mesh,
                                    const std::vector<std::vector<unsigned>>& pm_holes_indices,
                                    std::vector<std::vector<Surface_mesh::Halfedge_index>>& reconstruction_holes,
                                    const std::string& save_file)
{
  Surface_mesh reconstruction = make_surface_mesh_from_completed_mesh(out, out_mesh);
  std::ofstream os;

  Mark_map to_delete = get(Face_bool_tag(), reconstruction);
  std::vector<Surface_mesh::Face_index> to_handle;

  for(auto& hole : pm_holes_indices)
  {
    reconstruction_holes.emplace_back();
    unsigned size = hole.size();
    for(unsigned i = 0; i < size; ++i)
    {
      Surface_mesh::Vertex_index v1 = *(reconstruction.vertices().begin() + hole[i]);
      Surface_mesh::Vertex_index v2 = *(reconstruction.vertices().begin() + hole[(i + 1) % size]);
      Surface_mesh::Halfedge_index h = reconstruction.halfedge(v1, v2);
      reconstruction_holes.back().push_back(h);
      put(to_delete, reconstruction.face(reconstruction.opposite(h)), true);
    }
  }

  for(auto& hole : pm_holes_indices)
  {
    Surface_mesh::Vertex_index v1 = *(reconstruction.vertices().begin() + hole[0]);
    Surface_mesh::Vertex_index v2 = *(reconstruction.vertices().begin() + hole[1]);
    Surface_mesh::Halfedge_index h = reconstruction.halfedge(v1, v2);
    Surface_mesh::Halfedge_index h_ = reconstruction.opposite(reconstruction.next(reconstruction.opposite(h)));
    if(!get(to_delete, reconstruction.face(h_)))
    {
      put(to_delete, reconstruction.face(h_), true);
      to_handle.push_back(reconstruction.face(h_));
    }
    else
    {
      put(to_delete, reconstruction.face(reconstruction.opposite(reconstruction.next(reconstruction.opposite(h_)))), true);
      to_handle.push_back(reconstruction.face(reconstruction.opposite(reconstruction.next(reconstruction.opposite(h_)))));
    }
  }

  while(!to_handle.empty())
  {
    Surface_mesh::Face_index f = to_handle.back();
    to_handle.pop_back();
    Surface_mesh::Halfedge_index h = reconstruction.opposite(reconstruction.halfedge(f));
    for(unsigned i = 0; i < 3; ++i, h = reconstruction.opposite(reconstruction.next(reconstruction.opposite(h))))
    {
      if(get(to_delete, reconstruction.face(h)))
      {
        continue;
      }
      to_handle.push_back(reconstruction.face(h));
      put(to_delete, reconstruction.face(h), true);
    }
  }

  for(auto& f : reconstruction.faces())
  {
    if(get(to_delete, f))
    {
      CGAL::Euler::remove_face(halfedge(f, reconstruction), reconstruction);
    }
  }

  os = std::ofstream(save_file);
  CGAL::write_off(os, reconstruction);
  os.close();

  return reconstruction;
}

typedef std::pair<Surface_mesh::Halfedge_index, Surface_mesh::Halfedge_index> Halfedge_pair;
std::vector<Halfedge_pair> make_associations(const Surface_mesh& sm, const Surface_mesh& reconstruction,
                                             const std::vector<std::vector<unsigned>>& sm_holes_indices,
                                             const std::vector<std::vector<Surface_mesh::Halfedge_index>>& reconstruction_holes,
                                             const std::vector<Halfedge_pair>& h2h)
{
  std::vector<Halfedge_pair> associations;
  for(unsigned i = 0; i < sm_holes_indices.size(); ++i)
  {
    unsigned k = 0;
    for(unsigned j = 0; j < sm_holes_indices[i].size(); ++j)
    {
      Surface_mesh::Halfedge_index h1 = *(sm.halfedges().begin() + sm_holes_indices[i][j]);
      Surface_mesh::Halfedge_index h2 = reconstruction_holes[i][j];
      for(auto& it : h2h)
      {
        if(it.first == h2)
        {
          ++k;
          h2 = sm.opposite(it.second);
          break;
        }
      }
      associations.emplace_back(h1, h2);
    }
    std::cout << "trouvés : " << k << '/' << sm_holes_indices[i].size() - 1 << std::endl;
  }
  return associations;
}
//</editor-fold>


//<editor-fold desc="Normales de Mael pour la prothèse dentaire">
point_mesh make_mesh()
{
  std::string tobehealed2_file = "/home/felix/Bureau/Geo_Facto/PSR/tests-SGP13/jeux-de-test/RedCircleBox/RedCircleBoxOneCC.off";
  std::ifstream is(tobehealed2_file);
  Surface_mesh sm;
  CGAL::read_off(is, sm);

  auto vnormals = sm.add_property_map<vertex_descriptor, Vector_3>("v:normals", CGAL::NULL_VECTOR).first;
  auto fnormals = sm.add_property_map<face_descriptor, Vector_3>("f:normals", CGAL::NULL_VECTOR).first;

  PMP::compute_normals(sm, vnormals, fnormals);

  std::map<unsigned, unsigned> sm_indices_to_pm_indices;

  unsigned i = 0;
  std::vector<oriented_point> points;
  std::set<Surface_mesh::Vertex_index> treated;
  for(const auto& f : faces(sm))
  {
    for(const auto& v : sm.vertices_around_face(sm.halfedge(f)))
    {
      if(treated.count(v) > 0)
      {
        continue;
      }

      Point_3 p = sm.point(v);
      auto n = vnormals[v];
      Eigen::Vector3f p_(p.x(), p.y(), p.z());
      Eigen::Vector3f n_(n.x(), n.y(), n.z());

      points.emplace_back(p_, n_);
      sm_indices_to_pm_indices[v.idx()] = i;
      ++i;
      treated.insert(v);
    }
  }

  std::vector<unsigned> indices;
  for(const auto& f : faces(sm))
  {
    for(const auto& v : sm.vertices_around_face(sm.halfedge(f)))
    {
      indices.push_back(sm_indices_to_pm_indices[v.idx()]);
    }
  }

  return point_mesh(points, indices);
}

//</editor-fold>

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

  std::string mesh_file = "/home/felix/Bureau/Geo_Facto/PSR/tests-code/jeux-de-test/tests-couronnes/test4/cube-deux-trous.off";
  std::string holes_file = "/home/felix/Bureau/Geo_Facto/PSR/tests-code/jeux-de-test/tests-couronnes/test4/trous.ply";
  std::string guide_file = "";
  Surface_mesh sm;
  std::ifstream is(mesh_file);
  CGAL::read_off(is, sm);

  std::vector<std::vector<unsigned>> sm_holes_indices;
  std::vector<std::vector<unsigned>> pm_holes_indices;
  point_mesh mesh = make_point_mesh_for_in_meshing(sm, holes_file, guide_file, sm_holes_indices, pm_holes_indices, 6);

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

  //extrction de la partie intéressante

  std::vector<std::vector<Surface_mesh::Halfedge_index>> reconstruction_holes;
  std::string reconstruction_file = "/home/felix/Bureau/Geo_Facto/PSR/tests-code/couronnes/reconstruction.off";
  Surface_mesh reconstruction = extract_reconstruction(out, out_mesh, pm_holes_indices, reconstruction_holes, reconstruction_file);
  std::vector<Halfedge_pair> h2h;
  CGAL::copy_face_graph(reconstruction, sm, CGAL::parameters::halfedge_to_halfedge_output_iterator(std::back_inserter(h2h)));
  std::vector<Halfedge_pair> associations = make_associations(sm, reconstruction, sm_holes_indices, reconstruction_holes, h2h);


  CGAL::Polygon_mesh_processing::stitch_borders(sm, associations);


  std::string save_file = "/home/felix/Bureau/Geo_Facto/PSR/tests-code/couronnes/cube-deux-trous.off";
  std::ofstream os(save_file);
  CGAL::write_off(os, sm);
  os.close();
}
