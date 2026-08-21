#define CGAL_AW3_PROFILING

#include "output_helper.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Real_timer.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <iostream>
#include <string>

namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point_3>;

using Field = std::function<double(double,double,double)>;

// vertex_property and vertex_descriptor to store and query curvature
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
using VP = CGAL::dynamic_vertex_property_t<K::FT>;
using VPC = CGAL::dynamic_vertex_property_t<PMP::Principal_curvatures_and_directions<K>>;
using VP_PM = boost::property_map<Mesh, VP>::type;
using VPC_PM = boost::property_map<Mesh, VPC>::type;

// AABB_tree for closest point
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits_3<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

double diag_length(Mesh mesh) {
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  return std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                   CGAL::square(bbox.ymax() - bbox.ymin()) +
                   CGAL::square(bbox.zmax() - bbox.zmin()));
}

Field dist_to_center_field_generator (Mesh mesh) {
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  const double half_diag_length = diag_length(mesh) / 2.;
  Point_3 center = {
    (bbox.xmax()+bbox.xmin())/2.,
    (bbox.ymax()+bbox.ymin())/2.,
    (bbox.zmax()+bbox.zmin())/2.
  };
  double range = 1./8.;
  double start = 1./100.;
  Field field = [=](double x, double y, double z) {
    Point_3 point = {x,y,z};
    double dist = std::sqrt(squared_distance(point, center));
    double normalized = dist / half_diag_length;
    double alpha = half_diag_length * (start + range * normalized);
    return alpha;
  };
  return field;
}

Field dist_to_bbox_field_generator (Mesh mesh) {
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  const double max_length = std::max(bbox.xmax() - bbox.xmin(),
                                     std::max(bbox.ymax() - bbox.ymin(),
                                              bbox.zmax() - bbox.zmin()));
  Point_3 center = {
    (bbox.xmax()+bbox.xmin())/2.,
    (bbox.ymax()+bbox.ymin())/2.,
    (bbox.zmax()+bbox.zmin())/2.
  };
  double range = 1./8.;
  double start = 1./100.;
  Field field = [=](double x, double y, double z) {
    Point_3 point = {x-center.x(),y-center.y(),z-center.z()};
    double qx = std::abs(x) - bbox.xmax();
    double qy = std::abs(y) - bbox.ymax();
    double qz = std::abs(z) - bbox.zmax();
    double dist = - std::min(std::max(qx,std::max(qy,qz)),0.0);
	double normalized = dist / max_length;
	double alpha = max_length * (start + range * normalized);
	return alpha;
  };
  return field;
}

Field xaxis_field_generator(Mesh mesh) {  
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  const double x_length = bbox.xmax() - bbox.xmin();
  double range = 1./8.;
  double start = 1./100.;
  Field field = [=](double x, double y, double z) {
    double centered = x - bbox.xmin();
    double clamped = std::max(0.,centered);
    double normalized = clamped / x_length;
    double alpha = x_length * (start + range * normalized);
    return alpha;
  };
  return field;
}

Field curvature_field_generator(Mesh& mesh) {

  auto tree = std::make_shared<Tree>(mesh.faces_begin(), mesh.faces_end(), mesh);

  // define property map to store curvature value and directions
  auto mean_curvature_map = std::make_shared<VP_PM>(get(VP(), mesh));
  auto Gaussian_curvature_map = std::make_shared<VP_PM>(get(VP(), mesh));
  auto principal_curvatures_and_directions_map = std::make_shared<VPC_PM>(get(VPC(), mesh));
  // Compute curvatures
  PMP::interpolated_corrected_curvatures(mesh,
		CGAL::parameters::vertex_mean_curvature_map(*mean_curvature_map)
			.vertex_Gaussian_curvature_map(*Gaussian_curvature_map)
			.vertex_principal_curvatures_and_directions_map(*principal_curvatures_and_directions_map));

  Field field = [tree, &mesh, principal_curvatures_and_directions_map] (double x, double y, double z) {
	Point_3 point = {x,y,z};
	face_descriptor fd = tree->closest_point_and_primitive(point).second;
	Mesh::Halfedge_index h = halfedge(fd, mesh); //
	vertex_descriptor vd = target(h, mesh);
	auto PC = get(*principal_curvatures_and_directions_map, vd);
	double max_curvature = std::max(std::abs(PC.min_curvature), std::abs(PC.max_curvature));
	double alpha = 1. / max_curvature;
	return alpha;
  };
  return field;
}

std::function<Field(Mesh&)> field_generator = curvature_field_generator;

int main(int argc, char** argv)
{
  // Read the input
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/armadillo.off");
  std::cout << "Reading " << filename << "..." << std::endl;

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh) || is_empty(mesh) || !is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input:" << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input: " << num_vertices(mesh) << " vertices, " << num_faces(mesh) << " faces" << std::endl;

  // Compute the alpha and offset values
  const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;

  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));

  Field alpha_field = field_generator(mesh);
  const double offset = diag_length / relative_offset;
  std::cout << "offset: " << offset << std::endl;

  // Construct the wrap
  CGAL::Real_timer t;
  t.start();

  Mesh wrap;
  CGAL::alpha_wrap_3(mesh, alpha_field, offset, wrap);

  t.stop();
  std::cout << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces" << std::endl;
  std::cout << "Took " << t.time() << " s." << std::endl;

  // Save the result
  const std::string output_name = generate_output_name(filename, 0, relative_offset); // FIXME
  std::cout << "Writing to " << output_name << std::endl;
  CGAL::IO::write_polygon_mesh(output_name, wrap, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
