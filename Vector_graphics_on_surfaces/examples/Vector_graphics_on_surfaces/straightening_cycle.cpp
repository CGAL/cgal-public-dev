#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Vector_graphics_on_surfaces/locally_shortest_path.h>
#include <CGAL/Vector_graphics_on_surfaces/utils.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

namespace VGoS = CGAL::Vector_graphics_on_surfaces;
namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;

using Edge_location = PMP::Edge_location<Mesh, double>;

int main()
{
  std::string filename = CGAL::data_file_path("meshes/elephant.off");

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  std::vector<Mesh::Vertex_index> cycle;
  cycle.emplace_back(1093);
  cycle.emplace_back(712);
  cycle.emplace_back(2588);
  cycle.emplace_back(843);
  cycle.emplace_back(2577);
  cycle.emplace_back(1122);
  cycle.emplace_back(1123);

  std::ofstream debug("inital_path.polylines.txt");
  debug << cycle.size()+1;
  for (auto v : cycle)
    debug << " " << mesh.point(v);
  debug << " " << mesh.point(cycle[0]) << "\n";
  debug.close();

  std::vector<Edge_location> edge_locations;
  VGoS::straighten_cycle<double>(cycle[2], cycle, mesh, edge_locations);

  std::vector<K::Point_3> poly;
  poly.reserve(edge_locations.size()+2);
  auto loc = PMP::locate_vertex<double>(cycle[2], mesh);
  VGoS::convert_path_to_polyline(loc , edge_locations, loc, mesh, std::back_inserter(poly));

  std::ofstream out("straightened_cycle.polylines.txt");
  out << poly.size();
  for (auto p : poly)
    out << " " << p;
  out << "\n";
  out.close();

  std::vector<Edge_location> edge_locations_bis;
  VGoS::straighten_cycle<double>(0, edge_locations, mesh, edge_locations_bis);

  std::vector<K::Point_3> poly_bis;
  poly_bis.reserve(edge_locations_bis.size()+2);
  auto loc_bis = PMP::to_face_location(edge_locations[0], mesh);
  VGoS::convert_path_to_polyline(loc_bis , edge_locations_bis, loc_bis, mesh, std::back_inserter(poly_bis));

  out.open("straightened_cycle_bis.polylines.txt");
  out << poly_bis.size();
  for (auto p : poly_bis)
    out << " " << p;
  out << "\n";



  return 0;
}
