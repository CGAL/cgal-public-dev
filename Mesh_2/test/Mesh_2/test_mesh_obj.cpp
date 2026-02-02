#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/OBJ.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_no_edge_refinement_2.h>
#include <CGAL/Triangulation_conformer_2.h>

#include <fstream>
#include <iostream>
#include <vector>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = K::FT;
using Point_2 = K::Point_2;
using Vector_2 = K::Vector_2;

using Vb = CGAL::Triangulation_vertex_base_2<K>;
using Fb = CGAL::Delaunay_mesh_face_base_2<K>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, Tds, CGAL::Exact_predicates_tag>;

int main(int argc, char*argv[] )
{
  CDT cdt;
  // std::ifstream in(argv[1]);
  std::ifstream in(argc>1?argv[1]:"mini.obj");
  if(!in)
    return false;

  std::vector<Point_2> points;
  std::vector<std::vector<std::size_t> > id_polylines;
  std::vector<std::vector<std::size_t> > unused_id_polygons;
  bool success = CGAL::IO::internal::read_OBJ(in, points, id_polylines, unused_id_polygons);
  if(!success)
    return false;

  double squared_distance = (std::numeric_limits<double>::max)();
  for(const std::vector<std::size_t>& id_pl : id_polylines) {
    for(std::size_t pid = 1; pid < id_pl.size(); ++pid) {
      assert(points[id_pl[pid-1]] != points[id_pl[pid]]);
      if(CGAL::squared_distance(points[id_pl[pid-1]], points[id_pl[pid]]) < squared_distance)
        squared_distance = CGAL::squared_distance(points[id_pl[pid-1]], points[id_pl[pid]]);
      cdt.insert_constraint(points[id_pl[pid-1]], points[id_pl[pid]]);
    }
  }

  CDT::Face_circulator fc = cdt.incident_faces(cdt.infinite_vertex()), done = fc;
  do{
    fc->set_constraint(fc->index(cdt.infinite_vertex()), true);
    ++fc;
  } while (fc != done);

  std::cout << "Before conforming Gabriel: "
            << cdt.number_of_vertices() << " vertices.\n"
            << "  Smallest squared distance between constraint endpoints: "
            << squared_distance << "\n";

  CGAL::make_conforming_Gabriel_2(cdt);
  std::cout << "done" << std::endl;
  return true;
}
