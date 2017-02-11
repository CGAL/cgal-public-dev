#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
// #include <CGAL/Extended_cartesian.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/dual.h>
#include <CGAL/IO/print_OFF.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <utility>
#include <map>
#include <iterator>
#include <algorithm>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K>                   Delaunay;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;
typedef CGAL::Polyhedron_incremental_builder_3<Polyhedron::HalfedgeDS> Polyhedron_incremental_builder;
typedef CGAL::Bbox_3                                        Bbox;

typedef Delaunay::Vertex_handle            Vertex_handle;
typedef Delaunay::Finite_vertices_iterator Finite_vertices_iterator;
typedef Delaunay::Finite_facets_iterator   Finite_facets_iterator;
typedef Delaunay::Cell_handle              Cell_handle;
typedef Delaunay::Point                    Point;

namespace CGAL {

template <class HDS>
class Build_cube : public CGAL::Modifier_base<HDS> {
public:
  Build_cube(const Bbox_3& bbox): _bbox(bbox) {}
  void operator()(HDS& hds) {
    typedef typename HDS::Vertex   Vertex;
    typedef typename Vertex::Point Point;

    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
    B.begin_surface(8, 6);

    B.add_vertex( Point(_bbox.xmin(), _bbox.ymin(), _bbox.zmin()));
    B.add_vertex( Point(_bbox.xmax(), _bbox.ymin(), _bbox.zmin()));
    B.add_vertex( Point(_bbox.xmax(), _bbox.ymax(), _bbox.zmin()));
    B.add_vertex( Point(_bbox.xmin(), _bbox.ymax(), _bbox.zmin()));
    B.add_vertex( Point(_bbox.xmin(), _bbox.ymin(), _bbox.zmax()));
    B.add_vertex( Point(_bbox.xmax(), _bbox.ymin(), _bbox.zmax()));
    B.add_vertex( Point(_bbox.xmax(), _bbox.ymax(), _bbox.zmax()));
    B.add_vertex( Point(_bbox.xmin(), _bbox.ymax(), _bbox.zmax()));

    B.begin_facet();
    B.add_vertex_to_facet( 3);
    B.add_vertex_to_facet( 2);
    B.add_vertex_to_facet( 1);
    B.add_vertex_to_facet( 0);
    B.end_facet();

    B.begin_facet();
    B.add_vertex_to_facet( 4);
    B.add_vertex_to_facet( 5);
    B.add_vertex_to_facet( 6);
    B.add_vertex_to_facet( 7);
    B.end_facet();

    B.begin_facet();
    B.add_vertex_to_facet( 5);
    B.add_vertex_to_facet( 4);
    B.add_vertex_to_facet( 0);
    B.add_vertex_to_facet( 1);
    B.end_facet();

    B.begin_facet();
    B.add_vertex_to_facet( 6);
    B.add_vertex_to_facet( 5);
    B.add_vertex_to_facet( 1);
    B.add_vertex_to_facet( 2);
    B.end_facet();

    B.begin_facet();
    B.add_vertex_to_facet( 7);
    B.add_vertex_to_facet( 6);
    B.add_vertex_to_facet( 2);
    B.add_vertex_to_facet( 3);
    B.end_facet();

    B.begin_facet();
    B.add_vertex_to_facet( 4);
    B.add_vertex_to_facet( 7);
    B.add_vertex_to_facet( 3);
    B.add_vertex_to_facet( 0);
    B.end_facet();

    B.end_surface();
  }

private:
  Bbox_3 _bbox;
};

// Print a triangulation (off format) as a set of finite triangles
template <class Triang>
void
triangulation_print_OFF(std::ostream& out, const Triang& triang) {
  size_t v_nb = 0;
  size_t f_nb = 0;
  std::map<Vertex_handle, size_t> T_v_ids;
  std::stringstream triangulation_out_tris;
  for(
    Finite_facets_iterator
      fit = triang.finite_facets_begin(),
      fend = triang.finite_facets_end();
    fit != fend; ++fit, ++f_nb
  ) {
    Cell_handle c = fit->first;
    int i_op = fit->second;
    for(int i = 0; i < 4; ++i) {
      if(i == i_op)
        continue;
      Vertex_handle v = c->vertex(i);
      if(T_v_ids.find(v) == T_v_ids.end()) {
        T_v_ids.insert(std::make_pair(v, v_nb++));
      }
    }
    triangulation_out_tris << 3;
    for(int i = 0; i < 4; ++i) {
      if(i == i_op)
        continue;
      triangulation_out_tris << " " << T_v_ids[c->vertex(i)];
    }
    triangulation_out_tris << std::endl;
  }
  out << "OFF" << std::endl;
  out << v_nb << " " << f_nb << " " << 0 << std::endl << std::endl;
  std::map<size_t, Point> v_pts;
  for(
    Finite_vertices_iterator
      vit = triang.finite_vertices_begin(),
      vend = triang.finite_vertices_end();
    vit != vend; ++vit
  ) {
    v_pts.insert(std::make_pair(T_v_ids[vit], vit->point()));
  }
  for(std::map<size_t, Point>::iterator vpt_it = v_pts.begin(), vpt_end = v_pts.end();
    vpt_it != vpt_end; ++ vpt_it
  ) {
    out << vpt_it->second << std::endl;
  }
  out << triangulation_out_tris.rdbuf();
}

} // namespace CGAL

int main(int argc, char** argv)
{
  std::ofstream inner_voronoi_out("inner_voronoi.off");
  if(!inner_voronoi_out.is_open()) {
    std::cerr << "Can't write to Voronoi cell" << std::endl;
    return 1;
  }

  std::ofstream inner_voronoi_with_bbox_out("inner_voronoi_with_bbox.off");
  if(!inner_voronoi_with_bbox_out.is_open()) {
    std::cerr << "Can't write to Voronoi cell" << std::endl;
    return 1;
  }

  std::ofstream hull_voronoi_out("hull_voronoi.off");
  if(!hull_voronoi_out.is_open()) {
    std::cerr << "Can't write to Voronoi cell" << std::endl;
    return 1;
  }

  std::ofstream triangulation_out("triangulation.off");
  if(!triangulation_out.is_open()) {
    std::cerr << "Can't write to Triangulation" << std::endl;
    return 1;
  }

  std::ofstream bbox_out("bbox.off");
  if(!bbox_out.is_open()) {
    std::cerr << "Can't write to Bounding box" << std::endl;
    return 1;
  }

  Delaunay T;

  // Octahedron hull and 2 points inside
  Vertex_handle v0 = T.insert(Point(-1,-1,0));
  Vertex_handle v1 = T.insert(Point(-1,1,0));
  Vertex_handle v2 = T.insert(Point(1,-1,0));
  Vertex_handle v3 = T.insert(Point(1,1,0));
  Vertex_handle v4 = T.insert(Point(0,0,1));
  Vertex_handle v5 = T.insert(Point(0,0,-1));
  Vertex_handle v6 = T.insert(Point(0.25,0.25,0.25));
  Vertex_handle v7 = T.insert(Point(-0.25,-0.25,-0.25));

  Vertex_handle inner_v = v6;
  Vertex_handle hull_v = v3;

  Bbox bbox(-1, -1, -1, 1, 1, 1);
  Polyhedron bbox_poly;
  CGAL::Build_cube<Polyhedron::HalfedgeDS> cube(bbox);
  bbox_poly.delegate(cube);
  std::vector<Delaunay::Geom_traits::Plane_3> tmp; // TODO:remove after refactoring of dual interface
  Polyhedron inner_voronoi = CGAL::dual<Polyhedron>(T, inner_v, tmp.begin(), tmp.end());
  Polyhedron inner_voronoi_with_bbox = CGAL::dual<Polyhedron>(T, inner_v, bbox);
  Polyhedron hull_voronoi = CGAL::dual<Polyhedron>(T, hull_v, bbox);

  // Write triangulation
  CGAL::triangulation_print_OFF(triangulation_out, T);
  triangulation_out.close();

  // Write dual cells
  CGAL::print_OFF(inner_voronoi_out, inner_voronoi);
  inner_voronoi_out.close();
  CGAL::print_OFF(inner_voronoi_with_bbox_out, inner_voronoi_with_bbox);
  inner_voronoi_with_bbox_out.close();
  CGAL::print_OFF(hull_voronoi_out, hull_voronoi);
  hull_voronoi_out.close();

  // Write bounding box
  CGAL::print_OFF(bbox_out, bbox_poly);
  bbox_out.close();

  return 0;
}
