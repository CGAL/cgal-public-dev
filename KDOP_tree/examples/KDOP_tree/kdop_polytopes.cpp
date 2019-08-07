/*
 * An example to visualise the k-dops with polytopes.
 */

#include <iostream>
#include <fstream>
#include <list>

#include <boost/lexical_cast.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/KDOP_tree/KDOP_tree.h>
#include <CGAL/KDOP_tree/KDOP_traits.h>
#include <CGAL/KDOP_tree/KDOP_kdop.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/point_generators_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;
typedef K::Segment_3 Segment;

typedef K::Plane_3 Plane;

typedef CGAL::Surface_mesh<Point> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;

// prescribed number of directions for the k-dop
const unsigned int NUM_DIRECTIONS = 6;

typedef CGAL::KDOP_tree::KDOP_traits<NUM_DIRECTIONS, K, Primitive> Traits;
typedef CGAL::KDOP_tree::KDOP_tree<Traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

template<typename NamedParameters>
void write_ply(std::ostream& os, const Mesh& sm, const NamedParameters& np);

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/tetrahedron.off";

  std::ifstream input(filename);

  Mesh mesh;
  input >> mesh;

  // read the mesh into the k-dop tree
  Tree tree(faces(mesh).first, faces(mesh).second, mesh);

  // user-defined directions for the k-dop, must be compatible to NUM_DIRECTIONS.
  std::vector<Vector> directions;
  for (int i = 0; i < 3; ++i) {
    std::vector<double> v(3, 0.);
    v[i] = 1.;

    Vector direction(v[0], v[1], v[2]);
    directions.push_back(direction);
  }
  if (NUM_DIRECTIONS == 14 || NUM_DIRECTIONS == 26) {
    directions.push_back(Vector(1., 2., 1.));
    directions.push_back(Vector(-1., 1., 1.));
    directions.push_back(Vector(-1., -1., 1.));
    directions.push_back(Vector(1., -1., 1.));
  }
  if (NUM_DIRECTIONS == 18 || NUM_DIRECTIONS == 26) {
    directions.push_back(Vector(1., 1., 0.));
    directions.push_back(Vector(1., 0., 1.));
    directions.push_back(Vector(0., 1., 1.));
    directions.push_back(Vector(1., -1., 0.));
    directions.push_back(Vector(1., 0., -1.));
    directions.push_back(Vector(0., 1., -1.));
  }

  // the remaining directions are opposite to the first k/2 directions
  for (int i = 0; i < NUM_DIRECTIONS/2; ++i) {
    Vector direction = directions[i];

    Vector direction1(-direction[0], -direction[1], -direction[2]);
    directions.push_back(direction1);
  }

  // give the user-defined directions to the k-dop tree
  tree.set_kdop_directions(directions);

  // build the tree, including splitting primitives and computing k-dops with
  // the directions given
  tree.build();

  typedef typename Tree::Kdop Kdop;
  std::vector< typename Kdop::Array_height > heights;

  tree.kdop_heights(heights); // get all k-dops of the tree

  std::cout << "number of polytopes: " << heights.size() << std::endl;

  // retrieve directions for k-dops, necessary if the directions are not
  // explicitly given
  typedef CGAL::KDOP_tree::Construct_kdop<K, NUM_DIRECTIONS> Construct_kdop;
  Construct_kdop construct_kdop(directions);

  std::vector< Vector > kdop_directions = construct_kdop.kdop_directions();

  // create convex hulls by intersecting half planes
  for (int i = 0; i < heights.size(); ++i) {
    std::list<Plane> planes;
    typename Kdop::Array_height height = heights[i];

    for (int j = 0; j < NUM_DIRECTIONS; ++j) {
      Vector v = kdop_directions[j];
      const double v_length = std::sqrt(v.squared_length());
      v = v / v_length;

      Plane plane(v.x(), v.y(), v.z(), -height[j]);
      planes.push_back(plane);
    }

    Mesh chull;
    CGAL::halfspace_intersection_3(planes.begin(),
                                   planes.end(),
                                   chull);

    // output polytopes representing k-dops
    std::string meshFile = "out/polytope" + boost::lexical_cast< std::string >(i) + ".ply";
    std::ofstream meshf(meshFile.c_str());

    write_ply(meshf, chull, CGAL::parameters::all_default());
  }

  return 0;
}

template<typename NamedParameters>
void write_ply(std::ostream& os, const Mesh& sm, const NamedParameters& np)
{
  typedef typename Mesh::Vertex_index Vertex_index;
  typedef typename Mesh::Face_index Face_index;

  os << "ply\n" << "format ascii 1.0\n";
  os << "element vertex " << sm.number_of_vertices() << "\n";
  os << "property float x\n" << "property float y\n" << "property float z\n";
  os << "element face " << sm.number_of_faces() << "\n";
  os << "property list uchar int vertex_indices\n";
  os << "end_header\n";

  std::vector<int> reindex;
  typename CGAL::Polygon_mesh_processing::GetVertexPointMap<Mesh, NamedParameters>::const_type
          vpm = choose_param(get_param(np, CGAL::internal_np::vertex_point),
                             CGAL::get_const_property_map(CGAL::vertex_point, sm));
  reindex.resize(sm.num_vertices());
  int n = 0;
  for(Vertex_index v : sm.vertices()){
    os << get(vpm, v);
    os << '\n';
    reindex[v]=n++;
  }

  for(Face_index f : sm.faces()){
    os << sm.degree(f);
    for(Vertex_index v : CGAL::vertices_around_face(sm.halfedge(f),sm)){
      os << " " << reindex[v];
    }
    os << '\n';
  }

}



