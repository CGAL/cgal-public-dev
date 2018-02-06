#include <CGAL/Polygon_mesh_processing/internal/Hole_filling/island_triangulate_hole_polyline.h>
#include <CGAL/Polygon_mesh_processing/internal/Hole_filling/Triangulate_hole_polyline.h>


namespace CGAL {
namespace Polygon_mesh_processing{


template<typename PointRange>
using Domain = CGAL::internal::Domain<PointRange>;


template <typename PointRange, typename PolygonMesh>
std::size_t triangulate_hole_islands(PointRange boundary, PointRange hole, PolygonMesh& mesh)
{
  std::cout << "triangulate_hole_islands" << std::endl;

  // create domain from the boundary
  boundary.pop_back(); // remove the last(=first) stupid point

  // indices
  std::vector<int> b_indices;
  for(std::size_t i = 0; i < boundary.size(); ++i)
    b_indices.push_back(i);

  Domain<PointRange> domain(b_indices);


  // add hole if there is one
  std::vector<int> h_ids;
  if(!hole.empty())
  {
    hole.pop_back();
    std::size_t n_b =  b_indices.size();
    for(std::size_t i = n_b; i < n_b + hole.size(); ++i)
      h_ids.push_back(i);
    domain.add_hole(h_ids);
  }

  // assign access edge on the boundary
  const int i =  static_cast<int>(b_indices.size())-1; // todo: switch these
  const int k = 0;
  std::size_t count = 0;

  // weight calculator
  typedef CGAL::internal::Weight_min_max_dihedral_and_area      Weight;
  typedef CGAL::internal::Weight_calculator<Weight,
                CGAL::internal::Is_not_degenerate_triangle>  WC;

  // lookup tables
  typedef CGAL::internal::Lookup_table_map<Weight> WeightTable;
  typedef CGAL::internal::Lookup_table_map<int> LambdaTable;
  int n = static_cast<int>(b_indices.size() + h_ids.size()); // todo: function to return their sum
  WeightTable W(n, Weight::DEFAULT());
  LambdaTable lambda(n, -1);

  // put together points list
  PointRange Points;
  Points.reserve(boundary.size() + hole.size());
  Points.insert(Points.end(), boundary.begin(), boundary.end());
  if(!hole.empty())
  {
    Points.insert(Points.end(), hole.begin(), hole.end());
  }

  // output triangulation
  std::vector<std::vector<std::size_t>> triplets;

  CGAL::internal::Triangulate<PointRange, WC, WeightTable, LambdaTable>
      triangulation(domain, Points, W, lambda, WC());
  triangulation.do_triangulation(i, k, count);
  triangulation.collect_triangles(triplets, k, i);

  triangulation.visualize(Points, triplets, mesh);

  return count;
}



} //namespace Polygon_mesh_processing
} //namespace CGAL
