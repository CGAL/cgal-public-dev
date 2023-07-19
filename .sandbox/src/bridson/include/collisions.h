#ifndef COLLISIONS_H
#define COLLISIONS_H

#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <AABB_triangle_pair_primitive.h> 

template <class SimulationMesh, class K>
bool do_collide(
    SimulationMesh& mesh_A,
    SimulationMesh& mesh_B
){

  // Build AABB tree whose boxes bound the swept primitives of each mesh

  typedef typename K::Point_3 Point;
  typedef typename K::Plane_3 Plane;
  typedef typename K::Vector_3 Vector;
  typedef typename K::Segment_3 Segment;

  std::vector<triangle_trajectory<K>> triangle_trajectories_A;
  std::vector<triangle_trajectory<K>> triangle_trajectories_A;
  std::vector<Point> positions(3);
  std::vector<Point> next_positions(3);

  for( auto & current_mesh : meshes )
  {

  }

  typedef typename CGAL::AABB_traits<K, AABB_tree_triangle_trajectory_primitive<K>> AABB_traits;
  typedef typename CGAL::AABB_tree<AABB_traits> Tree;
  typedef Tree::Primitive_id Primitive_id;

  Tree tree(triangle_trajectories.begin(),triangle_trajectories.end());

  bool intersection = tree.do_intersect(triangle_trajectories[0].bounding_iso_cuboid);

  return intersection; //
}

#endif