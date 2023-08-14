// Copyright (c) 2023
// INRIA Sophia-Antipolis (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef COLLISIONS_3_H
#define COLLISIONS_3_H

#include <Collision_candidate_3.h>
#include <Collision_mesh_3.h>
#include <Collision_scene_3.h>
#include <vector>
#include <algorithm>
#include <iostream>

namespace CGAL{

  template <class K>
  bool do_collide(
      std::vector< Collision_mesh<K> >& meshes
  ){
    Collision_scene<K> scene = Collision_scene<K>(meshes);
    return has_collision(scene); //
  }

  template <class K>
  bool has_collision(
      const Collision_scene<K>& scene
  ){
    using Scene_face_index = Collision_scene<K>::Scene_face_index;
    using Candidate = Collision_candidate<Scene_face_index>;

    std::vector<Candidate> candidates = get_collision_candidates(scene);

    for( const auto& candidate : candidates )
    {
      if( candidate_has_collision(scene, candidate) ){
        return true;
      }
    }
    return false; //
  }

  template <class K>
  bool candidate_has_collision(
      const Collision_candidate<typename Collision_scene<K>::Trajectory>& candidate 
  ){

    using Trajectory = Collision_scene<K>::Trajectory;

    Trajectory* trajectory_0{candidate.first};
    Trajectory* trajectory_1{candidate.second};

    size_t collision_parity{0};
    for( const auto& t : scene.trajectories_ )
    {
      collision_parity += scene.tree().do_intersect(t.bounding_iso_cuboid);
    }

    // If the number of ray-intersections is odd, 
    // then a collision has occurred.
    return static_cast<bool>(collision_parity%2);                                      
  }

  

  template <class K>
  bool candidate_has_edge_edge_collision(
      const Collision_candidate<typename Collision_scene<K>::Trajectory>& candidate 
  ){

    using Trajectory = Collision_scene<K>::Trajectory;

    size_t collision_parity{0};
    for( const auto& t : scene.trajectories_ )
    {
      collision_parity += scene.tree().do_intersect(t.bounding_iso_cuboid);
    }

    // If the number of ray-intersections is odd, 
    // then a collision has occurred.
    return static_cast<bool>(collision_parity%2);                                      
  }

  template <class K>
  auto get_collision_candidates(
    const Collision_scene<K>& scene
  ) -> std::vector< Collision_candidate<typename Collision_scene<K>::Trajectory> >
  {      
      typedef Collision_candidate< 
        typename Collision_scene<K>::Trajectory 
      > Candidate;

      const auto candidate_primitive_pairs = scene.tree().all_self_intersections();

      std::vector<Candidate> candidates;
      candidates.reserve(candidate_primitive_pairs.size());

      for( const auto& primitive_pair : candidate_primitive_pairs ) {
        candidates.push_back(
          Candidate(primitive_pair)
        );
      }
      
      return candidates;
  }


} // end CGAL
#endif