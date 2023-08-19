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

#ifndef COLLISION_SCENE_3_HAS_COLLISION_H
#define COLLISION_SCENE_3_HAS_COLLISION_H

#include <Collision_candidate_3.h>
#include <Collision_scene_3.h>
#include <AABB_triangle_trajectory_primitive.h>
#include <Triangle_3_Triangle_3_do_collide.h>
#include <vector>
#include <iostream>

namespace CGAL{

  template <class K>
  std::vector<Collision_candidate<typename Collision_scene<K>::Trajectory>> 
  get_collisions(
      const Collision_scene<K>& scene
  ){
    using Candidate = Collision_candidate<Collision_scene<K>::Trajectory>;

    std::vector<Candidate> candidates = get_collision_candidates(scene);
    std::vector<Candidate> collisions;

    std::copy_if(
        candidates.begin(), 
        candidates.end(), 
        std::back_inserter(collisions), 
        [](const auto& candidate){ return candidate_has_collision(candidate);} 
    );

    return collisions; //
  }

  template <class K>
  bool has_collision(
      const Collision_scene<K>& scene
  ){
    using Candidate = Collision_candidate<Collision_scene<K>::Trajectory>;

    std::vector<Candidate> candidates = get_collision_candidates(scene);

    for( const auto& candidate : candidates )
    {
      if( candidate_has_collision(candidate) ){
        return true;
      }
    }
    return false; //
  }

  template <class CollisionCandidate>
  bool candidate_has_collision(
      const CollisionCandidate& candidate 
  ){

    using K     = typename CollisionCandidate::K;
    using Index = typename CollisionCandidate::Index;

    Triangle_3_trajectory<K> trajectory_0 = to_Triangle_3_trajectory<K, Index>(*candidate.first);
    Triangle_3_trajectory<K> trajectory_1 = to_Triangle_3_trajectory<K, Index>(*candidate.second);

    return do_collide(trajectory_0, trajectory_1);                                      
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

      // TODO: modify collision candidates to limit the testing
      //       of adjacent faces so that they ignore common vertices/edges.
      //       Once this is done, can consider self-intersections again.
      for( const auto& primitive_pair : candidate_primitive_pairs ) {
        if( primitive_pair.first->index.mesh_index() !=  primitive_pair.second->index.mesh_index()) {
            candidates.push_back(
                Candidate(primitive_pair)
            );
        }
      }

      return candidates;
  }


} // end CGAL
#endif