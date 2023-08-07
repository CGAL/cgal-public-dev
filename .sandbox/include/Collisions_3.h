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
    return do_collide(scene); //
  }

  template <class K>
  bool do_collide(
      const Collision_scene<K>& scene
  ){

    bool do_collide_{false};
    for( const auto& t : scene.trajectories_ )
    {
      do_collide_ = do_collide_ || scene.tree().do_intersect(t.bounding_iso_cuboid);
    }
    return do_collide_; //
  }

  template <class K>
  std::vector<
    Collision_candidate<
      typename Collision_scene<K>::Scene_face_index
    > 
  > get_collision_candidates(
    const Collision_scene<K>& scene
  ){      
      typedef Collision_candidate< 
        typename Collision_scene<K>::Scene_face_index 
      > Candidate;

      std::vector<Candidate> candidates;
      const auto candidate_primitive_pairs = scene.tree().all_self_intersections();
      candidates.reserve(candidate_primitive_pairs.size());

      std::transform(
        candidate_primitive_pairs.begin(), 
        candidate_primitive_pairs.end(), 
        std::back_inserter(candidates),
        [](const auto & primitive_pair){
          return Candidate(
            primitive_pair.first->index,
            primitive_pair.second->index
          );
        }
      );
      return candidates;
  }

  template <class K>
  std::vector<
    Collision_candidate<
      typename Collision_scene<K>::Scene_face_index
    > 
  > get_collision_candidates(
    std::vector<Collision_mesh<K>>& meshes
  ){      
      Collision_scene<K> scene = Collision_scene<K>(meshes);
      return get_collision_candidates(scene);
  }

} // end CGAL
#endif