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

#include <Collision_mesh_3.h>
#include <Collision_scene_3.h>

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
      do_collide_ = do_collide_ || scene.trajectory_tree.do_intersect(t.bounding_iso_cuboid);
    }
    return do_collide_; //
  }

  template <class K, class OutputIterator>
  OutputIterator get_collision_candidates(
    const Collision_scene<K>& scene
  ){

      OutputIterator primitives;
      for( const auto& t : scene.trajectories_ )
      {
        scene.trajectory_tree.all_intersected_primitives(t.bounding_iso_cuboid, std::back_inserter(primitives));
      }
      
      return primitives;

  }

} // end CGAL
#endif