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

#ifndef COLLISION_MESH_3_COLLISION_MESH_3_DO_COLLIDE_H
#define COLLISION_MESH_3_COLLISION_MESH_3_DO_COLLIDE_H

#include <CGAL/Collision_mesh_3.h>
#include <CGAL/Collision_scene_3.h>
#include <CGAL/Collision_scene_3_has_collision.h>

namespace CGAL{



  template <class K>
  bool do_collide(
      std::vector< Collision_mesh<K> >& meshes
  ){
    Collision_scene<K> scene = Collision_scene<K>(meshes);
    return has_collision(scene); //
  }



} // end CGAL
#endif