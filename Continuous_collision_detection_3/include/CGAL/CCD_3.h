// Copyright (c) 2023 GeometryFactory (France).
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef COLLISIONS_3_H
#define COLLISIONS_3_H

// New geometric objects and corresponding routines
#include <CGAL/Bilinear_patch_3.h>
#include <CGAL/Ray_3_Bilinear_patch_3_do_intersect.h>

// New collision objects
#include <CGAL/Collision_mesh_3.h>
#include <CGAL/Collision_scene_3.h>
#include <CGAL/Trajectories.h>

// New collision routines
#include <CGAL/Segment_3_Segment_3_do_collide.h>
#include <CGAL/Point_3_Triangle_3_do_collide.h>
#include <CGAL/Triangle_3_Triangle_3_do_collide.h>
#include <CGAL/Collision_scene_3_has_collision.h>
#include <CGAL/Collision_mesh_3_Collision_mesh_3_do_collide.h>

// Utility routine
#include <CGAL/draw_collision_scene.h>

// Internal routines and objects
#include <CGAL/Continuous_collision_detection_3/internal/Update_functors.h>
#include <CGAL/Continuous_collision_detection_3/internal/Segment_3_Segment_3_collision_test_boundary.h>
#include <CGAL/Continuous_collision_detection_3/internal/Point_3_Triangle_3_collision_test_boundary.h>
#include <CGAL/Continuous_collision_detection_3/internal/AABB_triangle_trajectory_primitive.h>
#include <CGAL/Continuous_collision_detection_3/internal/Cull_test_boundary.h>


#endif