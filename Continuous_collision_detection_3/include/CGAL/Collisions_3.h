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
#include <Bilinear_patch_3.h>
#include <Ray_3_Bilinear_patch_3_do_intersect.h>

// New collision objects
#include <Collision_mesh_3.h>
#include <Collision_scene_3.h>
#include <Trajectories.h>

// New collision routines
#include <Segment_3_Segment_3_do_collide.h>
#include <Point_3_Triangle_3_do_collide.h>
#include <Triangle_3_Triangle_3_do_collide.h>
#include <Collision_scene_3_has_collision.h>
#include <Collision_mesh_3_Collision_mesh_3_do_collide.h>

#endif