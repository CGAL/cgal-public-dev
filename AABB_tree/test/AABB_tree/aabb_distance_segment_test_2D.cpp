// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
//
// Author(s)     : Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#include <iostream>
#include <fstream>
#include <CGAL/Timer.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Segment_2.h>

#include "AABB_test_util_2D.h"


template<typename K, typename Tree, typename Primitive>
void test_impl(Tree& tree, const double duration)
{
  tree.accelerate_distance_queries();

  test_distance_speed<Tree,K>(tree,duration);
 // test_all_distance_query_types<Tree,K>(tree);
}

int main(void)
{
  std::cout << "AABB distance to triangle tests" << std::endl;
  const double duration = 0.1;
  test_kernels<SEGMENT>("./data/2D/france.seg",duration);
  test_kernels<SEGMENT>("./data/2D/corsica.seg",duration);
  test_kernels<SEGMENT>("./data/2D/estonia.seg",duration);
  test_kernels<SEGMENT>("./data/2D/island.seg",duration);
  return 0;
}
