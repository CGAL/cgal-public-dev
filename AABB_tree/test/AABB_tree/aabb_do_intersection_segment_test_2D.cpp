/*
 * aabb_do_intersection_segment_2D.h
 *
 *  Created on: 22 juil. 2013
 *      Author: Mehdious
 */

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

  test_do_intersection_speed<Tree,K>(tree,duration,RAY_QUERY);

  test_do_intersection_speed<Tree,K>(tree,duration,SEGMENT_QUERY);

  test_do_intersection_speed<Tree,K>(tree,duration,LINE_QUERY);

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