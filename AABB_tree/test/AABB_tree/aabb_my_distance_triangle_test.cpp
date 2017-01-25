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
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>

#include "AABB_test_util.h"

class InputParser{
    public:
        InputParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }
        /// @author iain
        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            return empty_string;
        }
        /// @author iain
        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
        std::string empty_string;
};

template <class K>
bool test(const char* filename, const double duration,
          bool verbose, bool kd_tree)
{
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;
  Polyhedron polyhedron;
  std::ifstream ifs(filename);
  ifs >> polyhedron;
  if(!ifs) return false;
  Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
  test_impl(tree, polyhedron, duration, verbose, kd_tree);
  return true;
}

template<class K, class Tree, class Polyhedron>
void test_impl(Tree& tree, Polyhedron& p, const double duration,
               bool verbose, bool kd_tree)
{
#ifndef NO_KDTREE
  tree.accelerate_distance_queries(p.points_begin(),p.points_end());
#endif

  test_distance_speed<Tree,K>(tree,duration);
  test_all_distance_query_types<Tree,K>(tree);
}

template<Primitive_type Primitive>
void test_kernels(const char *filename,
                  const double duration)
{
    std::cout << std::endl;
    std::cout << "Polyhedron " << filename << std::endl;
    std::cout << "============================" << std::endl;

    std::cout << std::endl;
    std::cout << "Simple cartesian float kernel" << std::endl;
    test<CGAL::Simple_cartesian<float>,Primitive>(filename,duration);

    std::cout << std::endl;
    std::cout << "Cartesian float kernel" << std::endl;
    test<CGAL::Cartesian<float>,Primitive>(filename,duration);

    std::cout << std::endl;
    std::cout << "Simple cartesian double kernel" << std::endl;
    test<CGAL::Simple_cartesian<double>,Primitive>(filename,duration);

    std::cout << std::endl;
    std::cout << "Cartesian double kernel" << std::endl;
    test<CGAL::Cartesian<double>,Primitive>(filename,duration);

    std::cout << std::endl;
    std::cout << "Epic kernel" << std::endl;
    test<CGAL::Exact_predicates_inexact_constructions_kernel,Primitive>(filename,duration);
}

int main(int argc, char** argv)
{
  Input_parser input(argc, argv);

  const char* filename = argc > 0 ? argv[1] : "./data/cube.off";
  std::cout << "AABB distance to triangle tests" << std::endl;
  double duration = 5;
  std::string duration_str = input.getCmdOption("-d");
  if(!duration_str.empty()) duration = std::sttod(duration_str);
  test_kernels<TRIANGLE>("./data/cube-fei-int.off",duration);
  return 0;
}
