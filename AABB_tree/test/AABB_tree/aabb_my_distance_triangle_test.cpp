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
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include "AABB_test_util.h"

class Input_parser{
    public:
        Input_parser(int &argc, char **argv){
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

template <class Tree, class K>
void test_distance_speed(Tree& tree,
                         const double duration,
                         bool verbose = false)
{
    // typedef typename K::FT FT;
    // typedef typename K::Ray_3 Ray;
    typedef typename K::Point_3 Point;
    // typedef typename K::Vector_3 Vector;
    if(verbose) std::cout.flush();

    CGAL::Timer timer;
    timer.start();
    unsigned int nb = 0;
    while(timer.time() < duration)
    {
            // picks a random point in the tree bbox
            Point query = random_point_in<K>(tree.bbox());
            Point closest = tree.closest_point(query);
	    (void) closest;
            if(verbose) {
              std::cerr << ".";
            }
            nb++;
    }
    if(verbose) {
      std::cerr << std::endl;
    }
    double speed = (double)nb / timer.time();
    std::cout << speed << " distance queries/s" << std::endl;
    timer.stop();
}

template<class K, class Tree, class Polyhedron>
void test_impl(Tree& tree, Polyhedron& p, const double duration,
               bool verbose, bool kd_tree)
{
  if(kd_tree) {
    tree.accelerate_distance_queries(p.points_begin(), p.points_end());
  }
  test_distance_speed<Tree,K>(tree, duration, verbose);
}

template <class K>
bool test(const char* filename, const double duration,
          bool verbose = false, bool kd_tree = true)
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
  test_impl<K>(tree, polyhedron, duration, verbose, kd_tree);
  return true;
}

int main(int argc, char** argv)
{
  Input_parser input(argc, argv);

  const char* filename = argc > 1 ? argv[1] : "./data/cube.off";
  std::cout << "AABB distance to triangle tests" << std::endl;
  double duration = 5;
  bool verbose = input.cmdOptionExists("-v");
  bool kd_tree = input.cmdOptionExists("-k");
  std::string duration_str = input.getCmdOption("-d");
  if(!duration_str.empty()) duration = std::stod(duration_str);

  typedef CGAL::Simple_cartesian<float>                       SCF;
  typedef CGAL::Cartesian<float>                              CF;
  typedef CGAL::Simple_cartesian<double>                      SCD;
  typedef CGAL::Cartesian<double>                             CD;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;

  bool at_least_one_kernel=false;
  if(input.cmdOptionExists("SCF")) {
    std::cout << "Simple cartesian float kernel" << std::endl;
    test<CGAL::Simple_cartesian<float>>(filename, duration, verbose, kd_tree);
    at_least_one_kernel=true;
  }
  if(input.cmdOptionExists("CF")) {
    std::cout << "Cartesian float kernel" << std::endl;
    test<CGAL::Cartesian<float>>(filename, duration, verbose, kd_tree);
    at_least_one_kernel=true;
  }
  if(input.cmdOptionExists("SCD")) {
    std::cout << "Simple cartesian double kernel" << std::endl;
    test<CGAL::Simple_cartesian<double>>(filename, duration, verbose, kd_tree);
    at_least_one_kernel=true;
  }
  if(input.cmdOptionExists("CD")) {
    std::cout << "Cartesian double kernel" << std::endl;
    test<CGAL::Cartesian<double>>(filename, duration, verbose, kd_tree);
    at_least_one_kernel=true;
  }
  if(!at_least_one_kernel || input.cmdOptionExists("Epick")) {
    std::cout << "Epic kernel" << std::endl;
    test<CGAL::Epick>(filename, duration, verbose, kd_tree);
  }
  return 0;
}
