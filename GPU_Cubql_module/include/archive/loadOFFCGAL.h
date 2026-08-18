#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <iostream>

#include "cuBQL/bvh.h" 
#include "cuBQL/math/vec.h"
#include <cuBQL/queries/triangleData/Triangle.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

namespace cuBQL {
namespace samples {

  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Point  = Kernel::Point_3;
  using Mesh   = CGAL::Surface_mesh<Point>;

  inline std::vector<Triangle>
  loadOFFCGAL(const std::string& offFile)
  {
    // Fixes the standard I/O synchronization bottleneck
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    Mesh mesh;
    std::ifstream in(offFile);
    if (!in.is_open()) {
      throw std::runtime_error("Could not open OFF file: " + offFile);
    }
    if (!(in >> mesh)) {
      throw std::runtime_error("Failed to read OFF file with CGAL: " + offFile);
    }

    // Safely cache point positions up front to eliminate high-frequency lookup costs
    std::vector<vec3f> cached_vertices;
    cached_vertices.reserve(mesh.number_of_vertices());
    for (auto v : mesh.vertices()) {
      const auto& p = mesh.point(v);
      cached_vertices.push_back(vec3f{(float)p.x(), (float)p.y(), (float)p.z()});
    }

    std::vector<Triangle> triangles;
    triangles.reserve(mesh.number_of_faces());

    // Loop through faces using safe, standard member functions that optimize down perfectly
    for (auto f : mesh.faces()) {
      auto h0 = mesh.halfedge(f);
      auto h1 = mesh.next(h0);
      auto h2 = mesh.next(h1);

      size_t i0 = (size_t)mesh.target(h0);
      size_t i1 = (size_t)mesh.target(h1);
      size_t i2 = (size_t)mesh.target(h2);

      triangles.push_back({cached_vertices[i0], cached_vertices[i1], cached_vertices[i2]});
    }

    return triangles;
  }

} // namespace samples
} // namespace cuBQL