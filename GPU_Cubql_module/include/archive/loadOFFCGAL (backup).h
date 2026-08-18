#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>

// Make sure we include the explicit primitive definitions from cuBQL
#include "cuBQL/bvh.h" 
#include "cuBQL/math/vec.h"
//#include "loadOBJ.h"
#include <cuBQL/queries/triangleData/Triangle.h>
//#include "cuBQL/math/triangle.h" // or wherever cuBQL defines its Triangle struct

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

namespace cuBQL {
namespace samples {

  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Point  = Kernel::Point_3;
  using Mesh   = CGAL::Surface_mesh<Point>;
  namespace PMP = CGAL::Polygon_mesh_processing;

  // Stripping the cuBQL:: prefix inside the signature since we are already 
  // inside namespace cuBQL. This prevents the compiler from getting confused 
  // by CGAL's deep lookups.
  inline std::vector<Triangle>
  loadOFFCGAL(const std::string& offFile)
  {
    Mesh mesh;
    std::ifstream in(offFile);
    if (!in.is_open()) {
      throw std::runtime_error("Could not open OFF file: " + offFile);
    }
    if (!(in >> mesh)) {
      throw std::runtime_error("Failed to read OFF file with CGAL: " + offFile);
    }

    PMP::triangulate_faces(mesh);

    std::vector<Triangle> triangles;
    triangles.reserve(mesh.number_of_faces());

    for (auto f : mesh.faces()) {
      auto h = mesh.halfedge(f);
      auto v0 = mesh.target(h);
      h = mesh.next(h);
      auto v1 = mesh.target(h);
      h = mesh.next(h);
      auto v2 = mesh.target(h);

      const auto& p0 = mesh.point(v0);
      const auto& p1 = mesh.point(v1);
      const auto& p2 = mesh.point(v2);

      // Using direct types since we are in the namespace scope
      triangles.emplace_back(
        vec3f{(float)p0.x(), (float)p0.y(), (float)p0.z()},
        vec3f{(float)p1.x(), (float)p1.y(), (float)p1.z()},
        vec3f{(float)p2.x(), (float)p2.y(), (float)p2.z()}
      );
    }
    return triangles;
  }

} // namespace samples
} // namespace cuBQL