#pragma once

#include <string>
#include <vector>
#include <array>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector_types.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "CgalDefinitions.h"

namespace PolyscopeBridge {

// Global cumulative offset tracker for Mesh B
inline glm::vec3 g_currentTranslationB(0.0f, 0.0f, 0.0f);

// Helper: Convert CUDA float3 to GLM vec3
inline glm::vec3 toGlm(float3 v) {
  return glm::vec3(v.x, v.y, v.z);
}

// Helper: Build 4x4 matrix with rotation around pivot center + translation
inline glm::mat4 buildTransformMatrix(glm::vec3 rotDeg, glm::vec3 trans, glm::vec3 center = glm::vec3(0.0f)) {
  glm::mat4 mat(1.0f);
  
  // 1. Apply translation and shift back from center
  mat = glm::translate(mat, trans + center);
  
  // 2. Apply rotations (R = Rz * Ry * Rx)
  mat = glm::rotate(mat, glm::radians(rotDeg.z), glm::vec3(0.0f, 0.0f, 1.0f));
  mat = glm::rotate(mat, glm::radians(rotDeg.y), glm::vec3(0.0f, 1.0f, 0.0f));
  mat = glm::rotate(mat, glm::radians(rotDeg.x), glm::vec3(1.0f, 0.0f, 0.0f));
  
  // 3. Shift to origin before rotation
  mat = glm::translate(mat, -center);
  
  return mat;
}

// 1. Multi-threaded CGAL Surface_mesh conversion via TBB with custom base color
inline void registerCgalMesh(const std::string& name, const Mesh& cgalMesh, glm::vec3 color = glm::vec3(0.8f, 0.8f, 0.8f)) {
  size_t nVerts = num_vertices(cgalMesh);
  size_t nFaces = num_faces(cgalMesh);

  std::vector<glm::vec3> points(nVerts);
  std::vector<std::array<size_t, 3>> faces(nFaces);

  auto coords = cgalMesh.points();

  // Parallel Vertex Extraction
  tbb::parallel_for(tbb::blocked_range<size_t>(0, nVerts), [&](const tbb::blocked_range<size_t>& r) {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      typename Mesh::Vertex_index v(static_cast<typename Mesh::size_type>(i));
      const auto& p = coords[v];
      points[i] = glm::vec3(static_cast<float>(p.x()), 
                            static_cast<float>(p.y()), 
                            static_cast<float>(p.z()));
    }
  });

  // Parallel Face Topology Extraction
  tbb::parallel_for(tbb::blocked_range<size_t>(0, nFaces), [&](const tbb::blocked_range<size_t>& r) {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      typename Mesh::Face_index f(static_cast<typename Mesh::size_type>(i));
      auto h0 = cgalMesh.halfedge(f);
      auto h1 = cgalMesh.next(h0);
      auto h2 = cgalMesh.next(h1);

      faces[i] = {
        static_cast<size_t>(cgalMesh.target(h0)),
        static_cast<size_t>(cgalMesh.target(h1)),
        static_cast<size_t>(cgalMesh.target(h2))
      };
    }
  });

  // Register with Polyscope
  auto* psMesh = polyscope::registerSurfaceMesh(name, points, faces);
  
  if (psMesh) {
    // Set base mesh color (RGB values from 0.0 to 1.0)
    psMesh->setSurfaceColor(color);

    // Turn off wireframes for >1M faces to maintain rendering performance
    if (nFaces > 1000000) {
      psMesh->setEdgeWidth(0.0);
    }
  }
}

// 2. Initialize Polyscope engine
inline void init() {
  polyscope::options::autocenterStructures = false; // Preserve explicit coordinate space
  polyscope::init();
}

// 3. Reset visual state & auto-fit camera with distinct mesh colors
inline void reset(const Mesh& meshA, const Mesh& meshB) {
  g_currentTranslationB = glm::vec3(0.0f, 0.0f, 0.0f);

  // Mesh A = Soft Blue, Mesh B = Soft Green
  registerCgalMesh("Mesh A", meshA, glm::vec3(0.2f, 0.5f, 0.9f));
  registerCgalMesh("Mesh B", meshB, glm::vec3(0.2f, 0.8f, 0.3f));

  // Reset transforms on registered Polyscope surface structures
  if (auto* psMeshA = polyscope::getSurfaceMesh("Mesh A")) {
    psMeshA->setTransform(glm::mat4(1.0f));
  }
  if (auto* psMeshB = polyscope::getSurfaceMesh("Mesh B")) {
    psMeshB->setTransform(glm::mat4(1.0f));
  }

  // Frames both meshes cleanly in camera view
  polyscope::view::resetCameraToHomeView();
}

// 4. Apply relative translation to Mesh B in Polyscope (PRESERVED FOR BACKWARD COMPATIBILITY)
inline void translateMeshB(float x, float y, float z) {
  g_currentTranslationB = glm::vec3(x, y, z);

  auto* meshB = polyscope::getSurfaceMesh("Mesh B");
  if (meshB) {
    glm::mat4 xf = glm::translate(glm::mat4(1.0f), g_currentTranslationB);
    meshB->setTransform(xf);
  }
}

// 4b. NEW: Apply rotation + translation to Mesh A (glm::vec3 parameters)
inline void transformMeshA(glm::vec3 rotDeg, glm::vec3 trans, glm::vec3 center = glm::vec3(0.0f)) {
  if (auto* meshA = polyscope::getSurfaceMesh("Mesh A")) {
    meshA->setTransform(buildTransformMatrix(rotDeg, trans, center));
  }
}

// 4c. NEW: Apply rotation + translation to Mesh B (glm::vec3 parameters)
inline void transformMeshB(glm::vec3 rotDeg, glm::vec3 trans, glm::vec3 center = glm::vec3(0.0f)) {
  g_currentTranslationB = trans;
  if (auto* meshB = polyscope::getSurfaceMesh("Mesh B")) {
    meshB->setTransform(buildTransformMatrix(rotDeg, trans, center));
  }
}

// 4d. NEW: CUDA float3 Overloads for transform
inline void transformMeshA(float3 rotDeg, float3 trans, float3 center = make_float3(0,0,0)) {
  transformMeshA(toGlm(rotDeg), toGlm(trans), toGlm(center));
}

inline void transformMeshB(float3 rotDeg, float3 trans, float3 center = make_float3(0,0,0)) {
  transformMeshB(toGlm(rotDeg), toGlm(trans), toGlm(center));
}

// 4e. NEW: Set transforms for both meshes simultaneously (matches KernelBVHController workflow)
inline void transformBoth(float3 rotDegA, float3 transA, float3 centerA,
                          float3 rotDegB, float3 transB, float3 centerB) {
  transformMeshA(rotDegA, transA, centerA);
  transformMeshB(rotDegB, transB, centerB);
}

// 5. Highlight intersecting faces using explicit RGB values
inline void highlightIntersections(const std::vector<int2>& exactPairs, size_t numFacesA, size_t numFacesB) {
  auto* meshA = polyscope::getSurfaceMesh("Mesh A");
  auto* meshB = polyscope::getSurfaceMesh("Mesh B");

  // Define distinct base colors and the yellow intersection color
  const glm::vec3 baseColorA(0.2f, 0.5f, 0.9f); // Blue
  const glm::vec3 baseColorB(0.2f, 0.8f, 0.3f); // Green
  const glm::vec3 highlightColor(1.0f, 0.9f, 0.0f); // Bright Yellow

  std::vector<glm::vec3> colorsA(numFacesA, baseColorA);
  std::vector<glm::vec3> colorsB(numFacesB, baseColorB);

  // Parallel face coloring
  tbb::parallel_for(tbb::blocked_range<size_t>(0, exactPairs.size()), [&](const tbb::blocked_range<size_t>& r) {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      const auto& pair = exactPairs[i];
      if (pair.x >= 0 && static_cast<size_t>(pair.x) < numFacesA) colorsA[pair.x] = highlightColor;
      if (pair.y >= 0 && static_cast<size_t>(pair.y) < numFacesB) colorsB[pair.y] = highlightColor;
    }
  });

  if (meshA) {
    meshA->removeQuantity("Intersections"); 
    meshA->addFaceColorQuantity("Intersections", colorsA)->setEnabled(true);
  }
  if (meshB) {
    meshB->removeQuantity("Intersections");
    meshB->addFaceColorQuantity("Intersections", colorsB)->setEnabled(true);
  }
}

// 6. Draw one frame without blocking the terminal prompt
inline void drawFrame() {
  polyscope::frameTick();
}

} // namespace PolyscopeBridge