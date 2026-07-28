#pragma once

#include <string>
#include <vector>
#include <array>
#include <functional>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/matrix_decompose.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <vector_types.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "imgui.h"

#include "CgalDefinitions.h"

namespace PolyscopeBridge {

// Global cumulative offset tracker for Mesh B
inline glm::vec3 g_currentTranslationB(0.0f, 0.0f, 0.0f);

// Stored rotation centers for extraction logic
inline float3 g_centerA{0.0f, 0.0f, 0.0f};
inline float3 g_centerB{0.0f, 0.0f, 0.0f};

// Interactive gizmo visibility state
inline bool g_enableGizmoA = true;
inline bool g_enableGizmoB = true;

// Base Surface Colors for Mesh A and Mesh B
inline glm::vec3 g_baseColorA(0.2f, 0.5f, 0.9f);
inline glm::vec3 g_baseColorB(0.2f, 0.8f, 0.3f);

// Custom Intersection Colors for Mesh A and Mesh B
inline glm::vec3 g_intersectColorA(1.0f, 0.9f, 0.0f); // Default Yellow
inline glm::vec3 g_intersectColorB(1.0f, 0.9f, 0.0f); // Default Yellow

// Color presets for ImGui dropdown choices
struct ColorPreset {
  const char* name;
  glm::vec3 color;
};

inline const ColorPreset g_presetColors[] = {
    {"Yellow (Default)", glm::vec3(1.0f, 0.9f, 0.0f)},
    {"Hot Red",          glm::vec3(1.0f, 0.1f, 0.1f)},
    {"Neon Cyan",        glm::vec3(0.0f, 1.0f, 1.0f)},
    {"Magenta",          glm::vec3(1.0f, 0.0f, 1.0f)},
    {"Electric Orange",  glm::vec3(1.0f, 0.5f, 0.0f)},
    {"Lime Green",       glm::vec3(0.2f, 1.0f, 0.2f)},
    {"White Hot",        glm::vec3(1.0f, 1.0f, 1.0f)}
};

inline int g_selectedPresetA = 0;
inline int g_selectedPresetB = 0;

// Intersection highlighting active state
inline bool g_showIntersections = true;

// Cache of last computed intersection pairs for dynamic UI re-coloring
inline std::vector<int2> g_lastExactPairs;
inline size_t g_lastNumFacesA = 0;
inline size_t g_lastNumFacesB = 0;

// Callback triggered when clicking "Fire Intersection" in UI
inline std::function<void()> g_onFireCallback = nullptr;

// Helper: Convert CUDA float3 to GLM vec3
inline glm::vec3 toGlm(float3 v) {
  return glm::vec3(v.x, v.y, v.z);
}

// Helper: Compare GLM vec3 values safely
inline bool isVec3Different(const glm::vec3& a, const glm::vec3& b) {
  return a.x != b.x || a.y != b.y || a.z != b.z;
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

// Decompose Polyscope model matrix into rotDeg & trans relative to pivot center
inline void
decomposeTransformMatrix(const glm::mat4& mat, glm::vec3 center, glm::vec3& outRotDeg, glm::vec3& outTrans) {
  glm::vec4 transformedCenter = mat * glm::vec4(center, 1.0f);
  outTrans = glm::vec3(transformedCenter) - center;

  glm::mat3 rotMat(mat);
  rotMat[0] = glm::normalize(rotMat[0]);
  rotMat[1] = glm::normalize(rotMat[1]);
  rotMat[2] = glm::normalize(rotMat[2]);

  float rotX, rotY, rotZ;
  // FIXED: Extracted in ZYX order to strictly match buildTransformMatrix (R = Rz * Ry * Rx)
  glm::extractEulerAngleZYX(glm::mat4(rotMat), rotZ, rotY, rotX);

  outRotDeg = glm::vec3(glm::degrees(rotX), glm::degrees(rotY), glm::degrees(rotZ));
}

// Multi-threaded CGAL Surface_mesh conversion via TBB
inline void
registerCgalMesh(const std::string& name, const Mesh& cgalMesh, glm::vec3 color = glm::vec3(0.8f, 0.8f, 0.8f)) {
  size_t nVerts = num_vertices(cgalMesh);
  size_t nFaces = num_faces(cgalMesh);

  std::vector<glm::vec3> points(nVerts);
  std::vector<std::array<size_t, 3>> faces(nFaces);

  auto coords = cgalMesh.points();

  tbb::parallel_for(tbb::blocked_range<size_t>(0, nVerts), [&](const tbb::blocked_range<size_t>& r) {
    for(size_t i = r.begin(); i != r.end(); ++i) {
      typename Mesh::Vertex_index v(static_cast<typename Mesh::size_type>(i));
      const auto& p = coords[v];
      points[i] = glm::vec3(static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()));
    }
  });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, nFaces), [&](const tbb::blocked_range<size_t>& r) {
    for(size_t i = r.begin(); i != r.end(); ++i) {
      typename Mesh::Face_index f(static_cast<typename Mesh::size_type>(i));
      auto h0 = cgalMesh.halfedge(f);
      auto h1 = cgalMesh.next(h0);
      auto h2 = cgalMesh.next(h1);

      faces[i] = {static_cast<size_t>(cgalMesh.target(h0)), static_cast<size_t>(cgalMesh.target(h1)),
                  static_cast<size_t>(cgalMesh.target(h2))};
    }
  });

  auto* psMesh = polyscope::registerSurfaceMesh(name, points, faces);
  if(psMesh) {
    psMesh->setSurfaceColor(color);
    if(nFaces > 1000000) {
      psMesh->setEdgeWidth(0.0);
    }
  }
}

// Forward declare highlighting helper function
inline void highlightIntersections(const std::vector<int2>& exactPairs, size_t numFacesA, size_t numFacesB);

// Programmatically show/hide gizmos from CLI or code
inline void setGizmosEnabled(bool showA, bool showB) {
  g_enableGizmoA = showA;
  g_enableGizmoB = showB;
  if(auto* meshA = polyscope::getSurfaceMesh("Mesh A")) {
    meshA->setTransformGizmoEnabled(g_enableGizmoA);
  }
  if(auto* meshB = polyscope::getSurfaceMesh("Mesh B")) {
    meshB->setTransformGizmoEnabled(g_enableGizmoB);
  }
}

// Updates mesh base colors dynamically
inline void updateMeshBaseColors() {
  if(auto* meshA = polyscope::getSurfaceMesh("Mesh A")) {
    meshA->setSurfaceColor(g_baseColorA);
  }
  if(auto* meshB = polyscope::getSurfaceMesh("Mesh B")) {
    meshB->setSurfaceColor(g_baseColorB);
  }
}

// GUI Callback rendered every frame in Polyscope window
inline void polyscopeUiCallback() {
  ImGui::Begin("Pipeline Controls");

  // Guard against startup state when no meshes are loaded
  if(!polyscope::hasSurfaceMesh("Mesh A") || !polyscope::hasSurfaceMesh("Mesh B")) {
    ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.2f, 1.0f), "Status: Waiting for meshes...");
    ImGui::TextUnformatted("Please run 'load' in the CLI terminal.");
    ImGui::End();
    return;
  }

  auto* meshA = polyscope::getSurfaceMesh("Mesh A");
  auto* meshB = polyscope::getSurfaceMesh("Mesh B");

  // FIXED: Sync global colors if updated via native Polyscope panels
  if (meshA && isVec3Different(g_baseColorA, meshA->getSurfaceColor())) {
    g_baseColorA = meshA->getSurfaceColor();
  }
  if (meshB && isVec3Different(g_baseColorB, meshB->getSurfaceColor())) {
    g_baseColorB = meshB->getSurfaceColor();
  }

  // Checkbox toggles for Transform Gizmos
  if(ImGui::Checkbox("Enable Mesh A Gizmo", &g_enableGizmoA)) {
    meshA->setTransformGizmoEnabled(g_enableGizmoA);
  }
  ImGui::SameLine();
  if(ImGui::Checkbox("Enable Mesh B Gizmo", &g_enableGizmoB)) {
    meshB->setTransformGizmoEnabled(g_enableGizmoB);
  }

  if(ImGui::Button("Reset Gizmos to Identity")) {
    meshA->setTransform(glm::mat4(1.0f));
    meshB->setTransform(glm::mat4(1.0f));
  }

  ImGui::Separator();
  ImGui::TextUnformatted("Mesh Base Colors:");

  bool baseColorChanged = false;
  if(ImGui::ColorEdit3("Mesh A Surface", &g_baseColorA.x)) {
    baseColorChanged = true;
  }
  if(ImGui::ColorEdit3("Mesh B Surface", &g_baseColorB.x)) {
    baseColorChanged = true;
  }

  if(baseColorChanged) {
    updateMeshBaseColors();
    if(g_showIntersections && !g_lastExactPairs.empty()) {
      highlightIntersections(g_lastExactPairs, g_lastNumFacesA, g_lastNumFacesB);
    }
  }

  ImGui::Separator();
  ImGui::TextUnformatted("Intersection Highlights:");

  if(ImGui::Checkbox("Show Intersection Highlights", &g_showIntersections)) {
    if(!g_showIntersections) {
      if(meshA) meshA->removeQuantity("Intersections");
      if(meshB) meshB->removeQuantity("Intersections");
    } else if(!g_lastExactPairs.empty()) {
      highlightIntersections(g_lastExactPairs, g_lastNumFacesA, g_lastNumFacesB);
    }
  }

  if(g_showIntersections) {
    bool intersectColorChanged = false;

    const char* presetNames[] = {
        "Yellow (Default)", "Hot Red", "Neon Cyan", 
        "Magenta", "Electric Orange", "Lime Green", "White Hot"
    };

    // Preset Dropdown for Mesh A
    if(ImGui::Combo("Mesh A Highlight Preset", &g_selectedPresetA, presetNames, IM_ARRAYSIZE(presetNames))) {
      g_intersectColorA = g_presetColors[g_selectedPresetA].color;
      intersectColorChanged = true;
    }
    if(ImGui::ColorEdit3("Mesh A Custom Highlight", &g_intersectColorA.x)) {
      intersectColorChanged = true;
    }

    ImGui::Spacing();

    // Preset Dropdown for Mesh B
    if(ImGui::Combo("Mesh B Highlight Preset", &g_selectedPresetB, presetNames, IM_ARRAYSIZE(presetNames))) {
      g_intersectColorB = g_presetColors[g_selectedPresetB].color;
      intersectColorChanged = true;
    }
    if(ImGui::ColorEdit3("Mesh B Custom Highlight", &g_intersectColorB.x)) {
      intersectColorChanged = true;
    }

    if(intersectColorChanged && !g_lastExactPairs.empty()) {
      highlightIntersections(g_lastExactPairs, g_lastNumFacesA, g_lastNumFacesB);
    }
  }

  ImGui::Separator();

  // Fire Intersection Pipeline Button
  if(ImGui::Button("FIRE INTERSECTION PIPELINE", ImVec2(-1, 40))) {
    if(g_onFireCallback) {
      g_onFireCallback();
    }
  }

  ImGui::End();
}

// Initialize Polyscope engine
inline void init() {
  polyscope::options::autocenterStructures = false;
  polyscope::init();
  polyscope::state::userCallback = polyscopeUiCallback;
}

// Reset visual state & auto-fit camera
inline void reset(const Mesh& meshA,
                  const Mesh& meshB,
                  float3 centerA = make_float3(0, 0, 0),
                  float3 centerB = make_float3(0, 0, 0)) {
  g_currentTranslationB = glm::vec3(0.0f, 0.0f, 0.0f);
  g_centerA = centerA;
  g_centerB = centerB;

  g_lastExactPairs.clear();
  g_lastNumFacesA = 0;
  g_lastNumFacesB = 0;

  registerCgalMesh("Mesh A", meshA, g_baseColorA);
  registerCgalMesh("Mesh B", meshB, g_baseColorB);

  setGizmosEnabled(g_enableGizmoA, g_enableGizmoB);

  polyscope::view::resetCameraToHomeView();
}

inline void translateMeshB(float x, float y, float z) {
  g_currentTranslationB = glm::vec3(x, y, z);
  if(auto* meshB = polyscope::getSurfaceMesh("Mesh B")) {
    glm::vec3 currentRotB(0.0f), dummyTrans(0.0f);
    // Extract existing rotation relative to g_centerB
    decomposeTransformMatrix(meshB->getTransform(), toGlm(g_centerB), currentRotB, dummyTrans);
    // Rebuild matrix retaining active rotation
    meshB->setTransform(buildTransformMatrix(currentRotB, g_currentTranslationB, toGlm(g_centerB)));
  }
}

inline void transformMeshA(glm::vec3 rotDeg, glm::vec3 trans, glm::vec3 center = glm::vec3(0.0f)) {
  if(auto* meshA = polyscope::getSurfaceMesh("Mesh A")) {
    meshA->setTransform(buildTransformMatrix(rotDeg, trans, center));
  }
}

inline void transformMeshB(glm::vec3 rotDeg, glm::vec3 trans, glm::vec3 center = glm::vec3(0.0f)) {
  g_currentTranslationB = trans;
  if(auto* meshB = polyscope::getSurfaceMesh("Mesh B")) {
    meshB->setTransform(buildTransformMatrix(rotDeg, trans, center));
  }
}

inline void transformMeshA(float3 rotDeg, float3 trans, float3 center = make_float3(0, 0, 0)) {
  transformMeshA(toGlm(rotDeg), toGlm(trans), toGlm(center));
}

inline void transformMeshB(float3 rotDeg, float3 trans, float3 center = make_float3(0, 0, 0)) {
  transformMeshB(toGlm(rotDeg), toGlm(trans), toGlm(center));
}

inline void
transformBoth(float3 rotDegA, float3 transA, float3 centerA, float3 rotDegB, float3 transB, float3 centerB) {
  transformMeshA(rotDegA, transA, centerA);
  transformMeshB(rotDegB, transB, centerB);
}

// Reads the interactive gizmo transformation back into CUDA float3 vectors
inline bool getCurrentTransforms(float3& outRotA, float3& outTransA, float3& outRotB, float3& outTransB) {
  auto* meshA = polyscope::getSurfaceMesh("Mesh A");
  auto* meshB = polyscope::getSurfaceMesh("Mesh B");

  if(!meshA || !meshB)
    return false;

  glm::mat4 matA = meshA->getTransform();
  glm::mat4 matB = meshB->getTransform();

  glm::vec3 rotA, transA, rotB, transB;

  decomposeTransformMatrix(matA, toGlm(g_centerA), rotA, transA);
  decomposeTransformMatrix(matB, toGlm(g_centerB), rotB, transB);

  outRotA = make_float3(rotA.x, rotA.y, rotA.z);
  outTransA = make_float3(transA.x, transA.y, transA.z);
  outRotB = make_float3(rotB.x, rotB.y, rotB.z);
  outTransB = make_float3(transB.x, transB.y, transB.z);

  return true;
}

inline void highlightIntersections(const std::vector<int2>& exactPairs, size_t numFacesA, size_t numFacesB) {
  auto* meshA = polyscope::getSurfaceMesh("Mesh A");
  auto* meshB = polyscope::getSurfaceMesh("Mesh B");

  // Save parameters for live dynamic color tweaking in ImGui
  g_lastExactPairs = exactPairs;
  g_lastNumFacesA = numFacesA;
  g_lastNumFacesB = numFacesB;

  if (!g_showIntersections) return;

  // FIXED: Sync base color variables with actual Polyscope surface colors
  if (meshA) g_baseColorA = meshA->getSurfaceColor();
  if (meshB) g_baseColorB = meshB->getSurfaceColor();

  std::vector<glm::vec3> colorsA(numFacesA, g_baseColorA);
  std::vector<glm::vec3> colorsB(numFacesB, g_baseColorB);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, exactPairs.size()), [&](const tbb::blocked_range<size_t>& r) {
    for(size_t i = r.begin(); i != r.end(); ++i) {
      const auto& pair = exactPairs[i];
      if(pair.x >= 0 && static_cast<size_t>(pair.x) < numFacesA)
        colorsA[pair.x] = g_intersectColorA;
      if(pair.y >= 0 && static_cast<size_t>(pair.y) < numFacesB)
        colorsB[pair.y] = g_intersectColorB;
    }
  });

  if(meshA) {
    meshA->removeQuantity("Intersections");
    meshA->addFaceColorQuantity("Intersections", colorsA)->setEnabled(true);
  }
  if(meshB) {
    meshB->removeQuantity("Intersections");
    meshB->addFaceColorQuantity("Intersections", colorsB)->setEnabled(true);
  }
}

inline void drawFrame() {
  polyscope::frameTick();
}

} // namespace PolyscopeBridge