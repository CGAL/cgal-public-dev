#pragma once

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <charconv>
#include <cmath>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

// TBB Includes
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

// CUDA Vector Types (for double3 and uint3)
#include <vector_types.h>
#include "CgalDefinitions.h"

// CGAL Surface Mesh
//#include <CGAL/Surface_mesh.h>

namespace ParallelIO {

/**
 * @brief Memory-maps an OFF file and parses vertex coordinates and triangle face topology
 *        directly into host vector buffers in parallel using TBB with full 64-bit double precision.
 */
template <typename Point3>
inline bool loadOffFastTBB(const std::string& filepath,
                           std::vector<double3>& outVerts,
                           std::vector<uint3>& outIndices) 
{
  int fd = open(filepath.c_str(), O_RDONLY);
  if (fd == -1) return false;

  struct stat sb;
  if (fstat(fd, &sb) == -1) {
    close(fd);
    return false;
  }
  size_t fileSize = sb.st_size;

  // 1. POSIX Memory-map the file (Zero-copy disk access)
  char* data = static_cast<char*>(mmap(NULL, fileSize, PROT_READ, MAP_PRIVATE, fd, 0));
  close(fd);
  if (data == MAP_FAILED) return false;

  std::string_view fileView(data, fileSize);
  size_t pos = 0;

  auto nextLine = [&](size_t start) {
    size_t end = fileView.find('\n', start);
    return (end == std::string_view::npos) ? fileSize : end;
  };

  // 2. Read OFF Header
  size_t lineEnd = nextLine(pos);
  if (fileView.substr(pos, lineEnd - pos).find("OFF") == std::string_view::npos) {
    munmap(data, fileSize);
    return false;
  }
  pos = lineEnd + 1;

  // Skip comments and empty lines
  while (pos < fileSize && (fileView[pos] == '#' || fileView[pos] == '\n' || fileView[pos] == '\r')) {
    pos = nextLine(pos) + 1;
  }

  // 3. Read Element Counts (numVerts, numFaces)
  size_t numVerts = 0, numFaces = 0;
  lineEnd = nextLine(pos);
  std::string_view countLine = fileView.substr(pos, lineEnd - pos);
  const char* ptr = countLine.data();
  const char* endPtr = ptr + countLine.size();

  ptr = std::from_chars(ptr, endPtr, numVerts).ptr + 1;
  ptr = std::from_chars(ptr, endPtr, numFaces).ptr;
  pos = lineEnd + 1;

  outVerts.resize(numVerts);
  outIndices.resize(numFaces);

  // 4. Pre-scan line start offsets
  std::vector<size_t> vertLineOffsets(numVerts);
  for (size_t i = 0; i < numVerts; ++i) {
    vertLineOffsets[i] = pos;
    pos = nextLine(pos) + 1;
  }

  std::vector<size_t> faceLineOffsets(numFaces);
  for (size_t i = 0; i < numFaces; ++i) {
    faceLineOffsets[i] = pos;
    pos = nextLine(pos) + 1;
  }

  // 5. Parallel Vertex Parsing (Direct double precision parsing)
  constexpr size_t GRAIN_SIZE = 4096;
  tbb::parallel_for(tbb::blocked_range<size_t>(0, numVerts, GRAIN_SIZE),
    [&](const tbb::blocked_range<size_t>& r) {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        const char* p = data + vertLineOffsets[i];
        double x = 0.0, y = 0.0, z = 0.0;

        p = std::from_chars(p, data + fileSize, x).ptr + 1;
        p = std::from_chars(p, data + fileSize, y).ptr + 1;
        std::from_chars(p, data + fileSize, z);

        outVerts[i] = double3{x, y, z};
      }
    });

  // 6. Parallel Face Topology Parsing
  tbb::parallel_for(tbb::blocked_range<size_t>(0, numFaces, GRAIN_SIZE),
    [&](const tbb::blocked_range<size_t>& r) {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        const char* p = data + faceLineOffsets[i];
        int nv = 0;
        uint32_t v0 = 0, v1 = 0, v2 = 0;

        p = std::from_chars(p, data + fileSize, nv).ptr + 1; // skip vertex count (e.g. '3')
        p = std::from_chars(p, data + fileSize, v0).ptr + 1;
        p = std::from_chars(p, data + fileSize, v1).ptr + 1;
        std::from_chars(p, data + fileSize, v2);

        outIndices[i] = make_uint3(v0, v1, v2);
      }
    });

  munmap(data, fileSize);
  return true;
}

/**
 * @brief Bulk-populates a CGAL::Surface_mesh in memory from flat pre-parsed double-precision
 *        vertex and index arrays.
 */
template <typename Mesh, typename Point3>
inline void buildCgalMeshFromBuffers(const std::vector<double3>& verts,
                                     const std::vector<uint3>& indices,
                                     Mesh& mesh) 
{
  mesh.clear();
  mesh.reserve(verts.size(), indices.size() * 3, indices.size());

  std::vector<typename Mesh::Vertex_index> vhandles(verts.size());
  for (size_t i = 0; i < verts.size(); ++i) {
    vhandles[i] = mesh.add_vertex(Point3(verts[i].x, verts[i].y, verts[i].z));
  }

  for (size_t i = 0; i < indices.size(); ++i) {
    const auto& tri = indices[i];
    mesh.add_face(vhandles[tri.x], vhandles[tri.y], vhandles[tri.z]);
  }
}

template <typename Mesh, typename Point3>
inline bool loadOffToCgalMesh(const std::string& filepath,
                              Mesh& mesh,
                              std::vector<uint3>& outIndices) 
{
  std::vector<double3> tempDoubleVerts;
  if (!loadOffFastTBB<Point3>(filepath, tempDoubleVerts, outIndices)) {
    return false;
  }
  buildCgalMeshFromBuffers<Mesh, Point3>(tempDoubleVerts, outIndices, mesh);
  return true; // tempDoubleVerts is deallocated right here!
}

} // namespace ParallelIO