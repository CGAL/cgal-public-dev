#include "CPU/MeshIntersectionCurves.h"

#include <variant>
#include <optional>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <CGAL/intersections.h>
#include <CGAL/Aff_transformation_3.h>

#include "CPU/RotationTools.h"
#include "CPU/YellowFilter.h"

using Segment3 = Kernel::Segment_3;

namespace detail {

inline ExactPoint3POD toPOD(const Point3& pt) {
  return { pt.x(), pt.y(), pt.z() };
}

} // namespace detail

std::vector<IntersectionSegment> computeExactIntersectionCurves(
    const Mesh& meshA,
    const Mesh& meshB,
    const std::vector<int2>& candidatePairs,
    const Point3& centerA, const double3& rotDegA, const double3& transA,
    const Point3& centerB, const double3& rotDegB, const double3& transB)
{
  if (candidatePairs.empty()) return {};

  const bool isIdentityA = (rotDegA.x == 0.0 && rotDegA.y == 0.0 && rotDegA.z == 0.0 &&
                            transA.x == 0.0 && transA.y == 0.0 && transA.z == 0.0);
  const bool isIdentityB = (rotDegB.x == 0.0 && rotDegB.y == 0.0 && rotDegB.z == 0.0 &&
                            transB.x == 0.0 && transB.y == 0.0 && transB.z == 0.0);

  const auto xformA = isIdentityA ? CGAL::Aff_transformation_3<Kernel>(CGAL::IDENTITY) 
                                  : createRigidTransformation(centerA, rotDegA, transA);
  const auto xformB = isIdentityB ? CGAL::Aff_transformation_3<Kernel>(CGAL::IDENTITY) 
                                  : createRigidTransformation(centerB, rotDegB, transB);

  tbb::enumerable_thread_specific<std::vector<IntersectionSegment>> localBuffers;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, candidatePairs.size()),
    [&](const tbb::blocked_range<size_t>& range) {
      auto& localVec = localBuffers.local();
      localVec.reserve(range.size()); // Pre-allocate local buffer space

      for (size_t i = range.begin(); i != range.end(); ++i) {
        const int2 pair = candidatePairs[i];

        face_descriptor fA(static_cast<Mesh::size_type>(pair.x));
        face_descriptor fB(static_cast<Mesh::size_type>(pair.y));

        // Directly reuse detail::getTransformedTriangle from YellowFilter.h
        Triangle3 triA = detail::getTransformedTriangle(meshA, fA, xformA, isIdentityA);
        Triangle3 triB = detail::getTransformedTriangle(meshB, fB, xformB, isIdentityB);

        auto result = CGAL::intersection(triA, triB);
        if (!result) continue;

        if (const Segment3* seg = std::get_if<Segment3>(&*result)) {
          localVec.push_back({detail::toPOD(seg->source()), detail::toPOD(seg->target())});
        } 
        else if (const Point3* pt = std::get_if<Point3>(&*result)) {
          ExactPoint3POD pod = detail::toPOD(*pt);
          localVec.push_back({pod, pod});
        }
        else if (const Triangle3* tri = std::get_if<Triangle3>(&*result)) {
          localVec.push_back({detail::toPOD(tri->vertex(0)), detail::toPOD(tri->vertex(1))});
          localVec.push_back({detail::toPOD(tri->vertex(1)), detail::toPOD(tri->vertex(2))});
          localVec.push_back({detail::toPOD(tri->vertex(2)), detail::toPOD(tri->vertex(0))});
        }
      }
    });

  size_t totalSegments = 0;
  for (const auto& vec : localBuffers) totalSegments += vec.size();

  std::vector<IntersectionSegment> globalSegments;
  globalSegments.reserve(totalSegments);

  for (const auto& vec : localBuffers) {
    globalSegments.insert(globalSegments.end(), vec.begin(), vec.end());
  }

  return globalSegments;
}