#include "CPU/MeshIntersectionCurves.h"

#include <variant>
#include <optional>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/Aff_transformation_3.h>

#include "CPU/RotationTools.h"

using ExactKernel = CGAL::Exact_predicates_exact_constructions_kernel;
using ExactPoint3 = ExactKernel::Point_3;
using ExactTriangle3 = ExactKernel::Triangle_3;
using ExactSegment3 = ExactKernel::Segment_3;

namespace detail {

inline ExactTriangle3 getTransformedExactTriangle(
    const Mesh& mesh,
    face_descriptor f,
    const CGAL::Aff_transformation_3<Kernel>& xform,
    bool isIdentity)
{
  auto h0 = mesh.halfedge(f);
  auto h1 = mesh.next(h0);
  auto h2 = mesh.next(h1);

  const auto& coords = mesh.points();
  
  const Point3& p0 = coords[mesh.target(h0)];
  const Point3& p1 = coords[mesh.target(h1)];
  const Point3& p2 = coords[mesh.target(h2)];

  Point3 tp0 = isIdentity ? p0 : xform(p0);
  Point3 tp1 = isIdentity ? p1 : xform(p1);
  Point3 tp2 = isIdentity ? p2 : xform(p2);

  return ExactTriangle3(
      ExactPoint3(tp0.x(), tp0.y(), tp0.z()),
      ExactPoint3(tp1.x(), tp1.y(), tp1.z()),
      ExactPoint3(tp2.x(), tp2.y(), tp2.z())
  );
}

inline ExactPoint3POD toPOD(const ExactPoint3& pt) {
  return { CGAL::to_double(pt.x()), CGAL::to_double(pt.y()), CGAL::to_double(pt.z()) };
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

      for (size_t i = range.begin(); i != range.end(); ++i) {
        const int2 pair = candidatePairs[i];

        face_descriptor fA(static_cast<Mesh::size_type>(pair.x));
        face_descriptor fB(static_cast<Mesh::size_type>(pair.y));

        ExactTriangle3 triA = detail::getTransformedExactTriangle(meshA, fA, xformA, isIdentityA);
        ExactTriangle3 triB = detail::getTransformedExactTriangle(meshB, fB, xformB, isIdentityB);

        auto result = CGAL::intersection(triA, triB);
        if (!result) continue;

        if (const ExactSegment3* seg = std::get_if<ExactSegment3>(&*result)) {
          localVec.push_back({detail::toPOD(seg->source()), detail::toPOD(seg->target())});
        } 
        else if (const ExactPoint3* pt = std::get_if<ExactPoint3>(&*result)) {
          ExactPoint3POD pod = detail::toPOD(*pt);
          localVec.push_back({pod, pod});
        }
        else if (const ExactTriangle3* tri = std::get_if<ExactTriangle3>(&*result)) {
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