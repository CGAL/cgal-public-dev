#pragma once

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <vector_types.h> // For int2, double3

#include "CgalDefinitions.h"
#include "RotationTools.h"

namespace detail {
    // Internal helper to reconstruct triangle points from face descriptors
    inline Triangle3 getTransformedTriangle(
        const Mesh& mesh, 
        face_descriptor f, 
        const CGAL::Aff_transformation_3<Kernel>& xform,
        bool isIdentity) 
    {
        const auto& coords = mesh.points();
        auto h0 = mesh.halfedge(f);
        auto h1 = mesh.next(h0);
        auto h2 = mesh.next(h1);

        const auto& p0 = coords[mesh.target(h0)];
        const auto& p1 = coords[mesh.target(h1)];
        const auto& p2 = coords[mesh.target(h2)];

        if (isIdentity) {
            return Triangle3(p0, p1, p2);
        }
        return Triangle3(xform(p0), xform(p1), xform(p2));
    }
}

// Transformed / General Overload
inline void filterYellowPairsTBB(
    const Mesh& meshA,
    const Mesh& meshB,
    const int2* yellowPairsPtr,
    size_t numYellows,
    tbb::concurrent_vector<int2>& finalExactPairs,
    const Point3& centerA,
    const double3& rotDegA,
    const double3& transA,
    const Point3& centerB,
    const double3& rotDegB,
    const double3& transB) 
{
    if (numYellows == 0 || !yellowPairsPtr) return;

    const bool isIdentityA = (rotDegA.x == 0.0 && rotDegA.y == 0.0 && rotDegA.z == 0.0 &&
                              transA.x == 0.0 && transA.y == 0.0 && transA.z == 0.0);
    const bool isIdentityB = (rotDegB.x == 0.0 && rotDegB.y == 0.0 && rotDegB.z == 0.0 &&
                              transB.x == 0.0 && transB.y == 0.0 && transB.z == 0.0);

    const auto xformA = isIdentityA ? CGAL::Aff_transformation_3<Kernel>(CGAL::IDENTITY) 
                                    : createRigidTransformation(centerA, rotDegA, transA);
    const auto xformB = isIdentityB ? CGAL::Aff_transformation_3<Kernel>(CGAL::IDENTITY) 
                                    : createRigidTransformation(centerB, rotDegB, transB);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, numYellows), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            const int2 pair = yellowPairsPtr[i];

            face_descriptor fA(static_cast<Mesh::size_type>(pair.x));
            face_descriptor fB(static_cast<Mesh::size_type>(pair.y));

            Triangle3 triA = detail::getTransformedTriangle(meshA, fA, xformA, isIdentityA);
            Triangle3 triB = detail::getTransformedTriangle(meshB, fB, xformB, isIdentityB);

            if (CGAL::do_intersect(triA, triB)) {
                finalExactPairs.push_back(pair);
            }
        }
    });
}

// Convenient overload for non-transformed meshes
inline void filterYellowPairsTBB(
    const Mesh& meshA,
    const Mesh& meshB,
    const int2* yellowPairsPtr,
    size_t numYellows,
    tbb::concurrent_vector<int2>& finalExactPairs) 
{
    const Point3 zeroPoint(0.0, 0.0, 0.0);
    const double3 zeroVec{0.0, 0.0, 0.0};

    filterYellowPairsTBB(
        meshA, meshB, yellowPairsPtr, numYellows, finalExactPairs,
        zeroPoint, zeroVec, zeroVec,
        zeroPoint, zeroVec, zeroVec
    );
}