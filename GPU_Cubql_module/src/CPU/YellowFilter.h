#pragma once

#include <cmath>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <vector_types.h> // For int2, double3

#include "CgalDefinitions.h"

// Helper function to build a CGAL Affine Transformation matching transformCgalMesh strategy:
// Degrees -> Radians conversion, rotation around custom center, and translation.
inline CGAL::Aff_transformation_3<Kernel> createRigidTransformation(
    const Point3& center,
    const double3& rotDeg, 
    const double3& trans) 
{
    constexpr double DEG_TO_RAD = 3.14159265358979323846 / 180.0;
    double radX = rotDeg.x * DEG_TO_RAD;
    double radY = rotDeg.y * DEG_TO_RAD;
    double radZ = rotDeg.z * DEG_TO_RAD;

    double cx = std::cos(radX), sx = std::sin(radX);
    double cy = std::cos(radY), sy = std::sin(radY);
    double cz = std::cos(radZ), sz = std::sin(radZ);

    // Rotation Matrix R (Rz * Ry * Rx)
    double r00 = cy * cz;
    double r01 = sx * sy * cz - cx * sz;
    double r02 = cx * sy * cz + sx * sz;

    double r10 = cy * sz;
    double r11 = sx * sy * sz + cx * cz;
    double r12 = cx * sy * sz - sx * cz;

    double r20 = -sy;
    double r21 = sx * cy;
    double r22 = cx * cy;

    double cx_p = center.x();
    double cy_p = center.y();
    double cz_p = center.z();

    // Mathematically equivalent to: R * (P - C) + C + T
    // Total translation component: T_total = C + T - (R * C)
    double tX = cx_p + trans.x - (r00 * cx_p + r01 * cy_p + r02 * cz_p);
    double tY = cy_p + trans.y - (r10 * cx_p + r11 * cy_p + r12 * cz_p);
    double tZ = cz_p + trans.z - (r20 * cx_p + r21 * cy_p + r22 * cz_p);

    // Construct 3x4 CGAL Affine Transformation matrix
    return CGAL::Aff_transformation_3<Kernel>(
        r00, r01, r02, tX,
        r10, r11, r12, tY,
        r20, r21, r22, tZ
    );
}

// Original non-transformed function
inline void filterYellowPairsTBB(
    const Mesh& meshA,
    const Mesh& meshB,
    const int2* yellowPairsPtr,
    size_t numYellows,
    tbb::concurrent_vector<int2>& finalExactPairs) 
{
    if (numYellows == 0 || !yellowPairsPtr) return;

    const auto& coords1 = meshA.points();
    const auto& coords2 = meshB.points();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, numYellows), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            const int2 pair = yellowPairsPtr[i];

            face_descriptor fA(static_cast<Mesh::size_type>(pair.x));
            face_descriptor fB(static_cast<Mesh::size_type>(pair.y));

            auto halfedgeA = meshA.halfedge(fA);
            auto vA0 = meshA.target(halfedgeA);
            auto vA1 = meshA.target(meshA.next(halfedgeA));
            auto vA2 = meshA.target(meshA.next(meshA.next(halfedgeA)));
            Triangle3 triA(coords1[vA0], coords1[vA1], coords1[vA2]);

            auto halfedgeB = meshB.halfedge(fB);
            auto vB0 = meshB.target(halfedgeB);
            auto vB1 = meshB.target(meshB.next(halfedgeB));
            auto vB2 = meshB.target(meshB.next(meshB.next(halfedgeB)));
            Triangle3 triB(coords2[vB0], coords2[vB1], coords2[vB2]);

            if (CGAL::do_intersect(triA, triB)) {
                finalExactPairs.push_back(pair);
            }
        }
    });
} 

// Overloaded function accepting 6 parameters per mesh (Center, Rotation in Degrees, Translation)
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
    // Fast-path: If all rotation and translation parameters are 0, delegate to original method
    bool isIdentityA = (rotDegA.x == 0.0 && rotDegA.y == 0.0 && rotDegA.z == 0.0 &&
                        transA.x == 0.0 && transA.y == 0.0 && transA.z == 0.0);
    bool isIdentityB = (rotDegB.x == 0.0 && rotDegB.y == 0.0 && rotDegB.z == 0.0 &&
                        transB.x == 0.0 && transB.y == 0.0 && transB.z == 0.0);

    if (isIdentityA && isIdentityB) {
        filterYellowPairsTBB(meshA, meshB, yellowPairsPtr, numYellows, finalExactPairs);
        return;
    }

    if (numYellows == 0 || !yellowPairsPtr) return;

    // Pre-calculate transformation matrices once outside the parallel loop
    auto xformA = createRigidTransformation(centerA, rotDegA, transA);
    auto xformB = createRigidTransformation(centerB, rotDegB, transB);

    const auto& coords1 = meshA.points();
    const auto& coords2 = meshB.points();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, numYellows), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            const int2 pair = yellowPairsPtr[i];

            face_descriptor fA(static_cast<Mesh::size_type>(pair.x));
            face_descriptor fB(static_cast<Mesh::size_type>(pair.y));

            // Extract vertices and conditionally transform Mesh A triangle
            auto halfedgeA = meshA.halfedge(fA);
            auto vA0 = meshA.target(halfedgeA);
            auto vA1 = meshA.target(meshA.next(halfedgeA));
            auto vA2 = meshA.target(meshA.next(meshA.next(halfedgeA)));

            Point3 pA0 = isIdentityA ? coords1[vA0] : xformA(coords1[vA0]);
            Point3 pA1 = isIdentityA ? coords1[vA1] : xformA(coords1[vA1]);
            Point3 pA2 = isIdentityA ? coords1[vA2] : xformA(coords1[vA2]);
            Triangle3 triA(pA0, pA1, pA2);

            // Extract vertices and conditionally transform Mesh B triangle
            auto halfedgeB = meshB.halfedge(fB);
            auto vB0 = meshB.target(halfedgeB);
            auto vB1 = meshB.target(meshB.next(halfedgeB));
            auto vB2 = meshB.target(meshB.next(meshB.next(halfedgeB)));

            Point3 pB0 = isIdentityB ? coords2[vB0] : xformB(coords2[vB0]);
            Point3 pB1 = isIdentityB ? coords2[vB1] : xformB(coords2[vB1]);
            Point3 pB2 = isIdentityB ? coords2[vB2] : xformB(coords2[vB2]);
            Triangle3 triB(pB0, pB1, pB2);

            if (CGAL::do_intersect(triA, triB)) {
                finalExactPairs.push_back(pair);
            }
        }
    });
}