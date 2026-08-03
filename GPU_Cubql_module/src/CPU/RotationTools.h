#ifndef ROTATION_TOOLS_H
#define ROTATION_TOOLS_H

#include <cmath>
#include <vector>
#include <algorithm>

// TBB Includes for CPU Parallelization
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

// Vector Types (CUDA host headers or standalone vector types defining double3/float3)
#include <vector_types.h>
#include <vector_functions.h>

// CGAL Definitions
#include "CgalDefinitions.h"

// Basic 3x3 Matrix Container for CPU Host Operations
struct Mat3x3 {
    double m[3][3];
};

// ============================================================================
// 1. ROTATION MATRIX GENERATOR (CPU)
// ============================================================================

inline Mat3x3 makeRotationMatrixDeg(double rxDeg, double ryDeg, double rzDeg) {
    constexpr double DEG_TO_RAD = 3.14159265358979323846 / 180.0;
    double radX = rxDeg * DEG_TO_RAD;
    double radY = ryDeg * DEG_TO_RAD;
    double radZ = rzDeg * DEG_TO_RAD;

    double cx = std::cos(radX), sx = std::sin(radX);
    double cy = std::cos(radY), sy = std::sin(radY);
    double cz = std::cos(radZ), sz = std::sin(radZ);

    Mat3x3 R;
    // Composite Rotation: R = Rz * Ry * Rx
    R.m[0][0] = cy * cz;
    R.m[0][1] = sx * sy * cz - cx * sz;
    R.m[0][2] = cx * sy * cz + sx * sz;

    R.m[1][0] = cy * sz;
    R.m[1][1] = sx * sy * sz + cx * cz;
    R.m[1][2] = cx * sy * sz - sx * cz;

    R.m[2][0] = -sy;
    R.m[2][1] = sx * cy;
    R.m[2][2] = cx * cy;

    return R;
}

inline Mat3x3 makeRotationMatrixDeg(const double3& rotDeg) {
    return makeRotationMatrixDeg(rotDeg.x, rotDeg.y, rotDeg.z);
}

// Helper to transform a single double3 point on CPU
inline double3 transformPoint(const double3& p, const Mat3x3& R, double3 center, double3 trans) {
    double x = p.x - center.x;
    double y = p.y - center.y;
    double z = p.z - center.z;

    double rx = R.m[0][0] * x + R.m[0][1] * y + R.m[0][2] * z;
    double ry = R.m[1][0] * x + R.m[1][1] * y + R.m[1][2] * z;
    double rz = R.m[2][0] * x + R.m[2][1] * y + R.m[2][2] * z;

    return make_double3(rx + center.x + trans.x, ry + center.y + trans.y, rz + center.z + trans.z);
}

// ============================================================================
// 2. CGAL AFFINE TRANSFORMATION BUILDER
// ============================================================================

inline CGAL::Aff_transformation_3<Kernel> createRigidTransformation(
    const Point3& center,
    const double3& rotDeg, 
    const double3& trans) 
{
    if (rotDeg.x == 0.0 && rotDeg.y == 0.0 && rotDeg.z == 0.0 &&
        trans.x == 0.0 && trans.y == 0.0 && trans.z == 0.0) 
    {
        return CGAL::Aff_transformation_3<Kernel>(CGAL::IDENTITY);
    }

    Mat3x3 R = makeRotationMatrixDeg(rotDeg);

    double cx_p = center.x();
    double cy_p = center.y();
    double cz_p = center.z();

    // T_total = C + T - (R * C)
    double tX = cx_p + trans.x - (R.m[0][0] * cx_p + R.m[0][1] * cy_p + R.m[0][2] * cz_p);
    double tY = cy_p + trans.y - (R.m[1][0] * cx_p + R.m[1][1] * cy_p + R.m[1][2] * cz_p);
    double tZ = cz_p + trans.z - (R.m[2][0] * cx_p + R.m[2][1] * cy_p + R.m[2][2] * cz_p);

    return CGAL::Aff_transformation_3<Kernel>(
        R.m[0][0], R.m[0][1], R.m[0][2], tX,
        R.m[1][0], R.m[1][1], R.m[1][2], tY,
        R.m[2][0], R.m[2][1], R.m[2][2], tZ
    );
}

// ============================================================================
// 3. METHOD 1: CGAL MESH TRANSFORMATION (CPU TBB)
// ============================================================================

// In-place mesh update relative to a original reference point list
inline void transformCgalMesh(
    Mesh* mesh, 
    const std::vector<Point3>& origPoints, 
    Point3 center, 
    double3 rotDeg, 
    double3 trans) 
{
    if (!mesh || origPoints.empty()) return;

    auto pmap = mesh->points();
    const size_t numPoints = mesh->num_vertices();
    Mat3x3 R = makeRotationMatrixDeg(rotDeg);

    double cx_p = center.x(), cy_p = center.y(), cz_p = center.z();
    double tx = trans.x, ty = trans.y, tz = trans.z;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, numPoints),
        [mesh, &pmap, &origPoints, R, cx_p, cy_p, cz_p, tx, ty, tz](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                typename Mesh::Vertex_index vd(static_cast<typename Mesh::size_type>(i));
                if (!mesh->is_removed(vd)) {
                    const Point3& pOrig = origPoints[i];
                    double x = pOrig.x() - cx_p;
                    double y = pOrig.y() - cy_p;
                    double z = pOrig.z() - cz_p;

                    double rx = R.m[0][0] * x + R.m[0][1] * y + R.m[0][2] * z;
                    double ry = R.m[1][0] * x + R.m[1][1] * y + R.m[1][2] * z;
                    double rz = R.m[2][0] * x + R.m[2][1] * y + R.m[2][2] * z;

                    pmap[vd] = Point3(rx + cx_p + tx, ry + cy_p + ty, rz + cz_p + tz);
                }
            }
        });
}

// Out-of-place mesh transformation overload (Mesh& meshIn -> Mesh& meshOut)
inline void transformCgalMesh(
    const Mesh& meshIn, 
    Mesh& meshOut, 
    Point3 center, 
    double3 rotDeg, 
    double3 trans) 
{
    meshOut = meshIn; // Copy topology
    
    std::vector<Point3> origPoints;
    origPoints.reserve(meshIn.num_vertices());
    auto pmapIn = meshIn.points();
    
    for (auto v : meshIn.vertices()) {
        origPoints.push_back(pmapIn[v]);
    }

    transformCgalMesh(&meshOut, origPoints, center, rotDeg, trans);
}

// ============================================================================
// 4. METHOD 2: CPU PARALLEL VERTEX BUFFER TRANSFORMATION (TBB)
// ============================================================================

inline void transformVerticesCPU(
    double3* dVertsOut, 
    const double3* dVertsOrig, 
    int numVerts, 
    Mat3x3 R, 
    double3 center, 
    double3 trans) 
{
    if (!dVertsOut || !dVertsOrig || numVerts <= 0) return;

    tbb::parallel_for(tbb::blocked_range<int>(0, numVerts),
        [dVertsOut, dVertsOrig, &R, center, trans](const tbb::blocked_range<int>& range) {
            for (int i = range.begin(); i != range.end(); ++i) {
                dVertsOut[i] = transformPoint(dVertsOrig[i], R, center, trans);
            }
        });
}

// Convenience overload accepting double3 rotation vector directly
inline void transformVerticesCPU(
    double3* dVertsOut, 
    const double3* dVertsOrig, 
    int numVerts, 
    double3 rotDeg, 
    double3 center, 
    double3 trans) 
{
    Mat3x3 R = makeRotationMatrixDeg(rotDeg);
    transformVerticesCPU(dVertsOut, dVertsOrig, numVerts, R, center, trans);
}

// ============================================================================
// 5. METHOD 3: CPU TYPE CONVERSION HELPERS (double3 -> float3)
// ============================================================================

// Converts a raw double3 CPU array to float3 in parallel using TBB
inline void convertDouble3ToFloat3(
    float3* outFloat, 
    const double3* inDouble, 
    size_t count) 
{
    if (!outFloat || !inDouble || count == 0) return;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, count),
        [outFloat, inDouble](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                outFloat[i] = make_float3(
                    static_cast<float>(inDouble[i].x),
                    static_cast<float>(inDouble[i].y),
                    static_cast<float>(inDouble[i].z)
                );
            }
        });
}

// std::vector helper overload
inline std::vector<float3> convertDouble3ToFloat3(const std::vector<double3>& inDouble) {
    std::vector<float3> outFloat(inDouble.size());
    convertDouble3ToFloat3(outFloat.data(), inDouble.data(), inDouble.size());
    return outFloat;
}

#endif // ROTATION_TOOLS_H