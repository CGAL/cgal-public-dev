#pragma once

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <vector_types.h> // For int2

#include "CgalDefinitions.h"

// Explicit non-templated function using hardcoded types from CgalDefinitions.h
// FIX: Added 'inline' to prevent "multiple definition" errors when included in both .cpp and .cu files
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