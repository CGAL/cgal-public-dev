#pragma once

namespace Tables {
    /*
     * Naming convention from "A parallel dual marching cubes approach to quad only surface reconstruction - Grosso & Zint"
     *
     *        ^ y
     *        |
     *       v2------e2------v3
     *       /|             /|
     *     e11|           e10|
     *     /  e3          /  e1
     *   v6------e6------v7  |
     *    |   |          |   |
     *    |  v0------e0--|---v1 --> x
     *    e7 /           e5 /
     *    | e8           | e9
     *    |/             |/
     *   v4------e4------v5
     *   /
     *  < z
     */

    constexpr int N_VERTICES = 8;
    constexpr int N_EDGES    = 12;

    // This table iterates around an edge of a voxel in positive direction, starting from the given voxel (0,0,0). The iteration is described in
    // coordinates relative to the given voxel. The last number is the local edge index.
    constexpr int edge_to_voxel_neighbor[N_EDGES][4][4] = {
        { { 0, 0, 0, 0 }, { 0, -1, 0, 2 }, { 0, -1, -1, 6 }, { 0, 0, -1, 4 } },      // e0
        { { 0, 0, 0, 1 }, { 1, 0, 0, 3 }, { 1, 0, -1, 7 }, { 0, 0, -1, 5 } },        // e1
        { { 0, 0, 0, 2 }, { 0, 0, -1, 6 }, { 0, 1, -1, 4 }, { 0, 1, 0, 0 } },        // e2
        { { 0, 0, 0, 3 }, { 0, 0, -1, 7 }, { -1, 0, -1, 5 }, { -1, 0, 0, 1 } },      // e3
        { { 0, 0, 0, 4 }, { 0, 0, 1, 0 }, { 0, -1, 1, 2 }, { 0, -1, 0, 6 } },        // e4
        { { 0, 0, 0, 5 }, { 0, 0, 1, 1 }, { 1, 0, 1, 3 }, { 1, 0, 0, 7 } },          // e5
        { { 0, 0, 0, 6 }, { 0, 1, 0, 4 }, { 0, 1, 1, 0 }, { 0, 0, 1, 2 } },          // e6
        { { 0, 0, 0, 7 }, { -1, 0, 0, 5 }, { -1, 0, 1, 1 }, { 0, 0, 1, 3 } },        // e7
        { { 0, 0, 0, 8 }, { -1, 0, 0, 9 }, { -1, -1, 0, 10 }, { 0, -1, 0, 11 } },    // e8
        { { 0, 0, 0, 9 }, { 0, -1, 0, 10 }, { 1, -1, 0, 11 }, { 1, 0, 0, 8 } },      // e9
        { { 0, 0, 0, 10 }, { 1, 0, 0, 11 }, { 1, 1, 0, 8 }, { 0, 1, 0, 9 } },        // e10
        { { 0, 0, 0, 11 }, { 0, 1, 0, 8 }, { -1, 1, 0, 9 }, { -1, 0, 0, 10 } }       // e11
    };

    /* The global edge index consists of the lexicographical index of the v0 vertex of a voxel, and an index that represents the axis. This table maps
     * from the axis index to the local edge index:
     * 0 = x-axis --> 0
     * 1 = y-axis --> 3
     * 2 = z-axis --> 8
     */
    constexpr int edge_store_index[3] = { 0, 3, 8 };

    // The local vertex indices of an edge. The indices are sorted by axis direction.
    constexpr int edge_to_vertex[N_EDGES][2] = {
        { 0, 1 },    // e0
        { 1, 3 },    // e1
        { 2, 3 },    // e2
        { 0, 2 },    // e3
        { 4, 5 },    // e4
        { 5, 7 },    // e5
        { 6, 7 },    // e6
        { 4, 6 },    // e7
        { 0, 4 },    // e8
        { 1, 5 },    // e9
        { 3, 7 },    // e10
        { 2, 6 }     // e11
    };

    // The local vertex coordinates within a voxel.
    constexpr int local_vertex_position[N_VERTICES][3] = {
        { 0, 0, 0 },    // v0
        { 1, 0, 0 },    // v1
        { 0, 1, 0 },    // v2
        { 1, 1, 0 },    // v3
        { 0, 0, 1 },    // v4
        { 1, 0, 1 },    // v5
        { 0, 1, 1 },    // v6
        { 1, 1, 1 }     // v7
    };

}    // namespace Tables