// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : AlciX
// File          : demos/webxalci/rasterizer.C
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:04 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#define NDEBUG 1

#include <iostream>
#include <vector>

// #define CGAL_NO_LEDA

#include "include/IPC.h"
#include "include/CGAL_includes.h"
#include "include/skeletonizer.h"

#include "include/SFA_interface.h"
#include "include/Surface_skeletonizer_impl.h"

extern pthread_mutex_t algorithms_mtx;

// Kernel_2 kernel_2_inst;

// NOTE NOTE: make sure both surf analysis and skeletonizer use the same
// kernel_2 instance: i.e. rely on get_static_instance() or smth..

// NOTE: possibly you can maintain a container of surf engines
// for each analysed surface..
static SFA_interface surf_engine;

extern pthread_mutex_t painter_mtx;
extern pthread_mutex_t arcs_cache_mtx;

Skeletonizer::Skeletonizer() {
}

//! analyzes new surface and makes it current in the surface engine
bool Skeletonizer::analyse_surface(const MD5_digest& surface_ID,
        const Polynomial_3& poly) {

    CGAL_precondition("switch off preconditions!!!");
    pthread_mutex_lock(&algorithms_mtx);
    bool ret = surf_engine.analyse_surface(surface_ID, poly);
    pthread_mutex_unlock(&algorithms_mtx);
    return ret;
}

//! loads already analysed surface 
bool Skeletonizer::load_surface(const MD5_digest& surface_ID) {

    CGAL_precondition("switch off preconditions!!!");
    pthread_mutex_lock(&algorithms_mtx);
    bool ret = surf_engine.load_surface(surface_ID);
    pthread_mutex_unlock(&algorithms_mtx);
    return ret;
}

//! computes surface triangulation
//! NOTE NOTE NOTE: \c data and \c triangle_data point to the same
//! location therefore parameters from \c data must be copied in the first
//! place to avoid overwriting !!
uint Skeletonizer::skeletonize(const MD5_digest& surface_ID,
    SHM_Data *data, void *triangle_data, uint max_data_sz,
         bool use_auto_bounds) {

    // NOTE: Surface_skeletonizer_impl should probably be put inside
// of Skeletonizer: one time use only
    CGAL::Surface_skeletonizer_impl< Kernel_2, SFA_interface >
        skeletonizer_impl(0);

    // NOTE skeletonizer parameters must be checked for sanity beforehand !!
    skeletonizer_impl.reset_parameters(data->sx, data->sy, data->en_left,
           data->en_right, data->en_bottom, data->en_top,
           data->z_below, data->z_above);

    pthread_mutex_lock(&algorithms_mtx);
    skeletonizer_impl.set_surface_engine(&surf_engine);
    skeletonizer_impl.skeletonize(use_auto_bounds);
    pthread_mutex_unlock(&algorithms_mtx);

    const Triangles_vector& tris = skeletonizer_impl.get_triangles();
    const Point_vec_3d& vd = skeletonizer_impl.get_vertices();

    Point_vec_3d nd;
    Point_3d centroid;
    bool inverse_normals = false;
    skeletonizer_impl.compute_normals(vd, nd, centroid, inverse_normals);
    
    const Triangles_vector& indices = skeletonizer_impl.get_triangles();
    const Triangles_vector& sil_indices =
            skeletonizer_impl.get_sil_triangles();
    const std::vector< unsigned >& zpatches =
        skeletonizer_impl.get_zlift_ranges();

    double z_silhouette = skeletonizer_impl.get_z_silhouette();
    
    uint n_verts = vd.size(), n_tris = indices.size(),
        n_sil_tris = sil_indices.size();
    CGAL_assertion(vd.size() == nd.size());

    Triangulation_info *pinfo = (Triangulation_info *)triangle_data;
    pinfo->n_verts = n_verts;
    pinfo->n_tris = n_tris;
    pinfo->n_sil_tris = n_sil_tris;
    pinfo->n_z_patches = zpatches.size();
    pinfo->cx = centroid[0], pinfo->cy = centroid[1], pinfo->cz = centroid[2];

    skeletonizer_impl.get_surface_bounds(pinfo->left,
           pinfo->right, pinfo->btm, pinfo->top);

    uint *puint = (uint *)(pinfo + 1);
    std::copy(zpatches.begin(), zpatches.end(), puint);
    float *pfloat = (float *)(puint + zpatches.size());

    // estimate the size of buffer required
    uint sz = sizeof(Triangulation_info) + zpatches.size()*sizeof(uint) +
            n_verts * 2 * 3 * sizeof(float) +
            (n_tris * 3 + n_sil_tris * 3) * sizeof(unsigned short);
    if(sz > max_data_sz) {
        fprintf(stderr, "Insufficient space for triangulation data: %d\n", sz);
        return -1u;
    }

    // clamp data to single precision before loading to VBO
    float *verts = pfloat;
    for(uint i = 0; i < n_verts; i++) {
        pfloat[0] = static_cast<float>(vd[i][0]);
        pfloat[1] = static_cast<float>(vd[i][1]);
        pfloat[2] = static_cast<float>(vd[i][2]);
        pfloat += 3;
    }
    fprintf(stderr, "Centroid: (%.5f %.5f %.5f)\n", pinfo->cx, pinfo->cy,
                pinfo->cz);

    float *normals = pfloat;
    for(uint i = 0; i < n_verts; i++) {
        pfloat[0] = static_cast<float>(nd[i][0]);
        pfloat[1] = static_cast<float>(nd[i][1]);
        pfloat[2] = static_cast<float>(nd[i][2]);
        pfloat += 3;
    }
    
    unsigned short *pidx_short = (unsigned short *)pfloat;
    for(uint i = 0; i < n_tris; i++) {
        pidx_short[0] = indices[i][0];
        pidx_short[1] = indices[i][1];
        pidx_short[2] = indices[i][2];
        pidx_short += 3;
    }

    for(uint i = 0; i < n_sil_tris; i++) {
        pidx_short[0] = sil_indices[i][0];
        pidx_short[1] = sil_indices[i][1];
        pidx_short[2] = sil_indices[i][2];
        pidx_short += 3;
    }
    
    fprintf(stderr, "Model (%d bytes) n_verts: %d; n_tris: %d; "
            "sil_tris: %d; z_patches: %d\n", sz, n_verts, n_tris, n_sil_tris,
                zpatches.size());
    return sz;
}
