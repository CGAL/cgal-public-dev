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
// Library       : QdX
// File          : demos/xsurface/math.C
// QdX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*!@file math.C
 * collects various math routines
 */

#include "include/math_routines.h"

void quaternion_from_axis(Quaternion& q, 
            const Vector_3d& axis, double angle) {
    angle *= M_PI / 360.0;
    q.w = cos(angle);
    q.v = axis; 
    q.copy_on_write();
    q.v.normalize(); 
    q.v *= sin(angle);
}
  
void quaternion_trackball(Quaternion& q, const Point_2d& p1,
        const Point_2d& p2) {

    double d = p1[0]*p1[0] + p1[1]*p1[1];
    Vector_3d v1(p1[0], p1[1], 0), v2(p2[0], p2[1], 0);

    if(d > 1.0) {
        d = sqrt(d);
        v1 /= d;
    } else 
        v1[2] = sqrt(1.0 - d);

    d = p2[0]*p2[0] + p2[1]*p2[1];
    if(d > 1.0) {
        d = sqrt(d);
        v2 /= d;
    } else 
        v2[2] = sqrt(1.0 - d);
    
    q.copy_on_write();
    q.v = -v1.cross_product(v2);
    q.w = 1.0 + v1.dot_product(v2);    
}

// converts a quaternion to 4x4 matrix
void quaternion_2_matrix(const Quaternion& q, double *mat) {
   
    mat[0] = 1.0 - 2.0*(q.v[1]*q.v[1]+q.v[2]*q.v[2]);
    mat[1] = 2.0*(q.v[0]*q.v[1]-q.v[2]*q.w);
    mat[2] = 2.0*(q.v[2]*q.v[0]+q.v[1]*q.w);
    mat[3] = 0.0;
    
    mat[4] = 2.0*(q.v[0]*q.v[1]+q.v[2]*q.w);
    mat[5] = 1.0-2.0*(q.v[2]*q.v[2]+q.v[0]*q.v[0]);
    mat[6] = 2.0*(q.v[1]*q.v[2]-q.v[0]*q.w);
    mat[7] = 0.0;
    
    mat[8] = 2.0*(q.v[2]*q.v[0]-q.v[1]*q.w);
    mat[9] = 2.0*(q.v[1]*q.v[2]+q.v[0]*q.w);
    mat[10] = 1.0-2.0*(q.v[1]*q.v[1]+q.v[0]*q.v[0]);
    mat[11] = 0.0;
    
    mat[12] = 0.0;
    mat[13] = 0.0;
    mat[14] = 0.0;
    mat[15] = 1.0;
}

void quaternion_normalize(Quaternion& q) {
 
    q.copy_on_write();
    double mag = sqrt(q.v.square_magnitude() + q.w*q.w);
    q.v /= mag;
    q.w /= mag;
}

void quaternion_multiply(Quaternion& res, 
        const Quaternion& p, const Quaternion& q) {
    
    res.copy_on_write();
    double w = p.w * q.w - p.v.dot_product(q.v);
    Vector_3d tmp = p.v.cross_product(q.v);
    res.v = tmp + q.v * p.w + p.v * q.w;
    res.w = w;
}
