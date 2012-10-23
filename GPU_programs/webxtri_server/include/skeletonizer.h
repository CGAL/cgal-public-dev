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
// File          : demos/webxalci/include/rasterizer.h
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:11 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*! \file skeletonizer.h
 *  \brief Defines class \c Skeletonizer 
 *  
 *  Algebraic surface analysis and triangulation routines for \c Xtri_server
 */

#ifndef XTRI_SKELETONIZER_H
#define XTRI_SKELETONIZER_H

//! \class Skeletonizer
//! provides the server with an interface to EXACUS curve analysis and
//! rendering
class Skeletonizer 
{
public:
    //!\name public methods
    //!@{    

    //! default constructor
    Skeletonizer();
    
    //! analyzes new surface and makes it current in the surface engine
    bool analyse_surface(const MD5_digest& surface_ID,
        const Polynomial_3& poly);

    //! loads surface from the cache
    bool load_surface(const MD5_digest& surface_ID);

    //! computes surface triangulation
    uint skeletonize(const MD5_digest& surface_ID,
            SHM_Data *data, void *triangle_data, uint max_data_sz,
            bool use_auto_bounds);

    //!@}
private:
    //!\name private data members
    //!@{
   
    //!@}
};

#endif // XTRI_SKELETONIZER_H

