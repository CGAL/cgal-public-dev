// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// this file is not part of any library ;-)
//
// ----------------------------------------------------------------------------
//
// Library       : CUDA MP
//
// File          : 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================
    
#ifndef _COMMON_DEFS_H_
#define _COMMON_DEFS_H_

namespace CGAL {

namespace internal {

typedef std::vector< __mpz_struct > MPZ_vector_1;
typedef std::vector< MPZ_vector_1 > MPZ_vector_2;

class GPU_algorithm_exception {
};

} //namespace internal

} // namespace CGAL 

#endif // _COMMON_DEFS_H_
