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
    
#ifndef _MOD_ENTRY_H_
#define _MOD_ENTRY_H_

namespace CGAL {

namespace internal {

struct MOD_entry {  //! modulus table entry
    unsigned m;     //! 24-/32-bit modulus m
    unsigned mu;    //! -m^(-1) mod 2^24 or mod 2^32 (depending on bitsize)
                    //! (needed for Montgomery inverse)
    union {
        struct {
            float invk;     //! 65536.0f / m
            float invm;     //! 1.0f / m
        };
        double invm_d;       //! 1.0 / m
    };
};

} //namespace internal

} // namespace CGAL 

#endif // _MOD_ENTRY_H_
