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

#ifndef _CONST_MEM_DEFS_H_
#define _CONST_MEM_DEFS_H_

enum ConstMemEntries { // change this to 'data_padding'
    NU = 0, NV, NU_OFS4, NV_OFS4,
    MAX_NU, DATA_PAD, GCD_SZ, N_BLOCKS,
    DATA_IN_, DATA_OUT_,
    MODS_PAD
};

const unsigned ConstMemDataSz = 11;

#endif // _CONST_MEM_DEFS_H_
