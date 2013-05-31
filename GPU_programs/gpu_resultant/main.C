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

#include <stdlib.h>
#include <stdio.h>

#include <CGAL/GPU_algorithms/resultant_algorithm.h>

int main(int argc, char** argv) {

    CGAL::internal::GPU_resultant& obj =
             CGAL::internal::GPU_resultant::instance();
    obj.debug_run(argc, argv);
  
    return 0;
}

