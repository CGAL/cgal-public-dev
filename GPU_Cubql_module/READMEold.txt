============================================================
GSOC Project: CGAL GPU Integration (cuBQL Module)
============================================================

This module implements a GPU-accelerated mesh intersection 
algorithm using cuBQL for BVH traversal.

DEPENDENCIES:
- CUDA Toolkit (12.0+ recommended for Blackwell support)
- cuBQL: https://github.com/ingowald/cuBQL

COMPILATION:
By default, the Makefile looks for cuBQL in a local path. 
To compile on your system, provide the path to your cuBQL 
installation folder:

    make CUBQL_ROOT=/path/to/your/cuBQL

HARDWARE NOTE:
The Makefile is currently configured for sm_120 (RTX 50-series).
If you are using an older GPU, update the GENCODE variable 
in the Makefile to match your architecture (e.g., sm_86).

USAGE:
./meshIntersectionRegistration <meshA.off> <meshB.off> <output.csv>
