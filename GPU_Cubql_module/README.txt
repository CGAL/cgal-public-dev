## 🚀 High-Performance Build Instructions (Crucial for Evaluation)

This module implements a hybrid GPU/CPU pipeline utilizing **cuBQL** for highly parallel bounding volume hierarchy (BVH) broad-phase pruning, and **CGAL (EPICK)** for exact geometric predicate resolution. 

Because CGAL relies heavily on arithmetic interval filtering and template-heavy multi-precision arithmetic, compiling in **Debug mode will result in a 10x–30x performance degradation**. To evaluate the true speed of the pipeline ($\approx$680ms for 200k-face meshes), the project **must** be built with maximum optimizations and safety assertions disabled (`-O3 -DNDEBUG`).

### Prerequisites
* CUDA Toolkit 12.x+
* CGAL 5.x+ 
* GMP & MPFR libraries

### Boosted Compilation Sequence

To preserve workspace cleanliness, use CMake's explicit source (`-S`) and build (`-B`) target flags. This ensures all compiler optimizations, loop unrolling, and native vector lanes are fully utilized:

```bash
# 1. Configure the project directly into the build directory with Release flags enabled
cmake -S ./cgal-public-dev/GPU_Cubql_module \
      -B ./cgal-build/GPU_Cubql_module \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS="-O3 -march=native -DNDEBUG"

# 2. Compile using all available CPU threads
cmake --build ./cgal-build/GPU_Cubql_module -j$(nproc)
