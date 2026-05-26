cmake_minimum_required(VERSION 3.18...3.31) # Bumped to 3.18 for sm_120 architecture support
project(meshIntersectionRegistration LANGUAGES CXX CUDA)

# 1. Locate the CGAL build configuration from your parent build directory
find_package(CGAL REQUIRED PATHS "${CMAKE_CURRENT_SOURCE_DIR}/../../cgal-build" NO_DEFAULT_PATH)

# 2. Define path variables for cuBQL and CUDA context
set(CUBQL_ROOT "/home/yury/Projects/CPPProj/Cubql" CACHE PATH "Path to cuBQL source root")
set(CUBQL_BUILD "${CUBQL_ROOT}/build")

# 3. Find threading dependencies (Threads for general use, TBB for CGAL parallelization)
find_package(Threads REQUIRED)
find_package(TBB REQUIRED)

# ====================================================================
# TARGET 1: Your Original Project (main.cpp)
# ====================================================================
add_executable(meshIntersectionRegistration main.cpp kernels.cu)

set_target_properties(meshIntersectionRegistration PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_ARCHITECTURES "120"
    CUDA_SEPARABLE_COMPILATION ON
)

# Enable CGAL parallel features via TBB macro
target_compile_definitions(meshIntersectionRegistration PRIVATE CGAL_LINKED_WITH_TBB)

target_include_directories(meshIntersectionRegistration PRIVATE
    include
    ${CUBQL_ROOT}
    ${CUBQL_ROOT}/samples/common
    ${CUDA_TOOLKIT_ROOT_DIR}/include
)

target_link_libraries(meshIntersectionRegistration PRIVATE
    CGAL::CGAL
    gmp
    mpfr
    Threads::Threads
    TBB::tbb
    "${CUBQL_BUILD}/libcuBQL_cuda_float2.so"
    "${CUBQL_BUILD}/libcuBQL_cuda_float3.so"
)


# ====================================================================
# TARGET 2: Your New Massive Test Runner (main_massive.cpp)
# ====================================================================
add_executable(meshIntersectionBatch main_massive.cpp kernels.cu)

set_target_properties(meshIntersectionBatch PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_ARCHITECTURES "120"
    CUDA_SEPARABLE_COMPILATION ON
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${CUBQL_BUILD}" # Automatically fixes dynamic .so loading at runtime
)

# Enable CGAL parallel features via TBB macro
target_compile_definitions(meshIntersectionBatch PRIVATE CGAL_LINKED_WITH_TBB)

# Force high-performance release flags specifically for this target's files
target_compile_options(meshIntersectionBatch PRIVATE 
    $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG>
    $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>
)

target_include_directories(meshIntersectionBatch PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    include
    ${CUBQL_ROOT}
    ${CUBQL_ROOT}/samples/common
    ${CUDA_TOOLKIT_ROOT_DIR}/include
)

target_link_libraries(meshIntersectionBatch PRIVATE
    CGAL::CGAL
    gmp         # FIXED: Using lowercase system fallback
    mpfr        # FIXED: Using lowercase system fallback
    Threads::Threads
    TBB::tbb
    "${CUBQL_BUILD}/libcuBQL_cuda_float2.so"
    "${CUBQL_BUILD}/libcuBQL_cuda_float3.so"
)