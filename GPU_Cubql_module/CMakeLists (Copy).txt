cmake_minimum_required(VERSION 3.15...3.31)
project(meshIntersectionRegistration LANGUAGES CXX CUDA)

# 1. Locate the CGAL build configuration from your parent build directory
# (We assume you will build CGAL in a directory named 'build' at the root)
find_package(CGAL REQUIRED PATHS "${CMAKE_CURRENT_SOURCE_DIR}/../../cgal-build" NO_DEFAULT_PATH)

# 2. Define path variables for cuBQL and CUDA context
set(CUBQL_ROOT "/home/yury/Projects/CPPProj/Cubql" CACHE PATH "Path to cuBQL source root")
set(CUBQL_BUILD "${CUBQL_ROOT}/build")

# 3. Define your target executable
add_executable(meshIntersectionRegistration main.cpp kernels.cu)

# 4. Set compilation and target options matching your RTX 5060 Ti architecture (sm_120)
set_target_properties(meshIntersectionRegistration PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_ARCHITECTURES "120"
    CUDA_SEPARABLE_COMPILATION ON
)

# 5. Provide target-specific include directories
target_include_directories(meshIntersectionRegistration PRIVATE
    include
    ${CUBQL_ROOT}
    ${CUBQL_ROOT}/samples/common
    ${CUDA_TOOLKIT_ROOT_DIR}/include
)

# 6. Bind dependencies (Links CGAL's dynamic header tree automatically)
target_link_libraries(meshIntersectionRegistration PRIVATE
    CGAL::CGAL
    gmp
    mpfr
    ${CUBQL_BUILD}/libcuBQL_cuda_float2.so
    ${CUBQL_BUILD}/libcuBQL_cuda_float3.so
)