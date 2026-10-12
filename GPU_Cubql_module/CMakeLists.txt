cmake_minimum_required(VERSION 3.18...3.31)
project(meshIntersectionRegistration LANGUAGES CXX CUDA)

# ====================================================================
# GLOBAL OPTIONS & BUILD CONFIGURATION
# ====================================================================

option(BUILD_ARCHIVED_TARGETS "Build archived legacy executables" OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile commands for IDEs" FORCE)

# Global language standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# ====================================================================
# DEPENDENCIES & PACKAGE LOCATIONS
# ====================================================================

find_package(CGAL REQUIRED PATHS "${CMAKE_CURRENT_SOURCE_DIR}/../../cgal-build" NO_DEFAULT_PATH)
find_package(CUDAToolkit REQUIRED)
find_package(Threads REQUIRED)
find_package(TBB REQUIRED)

# GNU readline is only used for the interactive CLI on POSIX systems.
# On Windows the CLI falls back to std::getline (see apps/mainApplication.cpp).
if(NOT WIN32)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(READLINE REQUIRED readline)
endif()

# CUDA architectures to build for (override with -DMESH_CUDA_ARCHITECTURES=...)
set(MESH_CUDA_ARCHITECTURES "86" CACHE STRING "CUDA architectures to build for")

# Path variables for cuBQL and CUDA context
set(CUBQL_ROOT "/home/yury/Projects/CPPProj/Cubql" CACHE PATH "Path to cuBQL source root")
set(CUBQL_BUILD "${CUBQL_ROOT}/build" CACHE PATH "Path to a cuBQL build directory (only used for prebuilt shared libs on non-Windows)")

# Use the same (shared) CUDA runtime for MainApp and the cuBQL libraries on Windows.
if(WIN32)
    set(CMAKE_CUDA_RUNTIME_LIBRARY Shared)
endif()

# cuBQL's precompiled builder instantiations.
#  * Linux/macOS: link the prebuilt shared libraries from CUBQL_BUILD (as before).
#  * Windows (or -DCUBQL_USE_SUBPROJECT=ON): build the static instantiation targets
#    straight from the cuBQL sources; no manual cuBQL build step and no .so files needed.
option(CUBQL_USE_SUBPROJECT "Build cuBQL instantiation libraries from CUBQL_ROOT as part of this project" ${WIN32})
if(CUBQL_USE_SUBPROJECT)
    set(CMAKE_CUDA_ARCHITECTURES "${MESH_CUDA_ARCHITECTURES}")
    add_subdirectory("${CUBQL_ROOT}" "${CMAKE_CURRENT_BINARY_DIR}/cubql" EXCLUDE_FROM_ALL)
    set(CUBQL_LIBS cuBQL_cuda_float2_static cuBQL_cuda_float3_static)
else()
    set(CUBQL_LIBS
        "${CUBQL_BUILD}/libcuBQL_cuda_float2.so"
        "${CUBQL_BUILD}/libcuBQL_cuda_float3.so")
endif()

# Shared Link Libraries List
set(SHARED_LINK_LIBS
    CGAL::CGAL
    Threads::Threads
    TBB::tbb
    ${CUBQL_LIBS}
)
if(NOT WIN32)
    # On Windows gmp/mpfr come in through CGAL::CGAL (vcpkg import libs)
    list(APPEND SHARED_LINK_LIBS gmp mpfr)
endif()

# Fetch Polyscope
include(FetchContent)
FetchContent_Declare(
    polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(polyscope)

# ====================================================================
# HELPER FUNCTION
# ====================================================================

function(configure_geometry_target TARGET_NAME)
    set_target_properties(${TARGET_NAME} PROPERTIES
        CUDA_ARCHITECTURES "${MESH_CUDA_ARCHITECTURES}"
        CUDA_SEPARABLE_COMPILATION ON
    )

    target_compile_definitions(${TARGET_NAME} PRIVATE CGAL_LINKED_WITH_TBB)

    if(MSVC)
        # MSVC: /bigobj is needed for CGAL's template-heavy translation units,
        # /utf-8 because the sources contain UTF-8 text, M_PI & co. need _USE_MATH_DEFINES.
        target_compile_definitions(${TARGET_NAME} PRIVATE
            _USE_MATH_DEFINES NOMINMAX _CRT_SECURE_NO_WARNINGS)
        # /FI force_std_size_t.h (CUDA TUs only): works around error C2039 "'size_t' is not a
        # member of 'std'" in nvcc's generated .cudafe1.cpp (see the header for details).
        target_compile_options(${TARGET_NAME} PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:/bigobj /utf-8>
            $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math -Xcompiler=/Zc:preprocessor,/bigobj,/utf-8,/FI${CMAKE_CURRENT_SOURCE_DIR}/src/common/force_std_size_t.h>
        )
    else()
        set_target_properties(${TARGET_NAME} PROPERTIES
            BUILD_WITH_INSTALL_RPATH TRUE
            INSTALL_RPATH "${CUBQL_BUILD}"
        )
        target_compile_options(${TARGET_NAME} PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG>
            $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>
        )
    endif()

    # NOTE: ${CMAKE_CURRENT_SOURCE_DIR} must NOT be an include directory: it contains a file
    # called "VERSION", which on case-insensitive file systems (Windows) shadows the standard
    # <version> header and breaks the CUDA/STL headers.
    target_include_directories(${TARGET_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CUBQL_ROOT}
        ${CUBQL_ROOT}/samples/common
        ${CUDAToolkit_INCLUDE_DIRS}
    )

    target_link_libraries(${TARGET_NAME} PRIVATE
        CUDA::cudart
        ${SHARED_LINK_LIBS}
    )
endfunction()

# ====================================================================
# PRIMARY TARGET: MainApp
# ====================================================================

set_source_files_properties(
    apps/mainApplication.cpp
    src/CPU/TestSuite.cpp
    PROPERTIES LANGUAGE CUDA
)

add_executable(MainApp
    apps/mainApplication.cpp
    src/CPU/MeshIntersectionCurves.cpp
    src/CPU/TestSuite.cpp
    # Top-Level Entry Points
    src/GPUIntersector/KernelBVHController.cu
    src/GPUIntersector/StandaloneBVHPipeline.cu

    # Traversal & Cross-Checking
    src/traversal/crossCheck.cu
    src/traversal/rapidDescendKernel.cu
    src/traversal/DualTreeStep.cu
    src/traversal/prune_pipeline.cu

    # Single-Tree Traversal & Predicate Dispatcher
    src/single_tree/SingleTreeBatchCommon.cu
    src/single_tree/SingleTreeBatchIntersector.cu

    # Core Infrastructure & Shared Utilities
    src/common/utils.cu
    src/Warmup/cuda_warmup.cu

    # Exact Geometric Predicates
    src/GPU_exact_predicates/TwoSumDoublePredicate.cu
    src/GPU_exact_predicates/InexactBigIntFilter.cu
    src/GPU_exact_predicates/ExactBigIntPredicate.cu
    src/GPU_exact_predicates/OriginShiftedFloatPredicate.cu
    src/GPU_exact_predicates/DoublePredicate.cu
)

configure_geometry_target(MainApp)

# Target-specific additions for MainApp
if(NOT WIN32)
    target_include_directories(MainApp PRIVATE ${READLINE_INCLUDE_DIRS})
    target_link_libraries(MainApp PRIVATE ${READLINE_LIBRARIES})
endif()
target_link_libraries(MainApp PRIVATE polyscope)

# Windows: copy the runtime DLLs (TBB, GMP, MPFR, ...) next to the executable.
if(WIN32)
    add_custom_command(TARGET MainApp POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E $<IF:$<BOOL:$<TARGET_RUNTIME_DLLS:MainApp>>,copy_if_different,true>
                $<TARGET_RUNTIME_DLLS:MainApp> $<TARGET_FILE_DIR:MainApp>
        COMMAND_EXPAND_LISTS)
endif()

# ====================================================================
# ARCHIVED / LEGACY TARGETS (OPTIONAL BUILD)
# ====================================================================

if(BUILD_ARCHIVED_TARGETS)
    message(STATUS "Building archived/legacy targets enabled.")

    # TARGET: meshIntersectionRegistration
    add_executable(meshIntersectionRegistration apps/archive/main.cpp src/registration/kernels.cu)
    configure_geometry_target(meshIntersectionRegistration)

    # TARGET: meshIntersectionBatch
    add_executable(meshIntersectionBatch apps/archive/main_massive.cpp src/registration/kernels.cu)
    configure_geometry_target(meshIntersectionBatch)

    # TARGET: meshIntersectionCustom
    add_executable(meshIntersectionCustom apps/archive/mainV2.cpp src/GPU_exact_predicates/kernelsV2.cu src/GPU_exact_predicates/InexactBigIntFilter.cu)
    configure_geometry_target(meshIntersectionCustom)

    # TARGET: meshVoxelizer
    add_executable(meshVoxelizer apps/archive/mainVoxel.cpp src/voxelizer/kernelsVoxel.cu)
    configure_geometry_target(meshVoxelizer)

    # TARGET: testBVH
    add_executable(testBVH apps/archive/mainTestBVH.cpp src/archive/kernelsTestBVH.cu src/traversal/rapidDescendKernel.cu src/traversal/crossCheck.cu src/archive/volumeSanityCheck.cu src/archive/batchedCrossIntersection.cu src/archive/batchedCrossIntersectionBruteV2.cu src/traversal/DualTreeStep.cu src/traversal/prune_pipeline.cu src/GPU_exact_predicates/InexactBigIntFilter.cu)
    configure_geometry_target(testBVH)

    # TARGET: dualTreeExecution
    add_executable(dualTreeExecution apps/archive/mainDualTreeTestLauncher.cpp src/traversal/crossCheck.cu src/GPU_exact_predicates/InexactBigIntFilter.cu src/DualTreeTraversal/DualTreeLauncher.cu src/DualTreeTraversal/DualTreeMeshIntersection.cu)
    configure_geometry_target(dualTreeExecution)

    # TARGET: testBVHV2
    add_executable(testBVHV2 apps/archive/mainTestBVHV2.cpp src/archive/kernelsTestBVHAlternative.cu src/traversal/rapidDescendKernel.cu src/archive/batchedCrossIntersectionV2.cu src/archive/batchedCrossIntersectionV3.cu src/single_tree/SingleTreeBatchCommon.cu src/traversal/DualTreeStep.cu src/traversal/prune_pipeline.cu src/Warmup/cuda_warmup.cu src/GPU_exact_predicates/ExactBigIntPredicate.cu src/common/utils.cu)
    configure_geometry_target(testBVHV2)

endif()