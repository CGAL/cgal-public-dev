cmake_minimum_required(VERSION 3.18...3.31) 
project(meshIntersectionRegistration LANGUAGES CXX CUDA)

# ====================================================================
# GLOBAL OPTIONS & BUILD CONFIGURATION
# ====================================================================

# Toggle to build archived/legacy targets (Default is OFF to keep builds clean)
option(BUILD_ARCHIVED_TARGETS "Build archived legacy executables" OFF)

# Automatically generate compile_commands.json for VS Code / IDEs
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile commands for IDEs" FORCE)

# ====================================================================
# DEPENDENCIES & PACKAGE LOCATIONS
# ====================================================================

find_package(CGAL REQUIRED PATHS "${CMAKE_CURRENT_SOURCE_DIR}/../../cgal-build" NO_DEFAULT_PATH)
find_package(CUDAToolkit REQUIRED)
find_package(Threads REQUIRED)
find_package(TBB REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(READLINE REQUIRED readline)

# Path variables for cuBQL and CUDA context
set(CUBQL_ROOT "/home/yury/Projects/CPPProj/Cubql" CACHE PATH "Path to cuBQL source root")
set(CUBQL_BUILD "${CUBQL_ROOT}/build")

# Shared Link Libraries List
set(SHARED_LINK_LIBS
    CGAL::CGAL
    gmp
    mpfr
    Threads::Threads
    TBB::tbb
    "${CUBQL_BUILD}/libcuBQL_cuda_float2.so"
    "${CUBQL_BUILD}/libcuBQL_cuda_float3.so"
)

# Fetch Polyscope
include(FetchContent)
FetchContent_Declare(
    polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(polyscope)

# ====================================================================
# PRIMARY TARGET: MainApp
# ====================================================================

# Specify that mainApplication.cpp uses CUDA language processing
set_source_files_properties(apps/mainApplication.cpp PROPERTIES LANGUAGE CUDA)

add_executable(MainApp 
    apps/mainApplication.cpp 
    src/GPUIntersector/KenelBVHController.cu
    src/testBVH/crossCheck.cu
    src/testBVH/rapidDescendKernel.cu
    src/testBVH/batchedCrossIntersectionCommon.cu
    src/testBVH/batchedCrossIntersectionDouble.cu
    src/testBVH/DualTreeStep.cu
    src/Warmup/cuda_warmup.cu
    src/custom_pipeline/GPUPredicatesCheckV2.cu
    src/custom_pipeline/GPUPredicatesCheckDouble.cu
    src/custom_pipeline/GPUPredicatesCheckBigInteger.cu
    src/custom_pipeline/GPUPredicatesCheckBigIntegerV2.cu
    src/custom_pipeline/GPUPredicatesCheckDoubleAssisted.cu
    src/custom_pipeline/GPUPredicatesCheckShewchukFloat.cu
    src/custom_pipeline/GPUPredicatesTwoLap.cu
    src/testBVH/utils.cu
    src/testBVH/prune_pipeline.cu
    src/testBVH/kernelsTestBVHV3.cu 
)

set_target_properties(MainApp PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_ARCHITECTURES "120"
    CUDA_SEPARABLE_COMPILATION ON
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${CUBQL_BUILD}"
)

target_compile_definitions(MainApp PRIVATE CGAL_LINKED_WITH_TBB)

target_compile_options(MainApp PRIVATE 
    $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG>
    $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>
)

target_include_directories(MainApp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/GPU_Cubql_module          
    ${CMAKE_CURRENT_SOURCE_DIR}/GPU_Cubql_module/include  
    include
    ${CUBQL_ROOT}
    ${CUBQL_ROOT}/samples/common
    ${CUDAToolkit_INCLUDE_DIRS}
    ${READLINE_INCLUDE_DIRS}                               
)

target_link_libraries(MainApp PRIVATE 
    CUDA::cudart 
    ${SHARED_LINK_LIBS}
    ${READLINE_LIBRARIES}                                   
    polyscope                                               
)

# ====================================================================
# ARCHIVED / LEGACY TARGETS (OPTIONAL BUILD)
# ====================================================================
if(BUILD_ARCHIVED_TARGETS)
    message(STATUS "Building archived/legacy targets enabled.")

    # TARGET: meshIntersectionRegistration
    add_executable(meshIntersectionRegistration apps/archive/main.cpp src/registration/kernels.cu)
    set_target_properties(meshIntersectionRegistration PROPERTIES CXX_STANDARD 17 CUDA_STANDARD 17 CUDA_ARCHITECTURES "120" CUDA_SEPARABLE_COMPILATION ON BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH "${CUBQL_BUILD}")
    target_compile_definitions(meshIntersectionRegistration PRIVATE CGAL_LINKED_WITH_TBB)
    target_compile_options(meshIntersectionRegistration PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG> $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>)
    target_include_directories(meshIntersectionRegistration PRIVATE include ${CUBQL_ROOT} ${CUBQL_ROOT}/samples/common ${CUDAToolkit_INCLUDE_DIRS})
    target_link_libraries(meshIntersectionRegistration PRIVATE ${SHARED_LINK_LIBS})

    # TARGET: meshIntersectionBatch
    add_executable(meshIntersectionBatch apps/archive/main_massive.cpp src/registration/kernels.cu)
    set_target_properties(meshIntersectionBatch PROPERTIES CXX_STANDARD 17 CUDA_STANDARD 17 CUDA_ARCHITECTURES "120" CUDA_SEPARABLE_COMPILATION ON BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH "${CUBQL_BUILD}")
    target_compile_definitions(meshIntersectionBatch PRIVATE CGAL_LINKED_WITH_TBB)
    target_compile_options(meshIntersectionBatch PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG> $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>)
    target_include_directories(meshIntersectionBatch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} include ${CUBQL_ROOT} ${CUBQL_ROOT}/samples/common ${CUDAToolkit_INCLUDE_DIRS})
    target_link_libraries(meshIntersectionBatch PRIVATE ${SHARED_LINK_LIBS})

    # TARGET: meshIntersectionCustom
    add_executable(meshIntersectionCustom apps/archive/mainV2.cpp src/custom_pipeline/kernelsV2.cu src/custom_pipeline/GPUPredicatesCheck.cu)
    set_target_properties(meshIntersectionCustom PROPERTIES CXX_STANDARD 17 CUDA_STANDARD 17 CUDA_ARCHITECTURES "120" CUDA_SEPARABLE_COMPILATION ON BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH "${CUBQL_BUILD}")
    target_compile_definitions(meshIntersectionCustom PRIVATE CGAL_LINKED_WITH_TBB)
    target_compile_options(meshIntersectionCustom PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG> $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG>)
    target_include_directories(meshIntersectionCustom PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} include ${CUBQL_ROOT} ${CUBQL_ROOT}/samples/common)
    target_link_libraries(meshIntersectionCustom PRIVATE CUDA::cudart ${SHARED_LINK_LIBS})

    # TARGET: meshVoxelizer
    add_executable(meshVoxelizer apps/archive/mainVoxel.cpp src/voxelizer/kernelsVoxel.cu)
    set_target_properties(meshVoxelizer PROPERTIES CXX_STANDARD 17 CUDA_STANDARD 17 CUDA_ARCHITECTURES "120" CUDA_SEPARABLE_COMPILATION ON BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH "${CUBQL_BUILD}")
    target_compile_definitions(meshVoxelizer PRIVATE CGAL_LINKED_WITH_TBB)
    target_compile_options(meshVoxelizer PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG> $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>)
    target_include_directories(meshVoxelizer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} include ${CUBQL_ROOT} ${CUBQL_ROOT}/samples/common)
    target_link_libraries(meshVoxelizer PRIVATE CUDA::cudart ${SHARED_LINK_LIBS})

    # TARGET: testBVH
    add_executable(testBVH apps/archive/mainTestBVH.cpp src/testBVH/kernelsTestBVH.cu src/testBVH/rapidDescendKernel.cu src/testBVH/crossCheck.cu src/testBVH/volumeSanityCheck.cu src/testBVH/batchedCrossIntersection.cu src/testBVH/batchedCrossIntersectionBruteV2.cu src/testBVH/DualTreeStep.cu src/testBVH/prune_pipeline.cu src/custom_pipeline/GPUPredicatesCheck.cu)
    set_target_properties(testBVH PROPERTIES CXX_STANDARD 17 CUDA_STANDARD 17 CUDA_ARCHITECTURES "120" CUDA_SEPARABLE_COMPILATION ON BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH "${CUBQL_BUILD}")
    target_compile_definitions(testBVH PRIVATE CGAL_LINKED_WITH_TBB)
    target_compile_options(testBVH PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG> $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>)
    target_include_directories(testBVH PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/GPU_Cubql_module ${CMAKE_CURRENT_SOURCE_DIR}/GPU_Cubql_module/include include ${CUBQL_ROOT} ${CUBQL_ROOT}/samples/common ${CUDAToolkit_INCLUDE_DIRS})
    target_link_libraries(testBVH PRIVATE CUDA::cudart ${SHARED_LINK_LIBS})

    # TARGET: dualTreeExecution
    add_executable(dualTreeExecution apps/archive/mainDualTreeTestLauncher.cpp src/testBVH/crossCheck.cu src/custom_pipeline/GPUPredicatesCheck.cu src/DualTreeTraversal/DualTreeLauncher.cu src/DualTreeTraversal/DualTreeMeshIntersection.cu)
    set_target_properties(dualTreeExecution PROPERTIES CXX_STANDARD 17 CUDA_STANDARD 17 CUDA_ARCHITECTURES "120" CUDA_SEPARABLE_COMPILATION ON BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH "${CUBQL_BUILD}")
    target_compile_definitions(dualTreeExecution PRIVATE CGAL_LINKED_WITH_TBB)
    target_compile_options(dualTreeExecution PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG> $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>)
    target_include_directories(dualTreeExecution PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/GPU_Cubql_module ${CMAKE_CURRENT_SOURCE_DIR}/GPU_Cubql_module/include include ${CUBQL_ROOT} ${CUBQL_ROOT}/samples/common ${CUDAToolkit_INCLUDE_DIRS})
    target_link_libraries(dualTreeExecution PRIVATE CUDA::cudart ${SHARED_LINK_LIBS})

    # TARGET: testBVHV2
    add_executable(testBVHV2 apps/archive/mainTestBVHV2.cpp src/testBVH/kernelsTestBVHAlternative.cu src/testBVH/rapidDescendKernel.cu src/testBVH/batchedCrossIntersectionV2.cu src/testBVH/batchedCrossIntersectionV3.cu src/testBVH/batchedCrossIntersectionCommon.cu src/testBVH/DualTreeStep.cu src/testBVH/prune_pipeline.cu src/Warmup/cuda_warmup.cu src/custom_pipeline/GPUPredicatesCheckV2.cu src/testBVH/utils.cu)
    set_target_properties(testBVHV2 PROPERTIES CXX_STANDARD 17 CUDA_STANDARD 17 CUDA_ARCHITECTURES "120" CUDA_SEPARABLE_COMPILATION ON BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH "${CUBQL_BUILD}")
    target_compile_definitions(testBVHV2 PRIVATE CGAL_LINKED_WITH_TBB)
    target_compile_options(testBVHV2 PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG> $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>)
    target_include_directories(testBVHV2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/GPU_Cubql_module ${CMAKE_CURRENT_SOURCE_DIR}/GPU_Cubql_module/include include ${CUBQL_ROOT} ${CUBQL_ROOT}/samples/common ${CUDAToolkit_INCLUDE_DIRS})
    target_link_libraries(testBVHV2 PRIVATE CUDA::cudart ${SHARED_LINK_LIBS})

endif()