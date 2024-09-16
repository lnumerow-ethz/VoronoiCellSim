cmake_minimum_required(VERSION 3.2)

if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt.")
endif ()

project(lnumerowCRL C CXX)

message(STATUS "${CMAKE_BUILD_TYPE} Build")

set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)
set(VERSION_PATCH 0)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake)

set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -O0 -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O0 -fopenmp -Wno-unused-variable -Wno-unused-but-set-variable -Wsign-compare -Wreorder")

find_package(TBB REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(SuiteSparse REQUIRED)
find_package(GEOGRAM REQUIRED)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

find_package(CGAL REQUIRED COMPONENTS Core)

if (CUDA_FOUND)
    add_definitions(-DENABLE_AMGCL_CUDA)
endif ()

set(MKLROOT "/opt/intel/oneapi/mkl/2024.2")
add_definitions(-DEIGEN_USE_MKL_ALL)
add_definitions(-DMKL_LP64)
add_definitions(-DUSE_CHOLMOD)

include_directories(${GEOGRAM_INCLUDE_DIR})
include_directories(${CHOLMOD_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${MKLROOT}/include)
link_directories(${MKLROOT}/lib)
add_subdirectory(CRLHelper)
add_subdirectory(ThirdParty)
add_subdirectory(Projects)

# libigl
option(LIBIGL_WITH_OPENGL "Use OpenGL" ON)
option(LIBIGL_WITH_TRIANGLE "Use Triangle" ON)
option(LIBIGL_WITH_TETGEN "Use Tetgen" On)
# option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW" ON)
# option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui" ON)
option(LIBIGL_WITH_PNG "Use PNG" ON)
option(LIBIGL_WITH_XML "Use XML" ON)
option(LIBIGL_WITH_MOSEK "Use MOSEK" OFF)

# Libigl
include(libigl)
# Enable the target igl::glfw
# igl_include(glfw)
# igl_include(imgui)
igl_include(predicates)
igl_include(restricted triangle)

message("**************************************************")
message("C++ Cmake Flags: ${CMAKE_CXX_FLAGS}")
message("**************************************************")
