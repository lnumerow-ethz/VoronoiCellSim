add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)
set(DEPS Half tbb tbbmalloc tbbmalloc_proxy mkl_intel_lp64 mkl_sequential mkl_core igl::core igl_restricted::triangle geogram imgui CGAL suitesparseconfig cholmod metis CRLHelper)

file(GLOB_RECURSE HEADERS "include/*.h" "include/*.hpp")
file(GLOB_RECURSE SOURCES "src/*.cpp")

add_executable(VoronoiFoam ${HEADERS} ${SOURCES})
target_link_libraries(VoronoiFoam ${DEPS})

configure_file("resource/cylinder.obj" "resource/cylinder.obj" COPYONLY)
configure_file("resource/ronan.jpg" "resource/ronan.jpg" COPYONLY)
configure_file("resource/ronan.txt" "resource/ronan.txt" COPYONLY)
