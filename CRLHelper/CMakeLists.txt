set(DEPS tbb tbbmalloc tbbmalloc_proxy igl::core igl_restricted::triangle imgui polyscope glfw glad CGAL cholmod)

file(GLOB_RECURSE HEADERS "*.h" "*.hpp")
file(GLOB_RECURSE SOURCES "*.cpp")

add_library(CRLHelper ${HEADERS} ${SOURCES})
target_link_libraries(CRLHelper ${DEPS})
