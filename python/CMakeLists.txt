add_custom_target(python)

add_dependencies(python NumericalMesh)


pybind11_add_module(NumericalMesh MODULE THIN_LTO OPT_SIZE numerical-mesh-py.cpp numerical-mesh-integration-py.cpp)
target_link_libraries(NumericalMesh PRIVATE numerical-mesh-lib)
