add_subdirectory(integration)

pybind11_add_module(python MODULE THIN_LTO numerical-mesh-py.cpp)
target_link_libraries(python PRIVATE pybind11::headers numerical-mesh-lib)

set_target_properties(python PROPERTIES
                                INTERPROCEDURAL_OPTIMIZATION ON
                                CXX_VISIBILITY_PRESET hidden
                                VISIBLITY_INLINES_HIDDEN ON
				                LIBRARY_OUTPUT_NAME "NumericalMesh"
                                LIBRARY_OUTPUT_SUFFIX ""
				)

add_dependencies(python integration)
