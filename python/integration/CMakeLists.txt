pybind11_add_module(integration MODULE THIN_LTO numerical-mesh-integration-py.cpp)
target_link_libraries(integration PRIVATE pybind11::headers numerical-mesh-lib)

set_target_properties(integration PROPERTIES
                                INTERPROCEDURAL_OPTIMIZATION ON
                                CXX_VISIBILITY_PRESET hidden
                                VISIBLITY_INLINES_HIDDEN ON
				)
