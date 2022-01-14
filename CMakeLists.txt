cmake_minimum_required(VERSION 3.14)

project(
	Numerical-mesh
	VERSION 0.1
	DESCRIPTION "C++ template library for numerical meshes"
	LANGUAGES CXX
)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})


message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} -------------")

add_library(numerical-mesh-lib INTERFACE)

target_include_directories(numerical-mesh-lib INTERFACE 
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
	)

if(BUILD_PYTHON)
	find_package(pybind11)
	if(pybind11_FOUND)
		add_subdirectory(python)
	else()
		message(STATUS "pybind11 not found, not building python bindings.")
	endif()
endif()

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	set(CMAKE_CXX_EXTENSIONS OFF)
	set_property(GLOBAL PROPERTY USE_FOLDERS ON)

	find_package(Doxygen)
	if(Doxygen_FOUND)
		add_subdirectory(docs)
	else()
		message(STATUS "Doxygen not found, not building docs")
	endif()

	set(WFLAGS -Werror -Wall -Wextra -pedantic -Wshadow -Wnon-virtual-dtor
		-Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic
		-Wconversion -Wsign-conversion -Wmisleading-indentation
		-Wnull-dereference -Wdouble-promotion -Wformat=2 -Weffc++)
	if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		set(WFLAGS ${WFLAGS} -Wduplicated-cond -Wduplicated-branches
			-Wlogical-op -Wuseless-cast)
	endif()
endif() # Main project only part ends here

# IDEs should put the headers in a nice place
source_group(
	TREE "${PROJECT_SOURCE_DIR}/include"
	PREFIX "Header Files"
	FILES ${HEADER_LIST}
)

if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING)
AND BUILD_TESTING)
	add_subdirectory(test)
endif()

install(TARGETS numerical-mesh-lib 
	EXPORT numerical-mesh-targets
	INCLUDES DESTINATION include/numerical-mesh)
install(EXPORT numerical-mesh-targets
	NAMESPACE Numerical-mesh::
	DESTINATION share/numerical-mesh
	)
install(DIRECTORY include/ DESTINATION include/numerical-mesh 
	FILES_MATCHING PATTERN "*.h")

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/numerical-mesh/numerical-mesh-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)
export(EXPORT numerical-mesh-targets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/numerical-mesh/numerical-mesh-targets.cmake"
  NAMESPACE Numerical-mesh::
)
configure_file(cmake/numerical-mesh-config.cmake
  "${CMAKE_CURRENT_BINARY_DIR}/numerical-mesh/numerical-mesh-config.cmake"
  COPYONLY
)

set(ConfigPackageLocation share/cmake/numerical-mesh)
install(EXPORT numerical-mesh-targets
  FILE
    numerical-mesh-targets.cmake
  NAMESPACE
  Numerical-mesh::
  DESTINATION
    ${ConfigPackageLocation}
)
install(
  FILES
    cmake/numerical-mesh-config.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/numerical-mesh/numerical-mesh-config-version.cmake"
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Devel
)
