set(TEST_FILES
	main_test.cpp
    test_mesh.cpp
    test_mesh_integration.cpp
)

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(numerical-mesh-test ${TEST_FILES})
target_link_libraries(numerical-mesh-test ${GTEST_LIBRARIES} pthread numerical-mesh-lib)
GTEST_ADD_TESTS(numerical-mesh-test "" ${TEST_FILES})
message(STATUS "CXX flags : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
