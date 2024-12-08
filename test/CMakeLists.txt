find_package(GTest)
if( NOT GTest_FOUND )
    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
    set(BUILD_GTEST ON CACHE BOOL "" FORCE)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    add_subdirectory(${CMAKE_SOURCE_DIR}/3rdparty/googletest 3rdparty/googletest)
    set_target_properties(gtest PROPERTIES FOLDER extern)
    set_target_properties(gtest_main PROPERTIES FOLDER extern)
endif ()

add_executable(test_all test_mat.cpp test_lu.cpp test_mat_view.cpp test_strassen_inversion.cpp test_strassen_multiplication.cpp)
target_link_libraries(test_all GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(test_all)

# To compile the code making comparisons with BLAS (MacOS only for now)
if ( ${APPLE} )
    target_link_libraries(test_all "-framework Accelerate")
endif ()

set (CMAKE_CXX_FLAGS_RELEASE "-Ofast")
