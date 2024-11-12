enable_testing()

add_executable(test_mat test_mat.cpp)
target_link_libraries(test_mat GTest::gtest_main)

find_package(GTest REQUIRED)
gtest_discover_tests(test_mat)
