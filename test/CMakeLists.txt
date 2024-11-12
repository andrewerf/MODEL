enable_testing()

add_executable(test_all test_mat.cpp test_lu.cpp)
target_link_libraries(test_all GTest::gtest_main)

find_package(GTest REQUIRED)
gtest_discover_tests(test_all)
