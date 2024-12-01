find_package(benchmark)

if ( NOT benchmark_FOUND )
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
    add_subdirectory(${CMAKE_SOURCE_DIR}/3rdparty/benchmark 3rdparty/benchmark)
    set_target_properties(benchmark PROPERTIES FOLDER extern)
endif ()

add_executable(benchmark benchmark_strassen.cpp)
target_link_libraries(benchmark benchmark::benchmark_main)