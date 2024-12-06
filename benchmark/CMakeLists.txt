find_package(benchmark)

if ( NOT benchmark_FOUND )
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
    add_subdirectory(${CMAKE_SOURCE_DIR}/3rdparty/benchmark 3rdparty/benchmark)
    set_target_properties(benchmark PROPERTIES FOLDER extern)
endif ()

add_executable(mat_benchmark benchmark_strassen.cpp)
target_link_libraries(mat_benchmark benchmark::benchmark_main)
target_compile_options(mat_benchmark PUBLIC -march=native -Ofast)

add_executable(benchmark_lu benchmark_lu.cpp)
target_link_libraries(benchmark_lu benchmark::benchmark_main)
target_compile_options(benchmark_lu PUBLIC -march=native -Ofast)
