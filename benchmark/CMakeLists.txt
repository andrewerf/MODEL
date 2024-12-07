find_package(benchmark)

if ( NOT benchmark_FOUND )
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
    add_subdirectory(${CMAKE_SOURCE_DIR}/3rdparty/benchmark 3rdparty/benchmark)
    set_target_properties(benchmark PROPERTIES FOLDER extern)
endif ()

add_executable(mat_benchmark benchmark_strassen.cpp benchmark_lu.cpp)
target_compile_definitions(mat_benchmark PUBLIC "BLOCK_SIZES={ 8, 16, 32, 64, 128, 256 }")
target_link_libraries(mat_benchmark benchmark::benchmark_main)
target_compile_options(mat_benchmark PUBLIC -march=native -Ofast)
