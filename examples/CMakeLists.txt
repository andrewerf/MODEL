project(MatExamples)

add_executable(pluq_simple pluq_simple.cpp)

add_executable(pluq_mpfr pluq_mpfr.cpp)
target_link_libraries(pluq_mpfr Boost::boost mpfr)
