// MU4IN901 - MODEL - Implementation project
//
// Authors:
// Andrei Aralov (andrei.aralov@etu.sorbonne-universite.fr)
// Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
// Ruchi  Thareja (ruchi.thareja@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Benchmarks of LU decomposition and LU inversion.

#include <benchmark/benchmark.h>
#include "Mat.hpp"
#include "LU.hpp" // Include LU and PLUQ functions
#include "random.hpp" // Include random matrix generation utilities

namespace {

using T = double;

// Benchmark for LU decomposition
void BM_LU(benchmark::State& state) {
    auto mat = M::generateRandomMatrix<T>(state.range( 0 ), state.range( 0 ));
    for (auto _ : state) {
        benchmark::DoNotOptimize(M::LU(mat));
    }
}

// Benchmark for PLUQ decomposition
void BM_PLUQ(benchmark::State& state) {
    auto mat = M::generateRandomMatrix<T>(state.range( 0 ), state.range( 0 ));
    for (auto _ : state) {
        benchmark::DoNotOptimize(M::PLUQ(mat));
    }
}

// Benchmark for LU inversion
void BM_LUInversion( benchmark::State& state )
{
    auto mat = M::generateRandomMatrix<T>( state.range( 0 ), state.range( 0 ) );
    for ( auto _ : state )
    {
        benchmark::DoNotOptimize(M::inverseLU( mat ));
    }
}

} // namespace

// Register benchmarks
BENCHMARK(BM_LU)->DenseRange( 100, 2600, 500 );
BENCHMARK(BM_PLUQ)->DenseRange( 100, 2600, 500 );
BENCHMARK(BM_LUInversion)->DenseRange( 100, 2600, 500 );

BENCHMARK_MAIN();
