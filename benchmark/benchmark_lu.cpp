//S1 2024
// Author: Ruchi Thareja M1 IQ
//-----------------------------------------------------------------

#include <benchmark/benchmark.h>
#include "Mat.hpp"
#include "LU.hpp" // Include LU and PLUQ functions
#include "random.hpp" // Include random matrix generation utilities

namespace {

using T = float;

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

} // namespace

// Register benchmarks
BENCHMARK(BM_LU)->DenseRange( 100, 2600, 500 );
BENCHMARK(BM_PLUQ)->DenseRange( 100, 2600, 500 );

BENCHMARK_MAIN();
