//S1 2024
// Author: Ruchi Thareja M1 IQ
//-----------------------------------------------------------------

#include <benchmark/benchmark.h>
#include "Mat.hpp"
#include "LU.hpp" // Include LU and PLUQ functions
#include "random.hpp" // Include random matrix generation utilities

namespace {

// Benchmark for LU decomposition
template <typename T, M::MatDim n>
void BM_LU(benchmark::State& state) {
    auto mat = M::generateRandomMatrix<T>(n.val, n.val);
    for (auto _ : state) {
        benchmark::DoNotOptimize(M::LU(mat));
    }
}

// Benchmark for PLUQ decomposition
template <typename T, M::MatDim n>
void BM_PLUQ(benchmark::State& state) {
    auto mat = M::generateRandomMatrix<T>(n.val, n.val);
    for (auto _ : state) {
        benchmark::DoNotOptimize(M::PLUQ(mat));
    }
}

} // namespace

// Register benchmarks
BENCHMARK_TEMPLATE(BM_LU, double, M::MatDim(10));
BENCHMARK_TEMPLATE(BM_LU, double, M::MatDim(50));
BENCHMARK_TEMPLATE(BM_LU, double, M::MatDim(100));

BENCHMARK_TEMPLATE(BM_PLUQ, double, M::MatDim(10));
BENCHMARK_TEMPLATE(BM_PLUQ, double, M::MatDim(50));
BENCHMARK_TEMPLATE(BM_PLUQ, double, M::MatDim(100));

BENCHMARK_MAIN();
