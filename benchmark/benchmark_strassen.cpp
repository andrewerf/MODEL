#include <StrassenInversion.hpp>
#include <StrassenMultiplication.hpp>
#include <random.hpp>
#include <benchmark/benchmark.h>


using namespace M;

constexpr auto cTimeUnit = benchmark::kMillisecond;


static void BM_NaiveStrassenInversion( benchmark::State& state )
{
    auto sz = state.range( 0 );
    auto mat = generateRandomMatrix<float>( sz, sz );
    for ( auto _ : state )
    {
        benchmark::DoNotOptimize( inverseStrassen( mat ) );
    }
}

BENCHMARK(BM_NaiveStrassenInversion)->DenseRange( 100, 2000, 200 )->Unit( cTimeUnit );

static void BM_StrassenMultiplication( benchmark::State& state )
{
    auto m = state.range( 0 );
    auto n = state.range( 1 );
    auto p = state.range( 2 );
    auto a = generateRandomMatrix(m, n);
    auto b = generateRandomMatrix(n, p);
    for ( auto _ : state )
    {
        benchmark::DoNotOptimize( multiplyStrassen( a, b ) );
    }
}

BENCHMARK(BM_StrassenMultiplication)
    ->Args(
        // 1.2s [BLAS: 68ms]
        { 2222, 2222, 2222 }
    )
    ->Args(
        // 3.9s [BLAS: 184ms]
        { 2222, 3333, 4444 }
    )
// The following test cases were used during development to validate the
// splitting of wide or tall matrices into squarish blocks
// in multiplySubmatrix and to set the constants kInnerSplitRatio
// and kOuterSplitRatio
    ->Args(
        // 1.1s (1.9s without split) [BLAS: 45ms]
        { 500, 4000, 4000 }
    )
    ->Args(
        // 1.1s (1.9s without split) [BLAS: 44ms]
        { 4000, 4000, 500 }
    )
    ->Args(
        // 1.3s [BLAS: 71ms]
        { 4000, 800, 4000 }
    )
    ->Args(
        // 0.8s (1.3s without split) [BLAS: 27ms]
        { 800, 8000, 800 }
    )->Unit( cTimeUnit );