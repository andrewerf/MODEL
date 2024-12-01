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

BENCHMARK(BM_NaiveStrassenInversion)->DenseRange( 100, 2000, 500 )->Unit( cTimeUnit );


static void BM_StrassenInversion( benchmark::State& state )
{
    auto sz = state.range( 0 );
    auto mat = generateRandomMatrix<float>( sz, sz );
    auto mult = [] ( const auto& a, const auto& b )
    {
        return multiplyStrassen( a, b );
    };
    for ( auto _ : state )
    {
        benchmark::DoNotOptimize( inverseStrassen( mat, mult ) );
    }
}

BENCHMARK(BM_StrassenInversion)->DenseRange( 100, 2000, 500 )->Unit( cTimeUnit );

static void BM_StrassenMultiplication( benchmark::State& state )
{
    auto n = state.range( 0 );
    auto a = generateRandomMatrix(n, n);
    auto b = generateRandomMatrix(n, n);
    for ( auto _ : state )
    {
        benchmark::DoNotOptimize( multiplyStrassen( a, b ) );
    }
}

BENCHMARK(BM_StrassenMultiplication)->DenseRange( 100, 2500, 500 )->Unit( cTimeUnit );