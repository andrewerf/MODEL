//
// Created by Andrey Aralov on 11/25/24.
//
#include <StrassenInversion.hpp>
#include <random.hpp>
#include <benchmark/benchmark.h>


using namespace M;

constexpr auto cTimeUnit = benchmark::kMillisecond;


static void BM_NaiveStrassenInversionBenchmark( benchmark::State& state )
{
    auto sz = state.range( 0 );
    auto mat = generateRandomMatrix<float>( sz, sz );
    for ( auto _ : state )
    {
        benchmark::DoNotOptimize( inverseStrassen( mat ) );
    }
}

BENCHMARK(BM_NaiveStrassenInversionBenchmark)->Arg( 100 )->Arg( 500 )->Arg( 1000 )->Unit( cTimeUnit );
