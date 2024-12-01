#include <StrassenInversion.hpp>
#include <StrassenMultiplication.hpp>
#include <random.hpp>
#include <benchmark/benchmark.h>


using namespace M;

constexpr auto cTimeUnit = benchmark::kMillisecond;
using T = float;

template <Index min_size>
Mat<T> multiplyStrassen_TE( const Mat<T>& a, const Mat<T>& b )
{
    return multiplyStrassen<min_size>( a, b );
}


auto getStrassenMult( int min_size )
{
    switch ( min_size )
    {
        case 8:
            return multiplyStrassen_TE<8>;
        case 16:
            return multiplyStrassen_TE<16>;
        case 32:
            return multiplyStrassen_TE<32>;
        case 64:
            return multiplyStrassen_TE<64>;
        case 128:
            return multiplyStrassen_TE<128>;
        case 256:
            return multiplyStrassen_TE<256>;
        default:
            throw std::runtime_error( "This should never happen" );
    }
}


template <int input_from, int input_to, int input_step>
static void blockSizeAndInputSizeArgs( benchmark::internal::Benchmark* b )
{
    for ( int bs = 8; bs <= 256; bs = 2 * bs )
    {
        for ( int i = input_from; i <= input_to; i += input_step )
            b->Args( { bs, i } );
    }
}


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
    auto sz = state.range( 1 );
    auto mat = generateRandomMatrix<float>( sz, sz );
    auto mult = getStrassenMult( state.range( 0 ) );
    for ( auto _ : state )
    {
        benchmark::DoNotOptimize( inverseStrassen( mat, mult ) );
    }
}

BENCHMARK(BM_StrassenInversion)->Apply( blockSizeAndInputSizeArgs<100, 2600, 500> )->Unit( cTimeUnit );

static void BM_StrassenMultiplication( benchmark::State& state )
{
    auto n = state.range( 1 );
    auto a = generateRandomMatrix<float>(n, n);
    auto b = generateRandomMatrix<float>(n, n);
    auto mult = getStrassenMult( state.range( 0 ) );
    for ( auto _ : state )
    {
        benchmark::DoNotOptimize( mult( a, b ) );
    }
}

BENCHMARK(BM_StrassenMultiplication)->Apply( blockSizeAndInputSizeArgs<100, 2600, 500> )->Unit( cTimeUnit );