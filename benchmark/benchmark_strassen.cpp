// MU4IN901 - MODEL - Implementation project
//
// Authors:
// Andrei Aralov (andrei.aralov@etu.sorbonne-universite.fr)
// Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
// Ruchi  Thareja (ruchi.thareja@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Benchmark of matrix multiplication algorithms.

#include <StrassenInversion.hpp>
#include <StrassenMultiplication.hpp>
#include <random.hpp>
#include <benchmark/benchmark.h>
#include <stdexcept>

using namespace M;

constexpr auto cTimeUnit = benchmark::kMillisecond;
using T = double;

// Returns a version of multiplyStrassen specialized for a given value of
// min_size (minimum size of a leaf product)
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

// Generates the combination of parameters used to benchmark multiplication.
template <int input_from, int input_to, int input_step>
static void multBlockSizeAndInputSizeArgs( benchmark::internal::Benchmark* b )
{
    for ( int bs : MULT_BLOCK_SIZES )
    {
        for ( int i = input_from; i <= input_to; i += input_step )
            b->Args( { bs, i } );
    }
}

// Generates the combination of parameters used to benchmark inversion.
template <int input_from, int input_to, int input_step>
static void invBlockSizeAndInputSizeArgs( benchmark::internal::Benchmark* b )
{
    for ( int inv_bs : INV_BLOCK_SIZES )
    for ( int mult_bs : MULT_BLOCK_SIZES )
    {
        for ( int i = input_from; i <= input_to; i += input_step )
            b->Args( { inv_bs, mult_bs, i } );
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
    auto inv_bs = state.range( 0 );
    auto sz = state.range( 2 );
    auto mat = generateRandomMatrix<float>( sz, sz );
    auto mult = getStrassenMult( state.range( 1 ) );
    for ( auto _ : state )
    {
        benchmark::DoNotOptimize( inverseStrassen( mat, mult ) );
    }
}

BENCHMARK(BM_StrassenInversion)->Apply( invBlockSizeAndInputSizeArgs<100, 2600, 500> )->Unit( cTimeUnit );

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

BENCHMARK(BM_StrassenMultiplication)->Apply( multBlockSizeAndInputSizeArgs<100, 2600, 500> )->Unit( cTimeUnit );