// MU4IN901 - MODEL - Implementation project
//
// Authors:
// Andrei Aralov (andrei.aralov@etu.sorbonne-universite.fr)
// Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
// Ruchi  Thareja (ruchi.thareja@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Unit tests for the implementation of Strassen's inversion algorithm.

#include <random.hpp>
#include <StrassenInversion.hpp>
#include <StrassenMultiplication.hpp>

#include "test_utils.hpp"

using namespace M;

TEST( MatFastInv, NaiveIdentity )
{
    constexpr Index sz = 100;
    auto orig = generateRandomMatrix( sz, sz );
    auto inv = inverseStrassen( orig ).value();
    auto prod = orig * inv;
    EXPECT_MATRIX_NEAR(prod, Mat<double>::Identity( sz, sz ), 1e-3);
}

TEST( MatFastInv, NaiveDoubleInv )
{
    constexpr Index sz = 100;
    auto orig = generateRandomMatrix( sz, sz );
    auto inv = inverseStrassen( orig ).value();
    auto invinv = inverseStrassen( inv ).value();
    EXPECT_MATRIX_NEAR(orig, invinv, 1e-3);
}

constexpr auto strassenMultLambda = [] ( const auto& x, const auto& y )
{
    return multiplyStrassen( x, y );
};

TEST( MatFastInv, Identity )
{
    constexpr Index sz = 100;
    auto orig = generateRandomMatrix( sz, sz );
    auto inv = inverseStrassen( orig, strassenMultLambda ).value();
    auto prod = orig * inv;
    EXPECT_MATRIX_NEAR(prod, Mat<double>::Identity( sz, sz ), 1e-3);
}

TEST( MatFastInv, DoubleInv )
{
    constexpr Index sz = 100;
    auto orig = generateRandomMatrix( sz, sz );
    auto inv = inverseStrassen( orig, strassenMultLambda ).value();
    auto invinv = inverseStrassen( inv, strassenMultLambda ).value();
    EXPECT_MATRIX_NEAR(orig, invinv, 1e-3);
}


TEST( MatFastInv, BigMatrix )
{
    constexpr Index sz = 2000;
    auto orig = generateRandomMatrix( sz, sz );
    auto inv = inverseStrassen( orig ).value();
}

TEST( MatFastInv, NonInvertible )
{
    Mat<double> m(2, 2);
    m(0, 0) = 1;
    auto inv = inverseStrassen( m );
    EXPECT_FALSE(inv.has_value());
}
