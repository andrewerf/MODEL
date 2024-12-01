//
// Created by Andrey Aralov on 11/23/24.
//

#include <random.hpp>
#include <StrassenInversion.hpp>

#include "test_utils.hpp"

using namespace M;

TEST( MatFastInv, NaiveIdentity )
{
    constexpr Index sz = 100;
    auto orig = generateRandomMatrix( sz, sz );
    auto inv = inverseStrassen( orig );
    auto prod = orig * inv;
    EXPECT_MATRIX_NEAR(prod, Mat<double>::Identity( sz, sz ), 1e-3);
}

TEST( MatFastInv, NaiveDoubleInv )
{
    constexpr Index sz = 100;
    auto orig = generateRandomMatrix( sz, sz );
    auto inv = inverseStrassen( orig );
    auto invinv = inverseStrassen( inv );
    EXPECT_MATRIX_NEAR(orig, invinv, 1e-3);
}

TEST( MatFastInv, BigMatrix )
{
    constexpr Index sz = 2000;
    auto orig = generateRandomMatrix( sz, sz );
    auto inv = inverseStrassen( orig );
}
