//
// Created by Andrey Aralov on 11/23/24.
//
#include <gtest/gtest.h>
#include <random>
#include <StrassenInversion.hpp>
#include <io.hpp>


using namespace M;


class MatFastInv : public ::testing::Test {
protected:
    const double epsilon = 1e-5;

    Mat<double> generateRandomMatrix( Index sz ) const
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis( -1.0, 1.0 );

        Mat<double> m( sz, sz );

        for ( Index i = 0; i < sz; ++i ) {
            for ( Index j = 0; j < sz; ++j ) {
                m( i, j ) = dis( gen );
            }
        }

        return m;
    }

    // Helper function to check if two matrices are equal within epsilon
    void expectMatricesNearlyEqual(const Mat<double>& a, const Mat<double>& b) const
    {
        ASSERT_EQ( a.rows(), a.cols() );
        ASSERT_EQ( b.rows(), b.cols() );
        ASSERT_EQ( a.rows(), b.rows() );

        for (Index i = 0; i < a.rows(); ++i) {
            for (Index j = 0; j < a.cols(); ++j) {
                EXPECT_NEAR( a(i, j), b(i, j), epsilon )
                                    << "Matrices differ at position (" << i << "," << j << ")";
            }
        }
    }
};


TEST_F( MatFastInv, NaiveIdentity )
{
    constexpr Index sz = 100;
    auto orig = generateRandomMatrix( sz );
    auto inv = inverseStrassen( orig );
    auto prod = orig * inv;
    expectMatricesNearlyEqual( prod, Mat<double>::Identity( sz, sz ) );
}

TEST_F( MatFastInv, NaiveDoubleInv )
{
    constexpr Index sz = 100;
    auto orig = generateRandomMatrix( sz );
    auto inv = inverseStrassen( orig );
    auto invinv = inverseStrassen( inv );
    expectMatricesNearlyEqual( orig, invinv );
}
