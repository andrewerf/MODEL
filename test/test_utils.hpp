// S1 2024, MODEL
//
// Author: Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Utility functions for comparing matrices in unit tests.

#include <Mat.hpp>

#include <gtest/gtest.h>

#include <random>

namespace M
{

// Helper function to compare two matrices up to a tolerance.
template<typename T, typename Impl1, typename Impl2, MatDim n, MatDim m>
bool isMatrixNear(const MatFacade<Impl1, T, n, m>& a,
                 const MatFacade<Impl2, T, n, m>& b,
                 T atol = 1e-8) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;

    for (Index i = 0; i < a.rows(); ++i)
        for (Index j = 0; j < a.cols(); ++j)
            if (std::abs(a(i, j) - b(i, j)) > atol)
                return false;
    return true;
}

// Compare two matrices (up to the given tolerance).
template<typename T>
void EXPECT_MATRIX_NEAR(const Mat<T>& a, const Mat<T>& b, T atol = 1e-8) {
    ASSERT_EQ( a.rows(), a.cols() );
    ASSERT_EQ( b.rows(), b.cols() );
    ASSERT_EQ( a.rows(), b.rows() );

    for (Index i = 0; i < a.rows(); ++i) {
        for (Index j = 0; j < a.cols(); ++j) {
            EXPECT_NEAR( a(i, j), b(i, j), atol ) << \
                "Matrices differ at position (" << i << "," << j << ")";
        }
    }
}

// A functor used by the GTest library to nicely print the matrix size
// instantiated in parametric tests.
struct ProductSizeToString {

template <typename T>
std::string operator()(const ::testing::TestParamInfo<T>& info) const {
    ::std::stringstream ss;
    auto [m, n, p] = info.param;
    ss << m << "x" << n << "_" << n << "x" << p;
    return ss.str();
}

};

}
