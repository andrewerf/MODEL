// MU4IN901 - MODEL - Implementation project
//
// Authors:
// Andrei Aralov (andrei.aralov@etu.sorbonne-universite.fr)
// Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
// Ruchi  Thareja (ruchi.thareja@etu.sorbonne-universite.fr)
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
template<typename M1, typename M2, typename T = typename M1::ElemT>
void EXPECT_MATRIX_NEAR(const M1& a, const M2& b, T atol = 1e-8) {
    ASSERT_EQ( a.rows(), b.rows() );
    ASSERT_EQ( a.cols(), b.cols() );

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
