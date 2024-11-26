#include <gtest/gtest.h>
#include <Mat.hpp>
#include <LU.hpp>

#include "test_utils.hpp"

using namespace M;

namespace
{

// Helper function to check if a matrix is lower triangular
template<typename T, MatDim n>
bool isLowerTriangular(const Mat<T, n, n>& mat, double epsilon = 1e-10) {
    for(Index i = 0; i < mat.rows(); ++i) {
        for(Index j = i + 1; j < mat.cols(); ++j) {
            if(std::abs(mat(i, j)) > epsilon)
                return false;
        }
    }
    return true;
}

// Helper function to check if a matrix is upper triangular
template<typename T, MatDim n>
bool isUpperTriangular(const Mat<T, n, n>& mat, double epsilon = 1e-10) {
    for(Index i = 1; i < mat.rows(); ++i) {
        for(Index j = 0; j < i; ++j) {
            if(std::abs(mat(i, j)) > epsilon)
                return false;
        }
    }
    return true;
}

// Helper function to check if a matrix is a permutation matrix
template<typename T, MatDim n>
bool isPermutation(const Mat<T, n, n>& mat) {
    // Check each row and column has exactly one 1 and rest 0s
    for(Index i = 0; i < mat.rows(); ++i) {
        T row_sum = 0;
        T col_sum = 0;
        bool found_row_one = false;
        bool found_col_one = false;

        for(Index j = 0; j < mat.cols(); ++j) {
            // Check row
            if(std::abs(mat(i, j)) > 0.5) {  // Use 0.5 to account for floating point
                if(std::abs(mat(i, j) - 1.0) > 1e-10 || found_row_one)
                    return false;
                found_row_one = true;
            }

            // Check column
            if(std::abs(mat(j, i)) > 0.5) {
                if(std::abs(mat(j, i) - 1.0) > 1e-10 || found_col_one)
                    return false;
                found_col_one = true;
            }
        }

        if(!found_row_one || !found_col_one)
            return false;
    }
    return true;
}

// Helper function to multiply matrices
template<typename T, MatDim n>
Mat<T, n, n> multiply4(const Mat<T, n, n>& P, const Mat<T, n, n>& L,
                       const Mat<T, n, n>& U, const Mat<T, n, n>& Q) {
    return P * L * U * Q;
}

}

// Test basic properties of LU decomposition
TEST(LU, LUBasicProperties) {
    // 2x2 matrix test
    {
        Mat<double, 2, 2> A;
        A(0, 0) = 4; A(0, 1) = 3;
        A(1, 0) = 6; A(1, 1) = 3;

        auto result = LU(A);

        // Test L properties
        EXPECT_DOUBLE_EQ(result.L(0, 0), 1.0);  // Diagonal should be 1
        EXPECT_DOUBLE_EQ(result.L(1, 1), 1.0);
        EXPECT_DOUBLE_EQ(result.L(0, 1), 0.0);  // Upper triangle should be 0

        // Test that A = L * U
        auto product = result.L * result.U;
        EXPECT_TRUE(matrix_near(product, A));
    }

    // 3x3 matrix test
    {
        Mat<double, 3, 3> A;
        // Using a well-conditioned matrix
        A(0, 0) = 4; A(0, 1) = -2; A(0, 2) = 1;
        A(1, 0) = -2; A(1, 1) = 4; A(1, 2) = -2;
        A(2, 0) = 1; A(2, 1) = -2; A(2, 2) = 4;

        auto result = LU(A);

        // Test L properties
        for (Index i = 0; i < 3; ++i) {
            EXPECT_DOUBLE_EQ(result.L(i, i), 1.0);  // Diagonal = 1
            for (Index j = i + 1; j < 3; ++j)
                EXPECT_DOUBLE_EQ(result.L(i, j), 0.0);  // Upper triangle = 0
        }

        // Test U properties
        for (Index i = 0; i < 3; ++i)
            for (Index j = 0; j < i; ++j)
                EXPECT_DOUBLE_EQ(result.U(i, j), 0.0);  // Lower triangle = 0

        // Test A = L * U
        auto product = result.L * result.U;
        EXPECT_TRUE(matrix_near(product, A));
    }
}

// Test special cases
TEST(LU, LUSpecialCases) {
    // Identity matrix
    {
        auto I = Mat<double, 3, 3>::Identity(3, 3);
        auto result = LU(I);

        // L should be identity
        EXPECT_TRUE(matrix_near(result.L, I));
        // U should be identity
        EXPECT_TRUE(matrix_near(result.U, I));
    }

    // Upper triangular matrix
    {
        Mat<double, 3, 3> A;
        for (Index i = 0; i < 3; ++i)
            for (Index j = 0; j < 3; ++j)
                A(i, j) = (j >= i) ? i + j + 1 : 0;

        auto result = LU(A);

        // L should be identity
        auto I = Mat<double, 3, 3>::Identity(3, 3);
        EXPECT_TRUE(matrix_near(result.L, I));
        // U should be the original matrix
        EXPECT_TRUE(matrix_near(result.U, A));
    }

    // Lower triangular matrix
    {
        Mat<double, 3, 3> A;
        for (Index i = 0; i < 3; ++i)
            for (Index j = 0; j < 3; ++j)
                A(i, j) = (j <= i) ? i + j + 1 : 0;

        auto result = LU(A);

        // Test A = L * U
        auto product = result.L * result.U;
        EXPECT_TRUE(matrix_near(product, A));
    }
}

// Test numerical stability
TEST(LU, LUNumericalStability) {
    // Test with a poorly conditioned matrix
    Mat<double, 3, 3> A;
    A(0, 0) = 1e-10; A(0, 1) = 1; A(0, 2) = 1;
    A(1, 0) = 1; A(1, 1) = 1; A(1, 2) = 1;
    A(2, 0) = 1; A(2, 1) = 1; A(2, 2) = 1;

    auto result = LU(A);

    // Test A = L * U with larger epsilon due to numerical instability
    auto product = result.L * result.U;
    EXPECT_TRUE(matrix_near(product, A, 1e-5));
}


// Test with different numeric types
TEST(LU, LUDifferentTypes) {
    // Test with float
    {
        Mat<float, 2, 2> A;
        A(0, 0) = 4.0f; A(0, 1) = 3.0f;
        A(1, 0) = 6.0f; A(1, 1) = 3.0f;

        auto result = LU(A);
        auto product = result.L * result.U;
        EXPECT_TRUE(matrix_near(product, A, 1e-6f));
    }

    // Test with double
    {
        Mat<double, 2, 2> A;
        A(0, 0) = 4.0; A(0, 1) = 3.0;
        A(1, 0) = 6.0; A(1, 1) = 3.0;

        auto result = LU(A);
        auto product = result.L * result.U;
        EXPECT_TRUE(matrix_near(product, A, 1e-10));
    }
}

TEST(PLUQ, SimpleMatrix2x2) {
    Mat<double, 2, 2> A;
    A(0,0) = 4; A(0,1) = 3;
    A(1,0) = 6; A(1,1) = 3;

    auto [P, L, U, Q] = PLUQ(A);

    // Test P and Q are permutation matrices
    EXPECT_TRUE(isPermutation(P));
    EXPECT_TRUE(isPermutation(Q));

    // Test L is lower triangular with 1s on diagonal
    EXPECT_TRUE(isLowerTriangular(L));
    for(Index i = 0; i < 2; ++i)
        EXPECT_NEAR(L(i,i), 1.0, 1e-10);

    // Test U is upper triangular
    EXPECT_TRUE(isUpperTriangular(U));

    // Test decomposition: A = P*L*U*Q
    auto reconstructed = multiply4(P, L, U, Q);
    for(Index i = 0; i < 2; ++i)
        for(Index j = 0; j < 2; ++j)
            EXPECT_NEAR(reconstructed(i,j), A(i,j), 1e-10);
}

TEST(PLUQ, Matrix3x3) {
    Mat<double, 3, 3> A;
    // Fill with a known invertible matrix
    A(0,0) = 2; A(0,1) = -1; A(0,2) = 0;
    A(1,0) = -1; A(1,1) = 2; A(1,2) = -1;
    A(2,0) = 0; A(2,1) = -1; A(2,2) = 2;

    auto [P, L, U, Q] = PLUQ(A);

    // Test properties of the decomposition
    EXPECT_TRUE(isPermutation(P));
    EXPECT_TRUE(isPermutation(Q));
    EXPECT_TRUE(isLowerTriangular(L));
    EXPECT_TRUE(isUpperTriangular(U));

    // Test reconstruction
    auto reconstructed = multiply4(P, L, U, Q);
    for(Index i = 0; i < 3; ++i)
        for(Index j = 0; j < 3; ++j)
            EXPECT_NEAR(reconstructed(i,j), A(i,j), 1e-10);
}

TEST(PLUQ, DiagonalMatrix) {
    Mat<double, 3, 3> A;
    // Create a diagonal matrix
    for(Index i = 0; i < 3; ++i)
        for(Index j = 0; j < 3; ++j)
            A(i,j) = (i == j) ? i + 1 : 0;

    auto [P, L, U, Q] = PLUQ(A);

    // For a diagonal matrix, U should be similar to original
    // and L should be mostly identity
    EXPECT_TRUE(isLowerTriangular(L));
    EXPECT_TRUE(isUpperTriangular(U));

    auto reconstructed = multiply4(P, L, U, Q);
    for(Index i = 0; i < 3; ++i)
        for(Index j = 0; j < 3; ++j)
            EXPECT_NEAR(reconstructed(i,j), A(i,j), 1e-10);
}

TEST(PLUQ, IdentityMatrix) {
    const Index size = 3;
    auto A = Mat<double>::Identity(size, size);

    auto [P, L, U, Q] = PLUQ(A);

    // For identity matrix, decomposition should be trivial
    EXPECT_TRUE(isPermutation(P));
    EXPECT_TRUE(isPermutation(Q));
    EXPECT_TRUE(isLowerTriangular(L));
    EXPECT_TRUE(isUpperTriangular(U));

    // Test reconstruction
    auto reconstructed = multiply4(P, L, U, Q);
    for(Index i = 0; i < size; ++i)
        for(Index j = 0; j < size; ++j)
            EXPECT_NEAR(reconstructed(i,j), A(i,j), 1e-10);
}

TEST(PLUQ, NearSingularMatrix) {
    Mat<double, 2, 2> A;
    // Create a nearly singular matrix
    A(0,0) = 1e-15; A(0,1) = 1;
    A(1,0) = 1;     A(1,1) = 1;

    auto [P, L, U, Q] = PLUQ(A);

    // The decomposition should still maintain basic properties
    EXPECT_TRUE(isPermutation(P));
    EXPECT_TRUE(isPermutation(Q));
    EXPECT_TRUE(isLowerTriangular(L));
    EXPECT_TRUE(isUpperTriangular(U));

    // Test reconstruction (with higher tolerance due to near-singularity)
    auto reconstructed = multiply4(P, L, U, Q);
    for(Index i = 0; i < 2; ++i)
        for(Index j = 0; j < 2; ++j)
            EXPECT_NEAR(reconstructed(i,j), A(i,j), 1e-10);
}