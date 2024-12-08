// MU4IN901 - MODEL - Implementation project
//
// Authors:
// Andrei Aralov (andrei.aralov@etu.sorbonne-universite.fr)
// Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
// Ruchi  Thareja (ruchi.thareja@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Unit tests for the Mat class.

#include <gtest/gtest.h>
#include <Mat.hpp>

using namespace M;

TEST(Mat, StaticCreationAndIndexing) {
    Mat<int, 3, 4> matrix;

    EXPECT_EQ(matrix.rows(), 3);
    EXPECT_EQ(matrix.cols(), 4);

    matrix(0, 0) = 1;
    matrix(2, 3) = 42;

    EXPECT_EQ(matrix(0, 0), 1);
    EXPECT_EQ(matrix(2, 3), 42);
}

TEST(Mat, DynamicCreationAndIndexing) {
    Mat<int, DynamicMatDim, DynamicMatDim> matrix(3, 4);

    EXPECT_EQ(matrix.rows(), 3);
    EXPECT_EQ(matrix.cols(), 4);

    matrix(0, 0) = 1;
    matrix(2, 3) = 42;

    EXPECT_EQ(matrix(0, 0), 1);
    EXPECT_EQ(matrix(2, 3), 42);
}

TEST(Mat, MixedCreationAndIndexing) {
    Mat<int, 3, DynamicMatDim> matrix( 3, 4 );

    EXPECT_EQ(matrix.rows(), 3);
    EXPECT_EQ(matrix.cols(), 4);

    matrix(0, 0) = 1;
    matrix(2, 3) = 42;

    EXPECT_EQ(matrix(0, 0), 1);
    EXPECT_EQ(matrix(2, 3), 42);
}

TEST(Mat, InitializerListCreation) {
    Mat<int, 3, 3> mat{
        { 1, 2, 3 },
        { 4, 5, 6 },
        { 7, 8, 9 }
    };

    int k = 0;
    for ( Index i = 0; i < 3; ++i )
        for ( Index j = 0; j < 3; ++j )
            EXPECT_EQ( mat( i, j ), (++k) );

//    Mat<int, 1, 2> mat2{
//        { 1, 2, 3 }
//    };
}

TEST(Mat, Addition) {
    using namespace M;
    using MatrixType = Mat<int, MatDim(3), MatDim(4)>;
    MatrixType matrix1, matrix2, matrix_sum;

    // Initialize the matrices
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            matrix1(i, j) = i * 4 + j + 1;
            matrix2(i, j) = 10 - (i * 4 + j + 1);
        }
    }

    matrix_sum = matrix1 + matrix2;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_EQ(matrix_sum(i, j), matrix1(i, j) + matrix2(i, j));
        }
    }

    // Test addition with a scalar
//    MatrixType matrix3 = matrix1 + 5;
//    for (int i = 0; i < 3; i++) {
//        for (int j = 0; j < 4; j++) {
//            EXPECT_EQ(matrix3(i, j), matrix1(i, j) + 5);
//        }
//    }
}

TEST(Mat, Multiplication) {
    using namespace M;
    using MatrixType1 = Mat<int, MatDim(3), MatDim(4)>;
    using MatrixType2 = Mat<int, MatDim(4), MatDim(2)>;
    using ResultType = Mat<int, MatDim(3), MatDim(2)>;

    MatrixType1 matrix1;
    MatrixType2 matrix2;
    ResultType matrix_product;

    // Initialize the matrices
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            matrix1(i, j) = i * 4 + j + 1;
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 2; j++) {
            matrix2(i, j) = 10 - (i * 2 + j + 1);
        }
    }

    matrix_product = matrix1 * matrix2;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            int expected_value = 0;
            for (int k = 0; k < 4; k++) {
                expected_value += matrix1(i, k) * matrix2(k, j);
            }
            EXPECT_EQ(matrix_product(i, j), expected_value);
        }
    }

    // Test multiplication with a scalar
//    ResultType matrix3 = matrix1 * 5;
//    for (int i = 0; i < 3; i++) {
//        for (int j = 0; j < 4; j++) {
//            EXPECT_EQ(matrix3(i, j), matrix1(i, j) * 5);
//        }
//    }
}


// Helper function to initialize a matrix with sequential values
template<typename MatType>
void fillSequential(MatType& matrix) {
    for (int i = 0; i < matrix.rows(); i++) {
        for (int j = 0; j < matrix.cols(); j++) {
            matrix(i, j) = i * matrix.cols() + j + 1;
        }
    }
}

// Helper function to initialize a matrix with reverse sequential values
template<typename MatType>
void fillReverseSequential(MatType& matrix) {
    for (int i = 0; i < matrix.rows(); i++) {
        for (int j = 0; j < matrix.cols(); j++) {
            matrix(i, j) = 10 - (i * matrix.cols() + j + 1);
        }
    }
}

TEST(Mat, StaticAddition) {
    using namespace M;
    using MatrixType = Mat<int, MatDim(3), MatDim(4)>;
    MatrixType matrix1, matrix2, matrix_sum;

    fillSequential(matrix1);
    fillReverseSequential(matrix2);

    matrix_sum = matrix1 + matrix2;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_EQ(matrix_sum(i, j), matrix1(i, j) + matrix2(i, j));
        }
    }

    // Test addition with a scalar
    MatrixType matrix3 = matrix1 + 5;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_EQ(matrix3(i, j), matrix1(i, j) + 5);
        }
    }
}

TEST(Mat, DynamicAddition) {
    using namespace M;
    using DynamicMatrix = Mat<int, DynamicMatDim, DynamicMatDim>;

    DynamicMatrix matrix1(3, 4);
    DynamicMatrix matrix2(3, 4);
    DynamicMatrix matrix_sum(3, 4);

    fillSequential(matrix1);
    fillReverseSequential(matrix2);

    matrix_sum = matrix1 + matrix2;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_EQ(matrix_sum(i, j), matrix1(i, j) + matrix2(i, j));
        }
    }

    // Test addition with a scalar
    DynamicMatrix matrix3 = matrix1 + 5;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_EQ(matrix3(i, j), matrix1(i, j) + 5);
        }
    }
}

TEST(Mat, MixedAddition) {
    using namespace M;
    using StaticMatrix = Mat<int, 3, 4>;
    using MixedMatrix = Mat<int, 3, DynamicMatDim>;
    using DynamicMatrix = Mat<int, DynamicMatDim, DynamicMatDim>;

    StaticMatrix static_matrix;
    MixedMatrix mixed_matrix( 3, 4 );
    DynamicMatrix dynamic_matrix( 3, 4 );

    fillSequential( static_matrix );
    fillReverseSequential( mixed_matrix );
    fillSequential( dynamic_matrix );

    // Test different combinations
    auto result1 = static_matrix + mixed_matrix;
    auto result2 = mixed_matrix + dynamic_matrix;
    auto result3 = static_matrix + dynamic_matrix;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_EQ(result1(i, j), static_matrix(i, j) + mixed_matrix(i, j));
            EXPECT_EQ(result2(i, j), mixed_matrix(i, j) + dynamic_matrix(i, j));
            EXPECT_EQ(result3(i, j), static_matrix(i, j) + dynamic_matrix(i, j));
        }
    }
}

TEST(Mat, StaticMultiplication) {
    using namespace M;
    using MatrixType1 = Mat<int, 3, 4>;
    using MatrixType2 = Mat<int, 4, 2>;
    using ResultType = Mat<int, 3, 2>;

    MatrixType1 matrix1;
    MatrixType2 matrix2;
    ResultType matrix_product;

    fillSequential(matrix1);
    fillReverseSequential(matrix2);

    matrix_product = matrix1 * matrix2;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            int expected_value = 0;
            for (int k = 0; k < 4; k++) {
                expected_value += matrix1(i, k) * matrix2(k, j);
            }
            EXPECT_EQ(matrix_product(i, j), expected_value);
        }
    }
}

TEST(Mat, DynamicMultiplication) {
    using namespace M;
    using DynamicMatrix = Mat<int, DynamicMatDim, DynamicMatDim>;

    DynamicMatrix matrix1(3, 4);
    DynamicMatrix matrix2(4, 2);
    DynamicMatrix matrix_product(3, 2);

    fillSequential(matrix1);
    fillReverseSequential(matrix2);

    matrix_product = matrix1 * matrix2;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            int expected_value = 0;
            for (int k = 0; k < 4; k++) {
                expected_value += matrix1(i, k) * matrix2(k, j);
            }
            EXPECT_EQ(matrix_product(i, j), expected_value);
        }
    }
}

TEST(Mat, MixedMultiplication) {
    using namespace M;
    using StaticMatrix1 = Mat<int, 3, 4>;
    using StaticMatrix2 = Mat<int, 4, 2>;
    using MixedMatrix1 = Mat<int, 3, DynamicMatDim>;
    using MixedMatrix2 = Mat<int, DynamicMatDim, 2>;
    using DynamicMatrix = Mat<int, DynamicMatDim, DynamicMatDim>;

    StaticMatrix1 static_matrix1;
    StaticMatrix2 static_matrix2;
    MixedMatrix1 mixed_matrix1( 3, 4 );
    MixedMatrix2 mixed_matrix2( 4, 2 );
    DynamicMatrix dynamic_matrix1( 3, 4 );
    DynamicMatrix dynamic_matrix2( 4, 2 );

    fillSequential(static_matrix1);
    fillReverseSequential(static_matrix2);
    fillSequential(mixed_matrix1);
    fillReverseSequential(mixed_matrix2);
    fillSequential(dynamic_matrix1);
    fillReverseSequential(dynamic_matrix2);

    // Test different combinations
    auto result1 = static_matrix1 * static_matrix2;
    auto result2 = mixed_matrix1 * mixed_matrix2;
    auto result3 = dynamic_matrix1 * dynamic_matrix2;
    auto result4 = static_matrix1 * mixed_matrix2;
    auto result5 = mixed_matrix1 * dynamic_matrix2;

    // Verify results for each combination
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            int expected_value = 0;
            for (int k = 0; k < 4; k++) {
                expected_value += static_matrix1(i, k) * static_matrix2(k, j);
            }
            EXPECT_EQ(result1(i, j), expected_value);

            expected_value = 0;
            for (int k = 0; k < 4; k++) {
                expected_value += mixed_matrix1(i, k) * mixed_matrix2(k, j);
            }
            EXPECT_EQ(result2(i, j), expected_value);

            expected_value = 0;
            for (int k = 0; k < 4; k++) {
                expected_value += dynamic_matrix1(i, k) * dynamic_matrix2(k, j);
            }
            EXPECT_EQ(result3(i, j), expected_value);

            expected_value = 0;
            for (int k = 0; k < 4; k++) {
                expected_value += static_matrix1(i, k) * mixed_matrix2(k, j);
            }
            EXPECT_EQ(result4(i, j), expected_value);

            expected_value = 0;
            for (int k = 0; k < 4; k++) {
                expected_value += mixed_matrix1(i, k) * dynamic_matrix2(k, j);
            }
            EXPECT_EQ(result5(i, j), expected_value);
        }
    }
}

TEST(Mat, ScalarOperations) {
    using namespace M;
    using StaticMatrix = Mat<int, 3, 4>;
    using MixedMatrix = Mat<int, 3, DynamicMatDim>;
    using DynamicMatrix = Mat<int, DynamicMatDim, DynamicMatDim>;

    StaticMatrix static_matrix;
    MixedMatrix mixed_matrix( 3, 4 );
    DynamicMatrix dynamic_matrix( 3, 4 );

    fillSequential(static_matrix);
    fillSequential(mixed_matrix);
    fillSequential(dynamic_matrix);

    const int scalar = 5;

    // Test multiplication with scalar for all matrix types
    auto static_result = static_matrix * scalar;
    auto mixed_result = mixed_matrix * scalar;
    auto dynamic_result = dynamic_matrix * scalar;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_EQ(static_result(i, j), static_matrix(i, j) * scalar);
            EXPECT_EQ(mixed_result(i, j), mixed_matrix(i, j) * scalar);
            EXPECT_EQ(dynamic_result(i, j), dynamic_matrix(i, j) * scalar);
        }
    }

    // Test scalar multiplication from both sides
    auto left_scalar = scalar * static_matrix;
    auto right_scalar = static_matrix * scalar;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_EQ(left_scalar(i, j), right_scalar(i, j));
        }
    }
}

TEST(Mat, CopyConstructor) {
    // Test static size matrix copy
    {
        Mat<double, MatDim(2), MatDim(2)> original;
        for(Index i = 0; i < 2; ++i)
            for(Index j = 0; j < 2; ++j)
                original(i, j) = i * 2 + j;

        Mat<double, MatDim(2), MatDim(2)> copied(original);

        // Verify values are correctly copied
        for(Index i = 0; i < 2; ++i)
            for(Index j = 0; j < 2; ++j)
                EXPECT_DOUBLE_EQ(copied(i, j), original(i, j));

        // Verify independence
        copied(0, 0) = 99.0;
        EXPECT_NE(copied(0, 0), original(0, 0));
    }

    // Test dynamic size matrix copy
    {
        Mat<double> original(3, 4);
        for(Index i = 0; i < 3; ++i)
            for(Index j = 0; j < 4; ++j)
                original(i, j) = i * 4 + j;

        Mat<double> copied(original);

        // Verify dimensions
        EXPECT_EQ(copied.rows(), original.rows());
        EXPECT_EQ(copied.cols(), original.cols());

        // Verify values
        for(Index i = 0; i < 3; ++i)
            for(Index j = 0; j < 4; ++j)
                EXPECT_DOUBLE_EQ(copied(i, j), original(i, j));
    }

    // Test copy between static and dynamic matrices
    {
        Mat<double, 2, 2> static_mat;
        for(Index i = 0; i < 2; ++i)
            for(Index j = 0; j < 2; ++j)
                static_mat(i, j) = i + j;

        Mat<double> dynamic_copy(static_mat);
        EXPECT_EQ(dynamic_copy.rows(), 2);
        EXPECT_EQ(dynamic_copy.cols(), 2);

        for(Index i = 0; i < 2; ++i)
            for(Index j = 0; j < 2; ++j)
                EXPECT_DOUBLE_EQ(dynamic_copy(i, j), static_mat(i, j));
    }
}

// Test Identity matrix creation
TEST(Mat, IdentityMatrix) {
    // Test static size identity matrix
    {
        auto identity = Mat<double, MatDim(3), MatDim(3)>::Identity(3, 3);

        // Check dimensions
        EXPECT_EQ(identity.rows(), 3);
        EXPECT_EQ(identity.cols(), 3);

        // Check diagonal elements are 1
        for(Index i = 0; i < 3; ++i)
            EXPECT_DOUBLE_EQ(identity(i, i), 1.0);

        // Check off-diagonal elements are 0
        for(Index i = 0; i < 3; ++i)
            for(Index j = 0; j < 3; ++j)
                if(i != j)
                    EXPECT_DOUBLE_EQ(identity(i, j), 0.0);
    }

    // Test identity matrix properties
    {
        auto I = Mat<double, MatDim(2), MatDim(2)>::Identity(2, 2);
        Mat<double, MatDim(2), MatDim(2)> A;

        // Fill A with some values
        A(0, 0) = 1; A(0, 1) = 2;
        A(1, 0) = 3; A(1, 1) = 4;

        // Test A * I = A
        auto AI = A * I;
        for(Index i = 0; i < 2; ++i)
            for(Index j = 0; j < 2; ++j)
                EXPECT_DOUBLE_EQ(AI(i, j), A(i, j));

        // Test I * A = A
        auto IA = I * A;
        for(Index i = 0; i < 2; ++i)
            for(Index j = 0; j < 2; ++j)
                EXPECT_DOUBLE_EQ(IA(i, j), A(i, j));
    }

    // Test with different types
    {
        auto identity_int = Mat<int, MatDim(2), MatDim(2)>::Identity(2, 2);
        EXPECT_EQ(identity_int(0, 0), 1);
        EXPECT_EQ(identity_int(1, 1), 1);
        EXPECT_EQ(identity_int(0, 1), 0);
        EXPECT_EQ(identity_int(1, 0), 0);
    }
}

// Test mixed operations with identity matrix
TEST(Mat, IdentityOperations) {
    const Index size = 3;
    auto I = Mat<double, MatDim(3), MatDim(3)>::Identity(size, size);
    Mat<double, MatDim(3), MatDim(3)> A;

    // Fill A with some values
    for(Index i = 0; i < size; ++i)
        for(Index j = 0; j < size; ++j)
            A(i, j) = i * size + j + 1;

    // Test A + I
    auto AI_sum = A + I;
    for(Index i = 0; i < size; ++i) {
        for(Index j = 0; j < size; ++j) {
            if(i == j)
                EXPECT_DOUBLE_EQ(AI_sum(i, j), A(i, j) + 1.0);
            else
                EXPECT_DOUBLE_EQ(AI_sum(i, j), A(i, j));
        }
    }

    // Test scalar multiplication
    auto I_scaled = 2.0 * I;
    for(Index i = 0; i < size; ++i) {
        for(Index j = 0; j < size; ++j) {
            if(i == j)
                EXPECT_DOUBLE_EQ(I_scaled(i, j), 2.0);
            else
                EXPECT_DOUBLE_EQ(I_scaled(i, j), 0.0);
        }
    }
}