// MU4IN901 - MODEL - Implementation project
//
// Authors:
// Andrei Aralov (andrei.aralov@etu.sorbonne-universite.fr)
// Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
// Ruchi  Thareja (ruchi.thareja@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Unit tests for the MatView class.

#include <gtest/gtest.h>
#include <MatView.hpp>

using namespace M;

TEST(MatView, ConstructionAndAccess) {
    // Test construction from a matrix
    Mat<int, 3, 3> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    MatView view( matrix );

    // Test accessing elements
    EXPECT_EQ( view(0, 0), 1 );
    EXPECT_EQ( view(1, 1), 5 );
    EXPECT_EQ( view(2, 2), 9 );

    // Test no copying (compare addresses of the element)
    EXPECT_EQ( &view(0, 0), &matrix(0, 0) );
    EXPECT_EQ( &view(2, 1), &matrix(2, 1) );
    EXPECT_EQ( &view(1, 2), &matrix(1, 2) );
}

TEST(MatView, ConstructionFromSubMatrix) {
    // Test construction from a submatrix
    Mat<int, 4, 4> matrix = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
    auto view = matrix.submatrix<2, 2>( 1, 1 );

    // Test accessing elements of the submatrix
    EXPECT_EQ(view(0, 0), 6);
    EXPECT_EQ(view(0, 1), 7);
    EXPECT_EQ(view(1, 0), 10);
    EXPECT_EQ(view(1, 1), 11);

    auto view2 = matrix.submatrix( 1, 1, 2, 2 );
    EXPECT_EQ(view2(0, 0), 6);
    EXPECT_EQ(view2(0, 1), 7);
    EXPECT_EQ(view2(1, 0), 10);
    EXPECT_EQ(view2(1, 1), 11);

    // Test no copying (compare addresses of the element)
    for ( Index i = 0; i < 2; ++i )
        for ( Index j = 0; j < 2; ++j )
            EXPECT_EQ( &view( i, j ), &view2( i, j ) );
}


TEST(MatView, ConstructionFromMatView) {
    Mat<int, 4, 4> matrix = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
    auto view = matrix.submatrix<3, 3>( 1, 1 );
    auto view2 = view.submatrix<2, 1>( 0, 1 );

    EXPECT_EQ( &view2( 0, 0 ), &matrix( 1, 2 ) );
    EXPECT_EQ( &view2( 1, 0 ), &matrix( 2, 2 ) );
}

TEST(MatView, ConstructionFromMatViewConst) {
    const Mat<int, 4, 4> matrix = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
    const auto view = matrix.submatrix<3, 3>( 1, 1 );
    const auto view2 = view.submatrix<2, 1>( 0, 1 );

    EXPECT_EQ( &view2( 0, 0 ), &matrix( 1, 2 ) );
    EXPECT_EQ( &view2( 1, 0 ), &matrix( 2, 2 ) );
}

TEST(MatView, ConstructionFromMatView2) {
    Mat<int, 4, 4> matrix = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
    auto view = matrix.submatrix<3, 3>( 1, 1 );

    auto g = [] <typename Impl, typename T, MatDim n, MatDim m> ( MatFacade<Impl, T, n, m>& o )
    {
        return o.template submatrix<2, 1>( 0, 0 );
    };

    auto view2 = view.submatrix<2, 1>( 0, 1 );
    auto view3 = g( view2 );

    // check that MatView does not nest
    static_assert( std::same_as<decltype( view2 )::BaseT::Impl, decltype( matrix )> );
    static_assert( std::same_as<decltype( view3 )::BaseT::Impl, decltype( matrix )> );

    EXPECT_EQ( &view2( 0, 0 ), &matrix( 1, 2 ) );
    EXPECT_EQ( &view2( 1, 0 ), &matrix( 2, 2 ) );
}


TEST(MatView, SwapColumns) {
    // Test swapping columns
    Mat<int, 3, 3> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    MatView view(matrix);

    view.swapCols(0, 1);

    EXPECT_EQ(matrix(0, 0), 2);
    EXPECT_EQ(matrix(0, 1), 1);
    EXPECT_EQ(matrix(1, 0), 5);
    EXPECT_EQ(matrix(1, 1), 4);
    EXPECT_EQ(matrix(2, 0), 8);
    EXPECT_EQ(matrix(2, 1), 7);
}


TEST(MatView, ConstAccess) {
    // Test const access
    const Mat<int, 3, 3> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    const MatView view( matrix );

    EXPECT_EQ(view(0, 0), 1);
    EXPECT_EQ(view(1, 1), 5);
    EXPECT_EQ(view(2, 2), 9);

    EXPECT_EQ( &view(0, 0), &matrix(0, 0) );
    EXPECT_EQ( &view(2, 1), &matrix(2, 1) );
    EXPECT_EQ( &view(1, 2), &matrix(1, 2) );
}


TEST(MatView, Addition) {
    Mat<int, 3, 3> mat{
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };
    auto submat = mat.submatrix( 1, 1, 2, 2 );

    EXPECT_EQ( submat.rows(), 2 );
    EXPECT_EQ( submat.cols(), 2 );

    Mat<int, 2, 2> other{
        {10, 15},
        {20, 25}
    };
    submat = submat + other;

    // check modification of the original matrix
    EXPECT_EQ(mat(1, 1), 15);
    EXPECT_EQ(mat(1, 2), 21);
    EXPECT_EQ(mat(2, 1), 28);
    EXPECT_EQ(mat(2, 2), 34);
}

TEST(MatView, CopyConstructor1) {
    Mat<int, 3, 3> mat{
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };
    auto submat = mat.submatrix( 1, 1, 2, 2 );
    auto submat2 = submat;

    for ( Index i = 0; i < submat.rows(); ++i )
        for ( Index j = 0; j < submat.cols(); ++j )
            EXPECT_EQ( &submat2( i, j ), &submat( i, j ) );
}

template<typename T>
auto submatrix_split( T& a )
{
    Index m = a.rows();
    Index n = a.cols();
    Index m2 = (m + 1) / 2;
    Index n2 = (n + 1) / 2;
    return std::tuple{
        a.submatrix(0, 0, m2, n2),
        a.submatrix(0, n2, m2, n - n2),
        a.submatrix(m2, 0, m - m2, n2),
        a.submatrix(m2, n2, m - m2, n - n2)
    };
}

TEST(MatView, CopyConstructor2)
{
    Mat<int> mat = Mat<int, 4, 4>{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
    auto [a, b, c, d] = submatrix_split( mat );
    EXPECT_EQ( &a(0, 0), &mat(0, 0) );
}