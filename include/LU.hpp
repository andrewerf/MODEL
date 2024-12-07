//
// Created by Andrey Aralov on 9/17/24.
//
#pragma once

#include <limits>
#include <Mat.hpp>

#include <optional>

namespace M
{


template <typename T, MatDim n>
struct LU_Result
{
    Mat<T, n, n> L, U;
};

template <typename T, typename Impl, MatDim n, MatDim m>
std::optional<LU_Result<T, n> > LU( const MatFacade<Impl, T, n, m>& mat )
{
    static_assert( n == m );
    assert( mat.rows() == mat.cols() );
    auto L = Mat<T, n, n>::Identity( mat.rows(), mat.cols() );
    Mat<T, n, n> U = mat;

    for ( Index k = 0; k < mat.rows(); ++k )
    {
        for ( Index j = k + 1; j < mat.rows(); ++j )
        {
            if (U( k, k ) == 0)
            {
                // A null pivot was encountered, abort decomposition.
                return {};
            }
            const T alpha = -U( j, k ) / U( k, k );
            L( j, k ) = -alpha;

            for ( Index i = k; i < mat.rows(); ++i )
                U( j, i ) = U( j, i ) + alpha * U( k, i );
        }
    }

    return { { L, U } };
}


template <typename T, MatDim n>
struct PLUQ_Result
{
    Mat<T, n, n> P, L, U, Q;
};

template <typename T, typename Impl, MatDim n, MatDim m>
std::optional<PLUQ_Result<T, n> > PLUQ( const MatFacade<Impl, T, n, m>& matOrig )
{
    static_assert( n == m );
    assert( matOrig.rows() == matOrig.cols() );
    using TM = Mat<T, n, n>;
    TM U = matOrig;

    auto argmax = [&U] ( Index k0 )
    {
        std::pair<Index, Index> ret;
        T val = std::numeric_limits<T>::lowest();
        for ( Index row = k0; row < U.rows(); ++row )
        {
            for ( Index col = k0; col < U.rows(); ++col )
            {
                if ( fabs(U( row, col )) > val )
                {
                    val = fabs(U( row, col ));
                    ret = { row, col };
                }
            }
        }

        return ret;
    };

    TM P = TM::Identity( matOrig.rows(), matOrig.cols() );
    TM L = P;
    TM Q = P;

    for ( Index k = 0; k < matOrig.rows() - 1; ++k )
    {
        // find the max element in the subMatrix
        auto [mxrow, mxcol] = argmax( k );
    
        // swap cols and rows in the original U, so that U(k, k) is the max element
        U.swapCols( k, mxcol );
        U.swapRows( k, mxrow );

        // the largest element, in absolute value, is 0. It cannot be used as
        // a pivot: abort the decomposition.
        if ( U(k, k) == 0 )
            return {};

        // compensate for swapping in the original U by swapping surroundings
        // Remember that:
        // - multiplying on the right == swapping cols and
        // - multiplying on the left == swapping rows
        Q.swapRows( k, mxcol );

        L.swapRows( k, mxrow );
        L.swapCols( k, mxrow );

        P.swapCols( k, mxrow );


        for ( Index j = k + 1; j < matOrig.rows(); ++j )
        {
            const T alpha = -U( j, k ) / U( k, k );
            L( j, k ) = -alpha;

            for ( Index i = 0; i < matOrig.rows(); ++i )
                U( j, i ) = U( j, i ) + alpha * U( k, i );
        }
    }

    return { { P, L, U, Q } };
}


/// Returns solution of the linear system L*x = y where L is a lower triangular matrix
/// @note If L is not lower triangular, returns garbage
template <typename T, typename Impl, MatDim n, MatDim m, MatDim n1, MatDim m1>
Mat<T, n, 1> solveLower( const MatFacade<Impl, T, n, m>& L, const MatFacade<Impl, T, n1, m1>& y )
{
    static_assert( n == m, "Only square matrices are supported" );
    assert( L.rows() == L.cols() );
    static_assert( m == n1, "Dimension of y should match that of L" );
    assert( L.cols() == y.rows() );
    static_assert( m1 == 1, "y is a column-vector" );
    assert( y.cols() == 1 );

    Mat<T, n, 1> ret( L.rows(), 1 );
    for ( Index k = 0; k < L.rows(); ++k )
    {
        auto& f = ret( k, 0 );
        f = y( k, 0 );
        for ( Index i = 0; i < k; ++i )
            f -= L( k, i ) * ret( i, 0 );
        f /= L( k, k );
    }

    return ret;
}

/// Returns solution of the linear system U*x = y where U is an upper triangular matrix
/// @note If U is not upper triangular, returns garbage
template <typename T, typename Impl, MatDim n, MatDim m, MatDim n1, MatDim m1>
Mat<T, n, 1> solveUpper( const MatFacade<Impl, T, n, m>& U, const MatFacade<Impl, T, n1, m1>& y )
{
    static_assert( n == m, "Only square matrices are supported" );
    assert( U.rows() == U.cols() );
    static_assert( m == n1, "Dimension of y should match that of U" );
    assert( U.cols() == y.rows() );
    static_assert( m1 == 1, "y is a column-vector" );
    assert( y.cols() == 1 );

    const auto rows = U.rows();
    Mat<T, n, 1> ret( rows, 1 );
    for ( Index k = 1; k <= rows; ++k )
    {
        auto& f = ret( rows - k, 0 );
        f = y( rows - k, 0 );
        for ( Index i = 1; i < k; ++i )
            f -= U( rows - k, rows - i ) * ret( rows - i, 0 );
        f /= U( rows - k, rows - k );
    }

    return ret;
}

}