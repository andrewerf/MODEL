// MU4IN901 - MODEL - Implementation project
//
// Authors:
// Andrei Aralov (andrei.aralov@etu.sorbonne-universite.fr)
// Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
// Ruchi  Thareja (ruchi.thareja@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// LU and PLUQ decompositions; matrix inversion by linear system solving.

#pragma once

#include <limits>
#include <Mat.hpp>
#include <MatView.hpp>

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
/// @return If L is not lower triangular, returns garbage. If L is degenerate, return nullopt
template <typename T, typename Impl1, typename Impl2, MatDim n, MatDim m, MatDim n1, MatDim m1>
std::optional<Mat<T, n, 1>> solveLower( const MatFacade<Impl1, T, n, m>& L, const MatFacade<Impl2, T, n1, m1>& y )
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
        T f = y( k, 0 );
        for ( Index i = 0; i < k; ++i )
            f -= L( k, i ) * ret( i, 0 );
        auto t = L( k, k );
        if ( t == 0 )
            return std::nullopt;
        f /= t;
        ret( k, 0 ) = f;
    }

    return ret;
}

/// Returns solution of the linear system U*x = y where U is an upper triangular matrix
/// @return If U is not upper triangular, returns garbage. If U is degenerate, return nullopt
template <typename T, typename Impl1, typename Impl2, MatDim n, MatDim m, MatDim n1, MatDim m1>
std::optional<Mat<T, n, 1>> solveUpper( const MatFacade<Impl1, T, n, m>& U, const MatFacade<Impl2, T, n1, m1>& y )
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
        T f = y( rows - k, 0 );
        for ( Index i = 1; i < k; ++i )
            f -= U( rows - k, rows - i ) * ret( rows - i, 0 );
        auto t = U( rows - k, rows - k );
        if ( t == 0 )
            return std::nullopt;
        f /= t;
        ret( rows - k, 0 ) = f;
    }

    return ret;
}


template <typename T, typename Impl, MatDim n_, MatDim m_>
std::optional<Mat<T, n_, m_>> inverseLU( const MatFacade<Impl, T, n_, m_>& mat )
{
    static_assert( n_ == m_, "Only square matrices can be inverted" );
    assert( mat.rows() == mat.cols() );

    const auto maybeLU = LU( mat );
    if ( !maybeLU )
        return {};

    const auto n = mat.rows();
    const auto& [L, U] = *maybeLU;

    // for each column yi of the identity matrix
    // L * U * x = yi, denote z = U*x
    // => L * z = yi => find z
    // and then U * x = z => find x

    Mat<T, n_, n_> ret( n, n );
    Mat<T, n_, 1> col( n, 1 );
    for ( Index i = 0; i < n; ++i )
    {
        if ( i > 0 )
            col( i - 1, 0 ) = 0;
        col( i, 0 ) = T( 1 );

        auto x = solveLower( L, col ).and_then( [&U] ( auto&& z ) { return solveUpper( U, z ); } );
        if ( !x )
            return std::nullopt;
        ret.submatrix( 0, i, n, 1 ) = *x;
    }

    return ret;
}

}