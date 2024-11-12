//
// Created by Andrey Aralov on 9/17/24.
//
#pragma once

#include <limits>
#include <Mat.hpp>


namespace M
{


template <typename T, MatDim n>
struct LU_Result
{
    Mat<T, n, n> L, U;
};

template <typename T, typename Impl, MatDim n, MatDim m>
LU_Result<T, n> LU( const MatFacade<Impl, T, n, m>& mat )
{
    static_assert( n == m );
    assert( mat.rows() == mat.cols() );
    auto L = Mat<T, n, n>::Identity( mat.rows(), mat.cols() );
    Mat<T, n, n> U = mat;

    for ( Index k = 0; k < mat.rows(); ++k )
    {
        for ( Index j = k + 1; j < mat.rows(); ++j )
        {
            const T alpha = -U( j, k ) / U( k, k );
            L( j, k ) = -alpha;

            for ( Index i = k; i < mat.rows(); ++i )
                U( j, i ) = U( j, i ) + alpha * U( k, i );
        }
    }

    return { L, U };
}


template <typename T, MatDim n>
struct PLUQ_Result
{
    Mat<T, n, n> P, L, U, Q;
};

template <typename T, typename Impl, MatDim n, MatDim m>
PLUQ_Result<T, n> PLUQ( const MatFacade<Impl, T, n, m>& matOrig )
{
    static_assert( n == m );
    assert( matOrig.rows() == matOrig.cols() );
    using TM = Mat<T, n, n>;
    TM U = matOrig;

    auto argmax = [&U] ( Index k0 )
    {
        std::pair<Index, Index> ret;
        T val = std::numeric_limits<T>::min();
        for ( Index row = k0; row < U.rows(); ++row )
        {
            for ( Index col = k0; col < U.rows(); ++col )
            {
                if ( U( row, col ) > val )
                {
                    val = U( row, col );
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
        // find the max element in the subUrix
        auto [mxrow, mxcol] = argmax( k );

        // swap cols and rows in the original U, so that U(k, k) is the max element
        U.swapCols( k, mxcol );
        U.swapRows( k, mxrow );

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

    return { P, L, U, Q };
}

}