//
// Created by Andrey Aralov on 9/17/24.
//
#pragma once

#include <limits>


template <typename T, Index n>
struct LU_Result
{
    Mat<T, n, n> L, U;
};

template <typename T, Index n>
LU_Result<T, n> LU( Mat<T, n, n> mat )
{
    auto L = Mat<T, n, n>::Identity();

    for ( Index k = 0; k < n; ++k )
    {
        for ( Index j = k + 1; j < n; ++j )
        {
            const T alpha = -mat( j, k ) / mat( k, k );
            L( j, k ) = -alpha;

            for ( Index i = 0; i < n; ++i )
                mat( j, i ) = mat( j, i ) + alpha * mat( k, i );
        }
    }

    return { L, mat };
}


template <typename T, Index n>
struct PLUQ_Result
{
    Mat<T, n, n> P, L, U, Q;
};

template <typename T, Index n>
PLUQ_Result<T, n> PLUQ( Mat<T, n, n> mat )
{
    using TM = Mat<T, n, n>;

    auto argmax = [&mat] ( Index k0 )
    {
        typename TM::Coord ret;
        T val = std::numeric_limits<T>::min();
        for ( Index row = k0; row < n; ++row )
        {
            for ( Index col = k0; col < n; ++col )
            {
                if ( mat( row, col ) > val )
                {
                    val = mat( row, col );
                    ret = { row, col };
                }
            }
        }

        return ret;
    };

    TM P = TM::Identity();
    TM L = TM::Identity();
    TM Q = TM::Identity();

    for ( Index k = 0; k < n - 1; ++k )
    {
        // find the max element in the submatrix
        auto [mxrow, mxcol] = argmax( k );

        // swap cols and rows in the original matrix, so that mat(k, k) is the max element
        mat.swapCols( k, mxcol );
        mat.swapRows( k, mxrow );

        // compensate for swapping in the original matrix, by swapping surroundings
        // Remember that:
        // - multiplying on the right == swapping cols and
        // - multiplying on the left == swapping rows
        Q.swapRows( k, mxcol );

        L.swapRows( k, mxrow );
        L.swapCols( k, mxrow );

        P.swapCols( k, mxrow );


        for ( Index j = k + 1; j < n; ++j )
        {
            const T alpha = -mat( j, k ) / mat( k, k );
            L( j, k ) = -alpha;

            for ( Index i = 0; i < n; ++i )
                mat( j, i ) = mat( j, i ) + alpha * mat( k, i );
        }
    }

    return { P, L, mat, Q };
}
