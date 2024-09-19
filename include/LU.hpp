//
// Created by Andrey Aralov on 9/17/24.
//
#pragma once



template <typename T, Index n>
struct LU_Result
{
    Mat<T, n, n> L, U;
};

template <typename T, Index n>
LU_Result<T, n> LU( Mat<T, n, n> mat )
{
    Mat<T, n, n> L;
    for ( Index i = 0; i < n; ++i )
        L( i, i ) = T( 1 );

    for ( Index k = 0; k < n; ++k )
    {
        Mat<T, n, n> L_inv_0;
        for ( Index j = 0; j < n; ++j )
            L_inv_0( j, j ) = 1.f;

        for ( Index j = k + 1; j < n; ++j )
        {
            const auto alpha = -mat( j, k ) / mat( k, k );
            L_inv_0( j, j ) = 1;
            L_inv_0( j, k ) = alpha;

            for ( Index i = 0; i < n; ++i )
                mat( j, i ) = mat( j, i ) + alpha * mat( k, i );
        }

        Mat<T, n, n> L_0;
        for ( Index i = 0; i < n; ++i )
            L_0( i, i ) = T( 1 );
        for ( Index i = 1; i < n; ++i )
        {
            for ( Index j = 0; j < i; ++j )
            {
                L_0( i, j ) = -L_inv_0( i, j );
            }
        }

        L = L * L_0;
    }

    return { L, mat };
}
