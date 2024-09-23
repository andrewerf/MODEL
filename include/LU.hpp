//
// Created by Andrey Aralov on 9/17/24.
//
#pragma once

template <typename T, Index n, Index m>
void print( const Mat<T, n, m>& mat )
{
    for ( Index i = 0; i < n; ++i )
    {
        for ( Index j = 0; j < m; ++j )
        {
            std::cout << mat( i, j ) << ' ';
        }
        std::cout << '\n';
    }
}


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
            const auto alpha = -mat( j, k ) / mat( k, k );
            L( j, k ) = -alpha;

            for ( Index i = 0; i < n; ++i )
                mat( j, i ) = mat( j, i ) + alpha * mat( k, i );
        }
    }

    return { L, mat };
}
