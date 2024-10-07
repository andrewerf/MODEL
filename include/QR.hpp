//
// Created by Andrey Aralov on 10/7/24.
//
#pragma once

#include "Mat.hpp"
#include <cmath>

template <typename T, Index m, Index n>
struct QR_FullResult
{
    Mat<T, m, m> Q;
    Mat<T, m, n> R;
};


template <typename T, Index m, Index n>
QR_FullResult<T, m, n> QR_Given( const Mat<T, m, n>& A )
{
    static_assert( m >= n, "Number of rows must be greater or equal to the number of columns" );

    auto Q = Mat<T, m, m>::Identity();
    auto R = A;

    for ( Index j = 0; j < n; ++j )
    {
        for ( Index i = j + 1; i < m; ++i )
        {
            auto den = std::sqrt( R(j, j)*R(j, j) + R(i, j)*R(i, j) );
            auto cos = R(j, j) / den;
            auto sin = R(i, j) / den;

            auto G = Mat<T, m, m>::Identity();
            G(i, i) = cos;
            G(j, j) = cos;
            G(i, j) = -sin;
            G(j, i) = sin;

            R = G * R;
            G.transpose();
            Q = Q * G;
        }
    }

    return { Q, R };
}
