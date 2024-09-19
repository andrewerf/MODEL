//
// Created by Andrey Aralov on 9/16/24.
//
#pragma once


using Index = int;
struct NoInit{};


template <typename T, Index n, Index m>
class Mat
{
    std::array<T, n*m> a; // a[i][j] == a[i * m + j]

public:

    Mat()
    {
        for ( Index i = 0; i < n*m; ++i )
            a[i] = T( 0 );
    }

    explicit Mat( NoInit )
    {}

    T& operator() ( Index row, Index col )
    {
        return a[row * m + col];
    }

    const T& operator() ( Index row, Index col ) const
    {
        assert( row < n );
        assert( col < m );
        return a[row * m + col];
    }
};


template <typename T, Index n1, Index m, Index n2>
Mat<T, n1, n2> operator*( const Mat<T, n1, m>& a, const Mat<T, m, n2>& b )
{
    Mat<T, n1, n2> res;
    for ( Index i = 0; i < n1; ++i )
    {
        for ( Index j = 0; j < n2; ++j )
        {
            for ( Index k = 0; k < m; ++k )
                res( i, j ) += a( i, k ) * b( k, j );
        }
    }
    return res;
}
