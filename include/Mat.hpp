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
    struct Coord
    {
        Index row;
        Index col;
    };


    Mat()
    {
        for ( Index i = 0; i < n*m; ++i )
            a[i] = T( 0 );
    }

    explicit Mat( NoInit )
    {}


    constexpr Index rows()
    {
        return n;
    }
    constexpr Index cols()
    {
        return m;
    }


    static Mat<T, n, n> Identity()
    {
        Mat<T, n, n> r;
        for ( Index i = 0; i < n; ++i )
            r( i, i ) = T( 1 );
        return r;
    }

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


    template <typename NewT>
    Mat<NewT, n, m> cast()
    {
        Mat<NewT, n, m> ret;
        for ( Index row = 0; row < n; ++row )
        {
            for ( Index col = 0; col < m; ++col )
                ret( row, col ) = static_cast<NewT>( (*this)( row, col ) );
        }
        return ret;
    }


    void swapRows( Index row1, Index row2 )
    {
        for ( Index c = 0; c < cols(); ++c )
            std::swap( (*this)( row1, c ), (*this)( row2, c ) );
    }

    void swapCols( Index col1, Index col2 )
    {
        for ( Index r = 0; r < rows(); ++r )
            std::swap( (*this)( r, col1 ), (*this)( r, col2 ) );
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
