//
// Created by Andrey Aralov on 9/16/24.
//
#pragma once

#include <array>
#include <cassert>
#include <cmath>


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

    const T& operator[]( Index i ) const
        requires ( n == 1 || m == 1 )
    {
        if constexpr ( n == 1 )
            return ( *this )( 0, i );
        else
            return ( *this ) ( i, 0 );
    }
    T& operator[]( Index i )
        requires ( n == 1 || m == 1 )
    {
        if constexpr ( n == 1 )
            return ( *this )( 0, i );
        else
            return ( *this ) ( i, 0 );
    }

    T norm2() const
    {
        T res = 0;
        for ( auto& v : a )
            res += v * v;
        return std::sqrt( res );
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

    void transpose()
        requires (n == m)
    {
        for ( Index row = 0; row < n; ++row )
            for ( Index col = row + 1; col < m; ++col )
                std::swap( (*this)( row, col ), (*this)( col, row ) );
    }

    Mat<T, m, n> transposed()
    {
        if constexpr ( n == m )
        {
            auto r = *this;
            r.transpose();
            return r;
        }
        else
        {
            Mat<T, m, n> r;
            for ( Index row = 0; row < n; ++row )
                for ( Index col = 0; col < m; ++col )
                    r( col, row ) = (*this)( row, col );
            return r;
        }
    }

    Mat<T, n, m> operator+=( const Mat<T, n, m>& other )
    {
        for ( Index i = 0; i < n; ++i )
            for ( Index j = 0; j < m; ++j )
                (*this)( i, j ) += other( i, j );
        return *this;
    }

    Mat<T, n, m> operator-=( const Mat<T, n, m>& other )
    {
        for ( Index i = 0; i < n; ++i )
            for ( Index j = 0; j < m; ++j )
                (*this)( i, j ) -= other( i, j );
        return *this;
    }

    Mat<T, n, m> operator*=( T x )
    {
        for ( Index i = 0; i < n; ++i )
            for ( Index j = 0; j < m; ++j )
                (*this)( i, j ) *= x;
        return *this;
    }
    Mat<T, n, m> operator/=( T x )
    {
        for ( Index i = 0; i < n; ++i )
            for ( Index j = 0; j < m; ++j )
                (*this)( i, j ) /= x;
        return *this;
    }

    Mat<T, 1, m> row( Index i ) const
    {
        Mat<T, 1, m> ret( NoInit{} );
        for ( Index k = 0; k < m; ++k )
            ret[k] = (*this)( i, k );
        return ret;
    }

    Mat<T, n, 1> col( Index i ) const
    {
        Mat<T, n, 1> ret( NoInit{} );
        for ( Index k = 0; k < n; ++k )
            ret[k] = (*this)( k, i );
        return ret;
    }

    void setCol( Index i, const Mat<T, n, 1>& col )
    {
        for ( Index k = 0; k < n; ++k )
            (*this)( k, i ) = col[k];
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

template <typename T, Index n, Index m>
Mat<T, n, m> operator*( Mat<T, n, m> mat, T x )
{
    mat *= x;
    return mat;
}
template <typename T, Index n, Index m>
Mat<T, n, m> operator*( T x, Mat<T, n, m> mat )
{
    mat *= x;
    return mat;
}
template <typename T, Index n, Index m>
Mat<T, n, m> operator/( Mat<T, n, m> mat, T x )
{
    mat /= x;
    return mat;
}


template <typename T, Index n, Index m>
Mat<T, n, m> operator+( Mat<T, n, m> a, const Mat<T, n, m>& b )
{
    a += b;
    return a;
}
template <typename T, Index n, Index m>
Mat<T, n, m> operator-( Mat<T, n, m> a, const Mat<T, n, m>& b )
{
    a -= b;
    return a;
}


template <typename T, Index n>
T dot( const Mat<T, 1, n>& a, const Mat<T, n, 1>& b )
{
    T res = 0;
    for ( Index k = 0; k < n; ++k )
        res += a[k] * b[k];
    return res;
}
