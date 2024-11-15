//
// Created by Andrey Aralov on 11/15/24.
//
#pragma once

#include "Mat.hpp"


namespace M
{


template <typename BaseT_, typename T, MatDim n = DynamicMatDim, MatDim m = DynamicMatDim>
class MatView : public MatFacade<MatView<BaseT_, T, n, m>, T, n, m>
{
public:
    using BaseT = std::remove_reference_t<BaseT_>;
private:
    BaseT& base_;
    Index r1, c1;
public:
    using ElemT = T;
    using Parent = MatFacade<MatView<BaseT_, T, n, m>, T, n, m>;
    using Parent::cols;
    using Parent::rows;

    /// If size of the view is known at compile time, this constructor is available
    template <typename Impl, MatDim n1, MatDim m1>
        requires ( !n.dynamic && !m.dynamic )
    explicit MatView( MatFacade<Impl, T, n1, m1>& other, Index r1_ = 0, Index c1_ = 0 ):
        r1( r1_ ), c1( c1_ ),
        base_( other )
    {
        static_assert( n <= n1 && m <= m1 );
    }
    template <typename Impl, MatDim n1, MatDim m1>
        requires ( !n.dynamic && !m.dynamic )
    explicit MatView( const MatFacade<Impl, T, n1, m1>& other, Index r1_ = 0, Index c1_ = 0 ):
        r1( r1_ ), c1( c1_ ),
        base_( other )
    {
        static_assert( n <= n1 && m <= m1 );
    }

    template <typename Impl, MatDim n1, MatDim m1>
    MatView( MatFacade<Impl, T, n1, m1>& other, Index r1_, Index c1_, Index nrows, Index ncols ):
        Parent( nrows, ncols ),
        r1( r1_ ), c1( c1_ ),
        base_( other )
    {
        static_assert( n <= n1 && m <= m1 );
    }
    template <typename Impl, MatDim n1, MatDim m1>
    MatView( const MatFacade<Impl, T, n1, m1>& other, Index r1_, Index c1_, Index nrows, Index ncols ):
        Parent( nrows, ncols ),
        r1( r1_ ), c1( c1_ ),
        base_( other )
    {
        static_assert( n <= n1 && m <= m1 );
    }

    template <typename Impl, MatDim n1, MatDim m1>
    MatView& operator=( const MatFacade<Impl, T, n1, m1>& other )
    {
        static_assert( n1 == n );
        static_assert( m1 == m );
        assert( other.rows() == rows() );
        assert( other.cols() == cols() );
        for ( Index i = 0; i < rows(); ++i )
            for ( Index j = 0; j < cols(); ++j )
                (*this)( i, j ) = other( i, j );
        return *this;
    }


    constexpr T& operator() ( Index i, Index j )
        { return base_( i + r1, j + c1 ); }
    constexpr const T& operator() ( Index i, Index j ) const
        { return base_( i + r1, j + c1 ); }
};


// Simple deduction guides to create a full matrix view
template <typename Impl, typename T, MatDim n, MatDim m>
MatView( const MatFacade<Impl, T, n, m>& ) -> MatView<const MatFacade<Impl, T, n, m>, T, n, m>;

template <typename Impl, typename T, MatDim n, MatDim m>
MatView( MatFacade<Impl, T, n, m>& ) -> MatView<MatFacade<Impl, T, n, m>, T, n, m>;

}