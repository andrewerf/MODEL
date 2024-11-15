//
// Created by Andrey Aralov on 9/16/24.
//
#pragma once

#include <vector>
#include <cassert>
#include <cmath>


namespace M
{


/// Type used to represent dimensions and indexes of matrices
using Index = int;

/// A single matrix dimension
/// It has two constructors:
///     - if ``Index`` is given, then it denotes a dimension which is known at compile time.
///     - the default constructor denotes a dimension which is only known at runtime
/// The usage of this class is to perform static (i.e. compile-time) shape inference where possible.
/// @note This struct is designed to be used only in compile time.
struct MatDim
{
    consteval MatDim( Index x ): val( x ), dynamic( false ) {}
    consteval MatDim(): val( -1 ), dynamic( true ) {}

    Index val;
    bool dynamic;
};
consteval bool operator==( MatDim a, MatDim b )
{
    if ( a.dynamic || b.dynamic )
        return true;
    else
        return a.val == b.val;
}
consteval bool operator<( MatDim a, MatDim b )
{
    if ( a.dynamic || b.dynamic )
        return true;
    else
        return a.val < b.val;
}
consteval bool operator<=( MatDim a, MatDim b )
{
    if ( a.dynamic || b.dynamic )
        return true;
    else
        return a.val <= b.val;
}


/// Just for convenience
constexpr MatDim DynamicMatDim;


/// Everything we need from the implementation of a Matrix
template <typename Mat>
concept MatImplC = requires( Mat mat, const Mat& cmat, Index i, Index j )
{
    typename Mat::ElemT;
    requires ( !std::is_reference_v<typename Mat::ElemT> );
    { mat( i, j ) } -> std::same_as<typename Mat::ElemT&>;
    { cmat( i, j ) } -> std::same_as<const typename Mat::ElemT&>;
};


// forward-declare MatView to use it in `submatrix` method of the facade
template <typename BaseT_, typename T, MatDim n, MatDim m>
class MatView;

/// A base class for the static polymorphism based on CRTP. Represents the interface of a matrix.
template <typename Impl_, typename ElemT_, MatDim n, MatDim m>
class MatFacade
{
protected:
    Index rowsDynamic_, colsDynamic_;
public:
    using ElemT = ElemT_;
    using Impl = Impl_;

    constexpr Impl& impl()
        { return static_cast<Impl&>( *this ); }
    constexpr const Impl& impl() const
        { return static_cast<const Impl&>( *this ); }

    constexpr ElemT& operator() ( Index i, Index j )
        { return impl()( i, j ); }
    constexpr const ElemT& operator() ( Index i, Index j ) const
        { return impl()( i, j ); }


    /// @brief Get submatrix view of the statically known size
    /// @note This method does not involve copying of elements
    /// @tparam nrows Number of rows
    /// @tparam ncols Number of columns
    /// @param r1 Index of the first row of the submatrix in the matrix
    /// @param c1 Index of the first column of the submatrix in the matrix
    template <MatDim nrows, MatDim ncols>
        requires ( !nrows.dynamic && !ncols.dynamic )
    constexpr MatView<MatFacade<Impl_, ElemT, n, m>, ElemT, nrows, ncols>
        submatrix( Index r1, Index c1 )
    {
        static_assert( nrows <= n );
        static_assert( ncols <= m );
        assert( r1 >= 0 );
        assert( c1 >= 0 );
        assert( nrows.val + r1 <= rows() );
        assert( ncols.val + c1 <= cols() );
        return MatView<MatFacade<Impl_, ElemT, n, m>, ElemT, nrows, ncols>( *this, r1, c1 );
    }

    /// @brief Get submatrix view of the dynamically known size
    /// @note This method does not involve copying of elements
    /// @param r1 Index of the first row of the submatrix in the matrix
    /// @param c1 Index of the first column of the submatrix in the matrix
    /// @param nrows_ Number of rows
    /// @param ncols_ Number of columns
    template <MatDim nrows = DynamicMatDim, MatDim ncols = DynamicMatDim>
    constexpr MatView<MatFacade<Impl_, ElemT, n, m>, ElemT, nrows, ncols>
        submatrix( Index r1, Index c1, Index nrows_, Index ncols_ )
    {
        static_assert( nrows <= n );
        static_assert( ncols <= m );
        assert( r1 >= 0 );
        assert( c1 >= 0 );
        assert( nrows_ + r1 < rows() );
        assert( ncols_ + c1 < cols() );
        return MatView<MatFacade<Impl_, ElemT, n, m>, ElemT, nrows, ncols>( *this, r1, c1, nrows_, ncols_ );
    }


    constexpr Index rows() const
    {
        if constexpr ( n.dynamic )
            return rowsDynamic_;
        else
            return n.val;
    }
    constexpr Index cols() const
    {
        if constexpr ( m.dynamic )
            return colsDynamic_;
        else
            return m.val;
    }

    constexpr void swapCols( Index c1, Index c2 )
    {
        for ( Index r = 0; r < rows(); ++r )
            std::swap( (*this)( r, c1 ), (*this)( r, c2 ) );
    }
    constexpr void swapRows( Index r1, Index r2 )
    {
        for ( Index c = 0; c < cols(); ++c )
            std::swap( (*this)( r1, c ), (*this)( r2, c ) );
    }

    /// This function is called to ensure correct usage of the CRTP
    constexpr void checkInvariants()
    {
        static_assert( MatImplC<Impl> );
        static_assert( std::same_as<typename Impl::ElemT, ElemT> );
        static_assert( std::derived_from<Impl_, MatFacade> );
        static_assert( ( n.dynamic || n.val >= 0 ) && ( m.dynamic || m.val >= 0 ) );
    };

    constexpr MatFacade()
        requires ( !n.dynamic && !m.dynamic )
    {
        checkInvariants();
    }
    constexpr MatFacade( Index rows, Index cols ):
        rowsDynamic_( rows ), colsDynamic_( cols )
    {
        if constexpr ( !n.dynamic )
            assert( rowsDynamic_ == n.val );
        if constexpr ( !m.dynamic )
            assert( colsDynamic_ == m.val );
        checkInvariants();
    }
};


/// Simple matrix implementation
template <typename T, MatDim n = DynamicMatDim, MatDim m = DynamicMatDim>
struct Mat : public MatFacade<Mat<T, n, m>, T, n, m>
{
    using ElemT = T;
    using Parent = MatFacade<Mat<T, n, m>, T, n, m>;
    using Parent::cols;
    using Parent::rows;

    /// If both dimensions are known at compile time, we know how many elements to allocate
    Mat()
        requires ( !n.dynamic && !m.dynamic )
    { a.resize( n.val * m.val ); }

    Mat( Index rows, Index cols ):
        Parent( rows, cols )
    { a.resize( rows * cols ); }

    /// Copy constructor from any Mat-like type
    template <typename Impl, MatDim n1, MatDim m1>
    Mat( const MatFacade<Impl, T, n1, m1>& other ):
        Parent( other.rows(), other.cols() )
    {
        static_assert( n1 == m1 );
        a.resize( rows() * cols() );
        for ( Index i = 0; i < rows(); ++i )
            for ( Index j = 0; j < cols(); ++j )
                (*this)( i, j ) = other( i, j );
    }

    /// Nested initializer list constructor
    /// @tparam cols_ Pack of sizes of rows of the matrix. Each element in the pack must be equal to m (enforced by types)
    /// @param init Pack of initializer lists, containing rows of the matrix being constructed
    /// For example, this code produces a 3x3 matrix, with rows {1, 2, 3}, {4, 5, 6} and {7, 8, 9} respectively.
    /// ```
    ///     Mat<int, 3, 3> mat{
    ///        { 1, 2, 3 },
    ///        { 4, 5, 6 },
    ///        { 7, 8, 9 }
    ///    };
    /// ```
    template <Index ...cols_>
        requires ( !n.dynamic && !m.dynamic )
    constexpr Mat( const ElemT (&...init)[cols_] )
    {
        static_assert( sizeof...( cols_ ) == n );
        static_assert( ( ( cols_ == m ) && ... ) );
        a.resize( rows() * cols() );
        Index i = 0;
        auto fillRow = [this, &i] ( auto row ) {
            for ( Index j = 0; j < cols(); ++j )
                (*this)( i, j ) = row[j];
            ++i;
        };
        ( fillRow( init ) , ... );
    }

    static Mat Identity( Index rows, Index cols )
        requires ( m == n )
    {
        assert( rows == cols );
        Mat res( rows, cols );
        for ( Index i = 0; i < rows; ++i )
            res( i, i ) = 1;
        return res;
    }

    constexpr T& operator() ( Index i, Index j )
        { return a[i * cols() + j]; }
    constexpr const T& operator() ( Index i, Index j ) const
        { return a[i * cols() + j]; }

    std::vector<T> a;
};


namespace detail
{

template <typename MatA, typename MatB, typename MatC>
void productImpl( const MatA& a, const MatB& b, MatC& c )
{
    assert( a.cols() == b.rows() );
    for ( Index i = 0; i < a.rows(); ++i )
    {
        for ( Index j = 0; j < b.cols(); ++j )
        {
            c( i, j ) = 0;
            for ( Index t = 0; t < a.cols(); ++t )
            {
                c( i, j ) += a( i, t ) * b( t, j );
            }
        }
    }
}

template <typename MatA, typename MatB, typename MatC>
void sumImpl( const MatA& a, const MatB& b, MatC& c )
{
    assert( a.cols() == b.cols() );
    assert( a.rows() == b.rows() );
    assert( a.cols() == c.cols() );
    assert( a.rows() == c.rows() );
    for ( Index i = 0; i < a.rows(); ++i )
    {
        for ( Index j = 0; j < a.cols(); ++j )
        {
            c( i, j ) = a( i, j ) + b( i, j );
        }
    }
}

}


template <typename T,
          typename Impl1, MatDim n1, MatDim m1,
          typename Impl2, MatDim n2, MatDim m2>
Mat<T, n1, m2> operator*( const MatFacade<Impl1, T, n1, m1>& a, const MatFacade<Impl2, T, n2, m2>& b )
{
    static_assert( m1 == n2 );
    Mat<T, n1, m2> ret( a.rows(), b.cols() );
    detail::productImpl( a, b, ret );
    return ret;
}

template <typename T,
        typename Impl1, MatDim n1, MatDim m1,
        typename Impl2, MatDim n2, MatDim m2>
Mat<T, n1, m1> operator+( const MatFacade<Impl1, T, n1, m1>& a, const MatFacade<Impl2, T, n2, m2>& b )
{
    static_assert( n1 == n2 );
    static_assert( m1 == m2 );
    Mat<T, n1, m1> ret( a.rows(), a.cols() );
    detail::sumImpl( a, b, ret );
    return ret;
}


template <typename T, typename Impl, MatDim n, MatDim m>
Mat<T, n, m> operator+( const MatFacade<Impl, T, n, m>& a, T val )
{
    Mat<T, n, m> ret( a.rows(), a.cols() );
    for ( size_t i = 0; i < a.rows(); ++i )
        for ( size_t j = 0; j < a.cols(); ++j )
            ret( i, j ) = a( i, j ) + val;
    return ret;
}
template <typename T, typename Impl, MatDim n, MatDim m>
Mat<T, n, m> operator+( T val, const MatFacade<Impl, T, n, m>& a )
{
    Mat<T, n, m> ret( a.rows(), a.cols() );
    for ( size_t i = 0; i < a.rows(); ++i )
        for ( size_t j = 0; j < a.cols(); ++j )
            ret( i, j ) = a( i, j ) + val;
    return ret;
}

template <typename T, typename Impl, MatDim n, MatDim m>
Mat<T, n, m> operator-( const MatFacade<Impl, T, n, m>& a, T val )
{
    Mat<T, n, m> ret( a.rows(), a.cols() );
    for ( size_t i = 0; i < a.rows(); ++i )
        for ( size_t j = 0; j < a.cols(); ++j )
            ret( i, j ) = a( i, j ) - val;
    return ret;
}
template <typename T, typename Impl, MatDim n, MatDim m>
Mat<T, n, m> operator-( T val, const MatFacade<Impl, T, n, m>& a )
{
    Mat<T, n, m> ret( a.rows(), a.cols() );
    for ( size_t i = 0; i < a.rows(); ++i )
        for ( size_t j = 0; j < a.cols(); ++j )
            ret( i, j ) = val - a( i, j );
    return ret;
}


template <typename T, typename Impl, MatDim n, MatDim m>
Mat<T, n, m> operator*( const MatFacade<Impl, T, n, m>& a, T val )
{
    Mat<T, n, m> ret( a.rows(), a.cols() );
    for ( size_t i = 0; i < a.rows(); ++i )
        for ( size_t j = 0; j < a.cols(); ++j )
            ret( i, j ) = a( i, j ) * val;
    return ret;
}
template <typename T, typename Impl, MatDim n, MatDim m>
Mat<T, n, m> operator*( T val, const MatFacade<Impl, T, n, m>& a )
{
    Mat<T, n, m> ret( a.rows(), a.cols() );
    for ( size_t i = 0; i < a.rows(); ++i )
        for ( size_t j = 0; j < a.cols(); ++j )
            ret( i, j ) = a( i, j ) * val;
    return ret;
}


}