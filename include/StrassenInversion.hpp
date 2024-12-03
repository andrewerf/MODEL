//
// Created by Andrey Aralov on 11/23/24.
//
#pragma once

#include "Mat.hpp"
#include "MatView.hpp"

#include <optional>

namespace M
{


template <typename T, typename Impl, MatDim n, MatDim m, typename F = std::multiplies<>>
std::optional<Mat<T, n, m> > inverseStrassen( const MatFacade<Impl, T, n, m>& mat, F mult = {} )
{
    static_assert( n == m );
    assert( mat.rows() == mat.cols() );

    if ( mat.rows() == 1 )
    {
        if (mat(0, 0) == 0) {
            return {};
        } else {
            Mat<T, n, m> ret = mat;
            ret( 0, 0 ) = T( 1 ) / ret( 0, 0 );
            return ret;
        }
    }

    Index r2 = mat.rows() / 2;
    Index c2 = mat.cols() / 2;
    const auto a = mat.submatrix( 0, 0, r2, c2 );
    const auto c = mat.submatrix( r2, 0, mat.rows() - r2, c2 );
    const auto b = mat.submatrix( 0, c2, r2, mat.cols() - c2 );
    const auto d = mat.submatrix( r2, c2, mat.rows() - r2, mat.cols() - c2 );

    const auto a_inv = inverseStrassen( a );
    if ( !a_inv.has_value() )
        return {};

    const auto e = a_inv.value();
    const auto Z = d - mult( mult( c, e ), b );
    const auto Z_inv = inverseStrassen( Z );
    if ( !Z_inv.has_value() )
        return {};
    const auto t = Z_inv.value();

    const auto y = -mult( mult( e, b ), t );
    const auto z = -mult( mult( t, c ), e );
    const auto x = e - mult( y, mult( c, e ) );

    Mat<T, n, m> ret( mat.rows(), mat.cols() );
    ret.submatrix( 0, 0, r2, c2 ) = x;
    ret.submatrix( r2, 0, mat.rows() - r2, c2 ) = z;
    ret.submatrix( 0, c2, r2, mat.cols() - c2 ) = y;
    ret.submatrix( r2, c2, mat.rows() - r2, mat.cols() - c2 ) = t;
    return ret;
}


}