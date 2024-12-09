// MU4IN901 - MODEL - Implementation project
//
// Authors:
// Andrei Aralov (andrei.aralov@etu.sorbonne-universite.fr)
// Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
// Ruchi  Thareja (ruchi.thareja@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Implementation of the Strassen algorithm for matrix inversion.

#pragma once

#include "Mat.hpp"
#include "MatView.hpp"
#include "LU.hpp"

#include <optional>

namespace M
{


template <typename T, typename Impl, MatDim n, MatDim m, typename F = std::multiplies<>>
std::optional<Mat<T, n, m> > inverseStrassen( const MatFacade<Impl, T, n, m>& mat, F mult = {}, Index leafSize = 512 )
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

    if ( mat.rows() <= leafSize )
        return inverseLU( mat );

    Index r2 = mat.rows() / 2;
    Index c2 = mat.cols() / 2;
    const auto a = mat.submatrix( 0, 0, r2, c2 );
    const auto c = mat.submatrix( r2, 0, mat.rows() - r2, c2 );
    const auto b = mat.submatrix( 0, c2, r2, mat.cols() - c2 );
    const auto d = mat.submatrix( r2, c2, mat.rows() - r2, mat.cols() - c2 );

    const auto a_inv = inverseStrassen( a, mult );
    if ( !a_inv )
        return {};

    const auto& e = *a_inv;
    const auto ce = mult( c, e );
    const auto Z = d - mult( ce, b );
    const auto Z_inv = inverseStrassen( Z, mult );
    if ( !Z_inv )
        return {};
    const auto& t = *Z_inv;

    const auto y = -mult( mult( e, b ), t );
    const auto z = -mult( t, ce );
    const auto x = e - mult( y, ce );

    Mat<T, n, m> ret( mat.rows(), mat.cols() );
    ret.submatrix( 0, 0, r2, c2 ) = x;
    ret.submatrix( r2, 0, mat.rows() - r2, c2 ) = z;
    ret.submatrix( 0, c2, r2, mat.cols() - c2 ) = y;
    ret.submatrix( r2, c2, mat.rows() - r2, mat.cols() - c2 ) = t;
    return ret;
}


}