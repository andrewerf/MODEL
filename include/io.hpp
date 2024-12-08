// MU4IN901 - MODEL - Implementation project
//
// Authors:
// Andrei Aralov (andrei.aralov@etu.sorbonne-universite.fr)
// Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
// Ruchi  Thareja (ruchi.thareja@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Serialization of a matrix, as text, to an output stream.

#pragma once

#include <iostream>

#include <Mat.hpp>

namespace M
{

template <typename Impl, typename T, MatDim n, MatDim m>
std::ostream& operator<<( std::ostream& os, const MatFacade<Impl, T, n, m>& mat )
{
    for ( Index i = 0; i < mat.rows(); ++i )
    {
        for ( Index j = 0; j < mat.cols(); ++j )
        {
            os << mat( i, j ) << ' ';
        }
        os << '\n';
    }
    return os;
}

}

