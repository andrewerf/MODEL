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

