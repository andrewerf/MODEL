//
// Created by Andrey Aralov on 12/1/24.
//
#pragma once
#include <random>
#include "Mat.hpp"

namespace M
{

template <typename T>
Mat<T> generateRandomMatrix( Index rows, Index cols )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_real_distribution<> dis( -1.0, 1.0 );

    Mat<T> m( rows, cols );

    for ( Index i = 0; i < rows; ++i ) {
        for ( Index j = 0; j < cols; ++j ) {
            m( i, j ) = dis( gen );
        }
    }

    return m;
}

}