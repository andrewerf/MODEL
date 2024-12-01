//
// Created by Andrey Aralov on 12/1/24.
//
#pragma once
#include <random>
#include "Mat.hpp"

namespace M
{

// Fill a random matrix with uniformly distributed values in the range [-1, 1]
template<typename T, typename Impl, MatDim m, MatDim n>
void fillWithRandomValues(MatFacade<Impl, T, m, n>& a) {
    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_real_distribution<> distribution(-1.0, 1.0);
    for ( Index i = 0; i < a.rows(); ++i ) {
        for ( Index j = 0; j < a.cols(); ++j ) {
            a( i, j ) = distribution(generator);
        }
    }
}

// Create a random matrix of size (m, n).
template<typename T=double>
Mat<T> generateRandomMatrix(Index m, Index n) {
    Mat<T> a(m, n);
    fillWithRandomValues(a);
    return a;
}

}