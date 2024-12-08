// MU4IN901 - MODEL - Implementation project
//
// Authors:
// Andrei Aralov (andrei.aralov@etu.sorbonne-universite.fr)
// Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
// Ruchi  Thareja (ruchi.thareja@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Generation of matrices with random values.

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