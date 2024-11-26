// S1 2024, MODEL
//
// Author: Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Utility functions for comparing matrices in unit tests.

#include <Mat.hpp>

namespace M
{

// Helper function to check if two matrices are approximately equal
template<typename T, typename Impl1, typename Impl2, MatDim n, MatDim m>
bool matrix_near(const MatFacade<Impl1, T, n, m>& a,
                 const MatFacade<Impl2, T, n, m>& b,
                 T epsilon = 1e-10) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;

    for (Index i = 0; i < a.rows(); ++i)
        for (Index j = 0; j < a.cols(); ++j)
            if (std::abs(a(i, j) - b(i, j)) > epsilon)
                return false;
    return true;
}

}
