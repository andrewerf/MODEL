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

#include <StrassenMultiplication.hpp>

namespace M
{

template<typename T>
struct StrassenInversionIterationWorkspace {
    Mat<T> ce;
    Mat<T> eb;
    Mat<T> Z_ebz;
};

template<typename T>
using StrassenInversionWorkspace = std::vector<
    StrassenInversionIterationWorkspace<T> >;

// Create the workspace required to invert an n * n matrix
// using the Strassen algorithm.
template<typename T>
StrassenInversionWorkspace<T> makeStrassenInversionWorkspace(Index n) {
    StrassenInversionWorkspace<T> workspace;
    while (n > 1) {
        Index n2 = (n + 1) / 2;
        workspace.push_back( {
            Mat<T>(n - n2, n2),  // ce
            Mat<T>(n2, n - n2),  // eb
            Mat<T>(n2, n2),  // Z and ebz
        } );
        n = n2;
    }
    return workspace;
}

template <
    typename T,
    typename ImplM, MatDim n_m, MatDim m_m, 
    typename ImplMInv, MatDim n_m_inv, MatDim m_m_inv>
bool invertSubmatrix(
        MatFacade<ImplM, T, n_m_inv, m_m_inv>& m_inv,
        const MatFacade<ImplMInv, T, n_m, m_m>& m,
        typename StrassenInversionWorkspace<T>::iterator workspace) {
    static_assert( n_m == m_m && n_m_inv == m_m_inv && n_m == m_m_inv);
    Index n = m.rows();
    assert( m.cols() == n );
    assert( m_inv.cols() == n && m_inv.rows() == n );

    if (n == 1) {
        if (m(0, 0) == 0.0) {
            return false;
        } else {
            m_inv(0, 0) = 1.0 / m(0, 0);
            return true;
        }
    } else if (n == 2) {
        auto a = m(0, 0), b = m(0, 1), c = m(1, 0), d = m(1, 1);
        auto det = a * d - b * c;
        if (det == 0) {
            return false;
        }
        det = 1 / det;
        m_inv(0, 0) = d * det;
        m_inv(0, 1) = -b * det;
        m_inv(1, 0) = -c * det;
        m_inv(1, 1) = a * det;
        return true;
    }
    
    Index n2 = (n + 1) / 2;
    
    // Split into 4 blocks.
    const auto [a, b, c, d] = splitSubmatrix(m, n2, n2);
    auto [e, y, z, t] = splitSubmatrix(m_inv, n2, n2);

    // Use pre-allocated space for blocks.
    auto ce = workspace->ce.submatrix(0, 0, n - n2, n2);
    auto eb = workspace->eb.submatrix(0, 0, n2, n - n2);
    auto Z = workspace->Z_ebz.submatrix(0, 0, n - n2, n - n2);
    auto ebz = workspace->Z_ebz.submatrix(0, 0, n2, n2);
    
    // e = a^{-1}
    if (!invertSubmatrix(e, a, workspace + 1)) {
        return false;
    }

    // HACK!!!
    #define multiply multiplySubmatrix<OVERWRITE, 64>
    //#define multiply multiplyTiled
    
    // Z = d - ceb
    multiply(ce, c, e);
    multiply(Z, ce, b);
    accumulateSubmatrix<-1>(Z, d, Z);
    
    // t = Z^{-1}
    if (!invertSubmatrix(t, Z, workspace + 1)) {
        return false;
    }
    
    // y = -ebt
    multiply(eb, e, b);
    multiply(y, eb, t);
    y = -y;
    
    // z = -tce
    multiply(z, t, ce);
    z = -z;

    // x = e + ebtce = e - ebz
    multiply(ebz, eb, z);
    accumulateSubmatrix<-1>(e, e, ebz);
    return true;
}

template <
    typename T,
    typename Impl, MatDim n_m, MatDim m_m>
std::optional<Mat<T, n_m, m_m> > inverseStrassenFast(
        const MatFacade<Impl, T, n_m, m_m>& m) {
    static_assert( n_m == m_m );
    assert( m.rows() == m.cols() );
    const Index n = m.rows();
    auto workspace = makeStrassenInversionWorkspace<T>(n);
    Mat<T> m_inv(n, n);
    if (!invertSubmatrix(m_inv, m, workspace.begin())) {
        return {};
    } else {
        return m_inv;
    }
}

}