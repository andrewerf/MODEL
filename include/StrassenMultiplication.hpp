// S1 2024, MODEL
//
// Author: Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Implementation of the Strassen algorithm for matrix multiplication.
//
// This implementation does not use any particular data structure or
// representation, and directly operates on Mat objects.

#pragma once

#include <algorithm>

#include <Mat.hpp>
#include <MatView.hpp>

namespace M
{

using namespace std;

// Addition and multiplication routines are instantiated in two variants:
// - the OVERWRITE variant which performs C <- op
// - the ACCUMULATE variant which performs C <- C + op
enum OpMode {
    OVERWRITE,
    ACCUMULATE
};

// Perform the operation C <- A x B or C <- C + A x B on a submatrix (view).
// This will be called when the Strassen iterations has divided the
// multiplication into small enough matrices.
template<OpMode op=OVERWRITE, typename ViewC, typename ViewA, typename ViewB>
void multiplySubmatrixLeaf(ViewC& c, const ViewA& a, const ViewB& b) {
    Index m = a.rows();
    Index n = a.cols();
    Index p = b.cols();

    assert(b.rows() == n);
    assert(c.rows() == m && c.cols() == p);
    
    // Copy on the stack the transpose of B to speed up
    // the subsequent multiplication.
    typename ViewB::ElemT b_transpose[p * n];
    for (Index i = 0; i < n; ++i) {
        for (Index j = 0; j < p; ++j) {
            b_transpose[j * n + i] = b(i, j);
        }
    }
    
    for (Index i = 0; i < m; ++i) {
        for (Index j = 0; j < p; ++j) {
            typename ViewC::ElemT s = 0;
            for (Index k = 0; k < n; ++k) {
                s += a(i, k) * b_transpose[j * n + k];
                // Try to uncomment this for an awful slow-down.
                //s += a(i, k) * b(k, j);
            }
            if (op == ACCUMULATE) {
                c(i, j) += s;
            } else {
                c(i, j) = s;
            }
        }
    }
}

// Perform the operation C <- A + B * coefficient or
// C <- C + A + B * coefficient on a submatrix (view) without copy.
// Handle mismatches in sizes of C, A and B by padding with zeros.
template<
    OpMode op=OVERWRITE,
    typename ViewC,
    typename ViewA,
    typename ViewB>
void addSubmatrix(
        ViewC& c,
        const ViewA& a,
        const ViewB& b,
        typename ViewC::ElemT b_coefficient) {
    Index a_m = a.rows();
    Index a_n = a.cols();
    Index b_m = b.rows();
    Index b_n = b.cols();
    Index c_m = c.rows();
    Index c_n = c.cols();
    
    Index m = min(min(a_m, b_m), c_m);
    Index n = min(min(a_n, b_n), c_n);

    // Handle the common domain: no boundary checks.
    for (Index i = 0; i < m; ++i) {
        for (Index j = 0; j < n; ++j) {
            if (op == ACCUMULATE) {
                c(i, j) += a(i, j) + b(i, j) * b_coefficient;
            } else {
                c(i, j) = a(i, j) + b(i, j) * b_coefficient;
            }
        }
    }
    
    // Handle the edges, with boundary checks and zero-padding.
    for (Index i = 0; i < c_m; ++i) {
        Index start = i < m ? n : 0;
        for (Index j = start; j < c_n; ++j) {
            auto a_v = (i < a_m && j < a_n) ? a(i, j) : 0;
            auto b_v = (i < b_m && j < b_n) ? b(i, j) : 0;
            if (op == ACCUMULATE) {
                c(i, j) += a_v + b_v * b_coefficient;
            } else {
                c(i, j) = a_v + b_v * b_coefficient;
            }
        }
    }
}

// The pre-allocated workspace passed to the Strassen recursion.
// For a given recursion level, the workspace consists of 6 matrices:
// - one matrix to store the differences of elements from A.
// - one matrix to store the differences of elements from B.
// - four matrices storing sub-products M1 to M7 (with reuse).
template<typename T>
struct StrassenIterationWorkspace {
    Mat<T> a;
    Mat<T> b;
    Mat<T> m_1;
    Mat<T> m_4_m_6;
    Mat<T> m_5_m_2;
    Mat<T> m_7_m_3;
};

template<typename T>
using StrassenWorkspace = vector<StrassenIterationWorkspace< T> >;

// Create the workspace required to multiply an m * n and n * p matrix
// using the Strassen algorithm and regular multiplication for
// sizes <= min_size).
template<typename T>
StrassenWorkspace<T> makeStrassenWorkspace(
    Index m, Index n, Index p, Index min_size) {
    StrassenWorkspace<T> workspace;
    while (m > min_size && n > min_size && p > min_size) {
        m = (m + 1) / 2;
        n = (n + 1) / 2;
        p = (p + 1) / 2;
        workspace.push_back( {
            Mat<T>(m, n),  // a
            Mat<T>(n, p),  // b
            Mat<T>(m, p),  // m
            Mat<T>(m, p),  // m
            Mat<T>(m, p),  // m
            Mat<T>(m, p)   // m
        } );
    }
    return workspace;
}

// Split a submatrix into 4 smaller blocks of approximately equal size.
template<typename T>
auto splitSubmatrix(T& a) {
    Index m = a.rows();
    Index n = a.cols();
    Index m2 = (m + 1) / 2;
    Index n2 = (n + 1) / 2;
    return std::tuple{
        a.submatrix(0, 0, m2, n2),
        a.submatrix(0, n2, m2, n - n2),
        a.submatrix(m2, 0, m - m2, n2),
        a.submatrix(m2, n2, m - m2, n - n2)
    };
}

// Multiply two submatrices using the Strassen algorithm, which performs
// optimally when a and b are "squarish" matrices
// (0.5 < ratio of dimensions < 2).
template<
    OpMode op=OVERWRITE,
    Index min_size=16,
    typename T,
    MatDim mm, MatDim nn, MatDim pp,
    typename ImplA, typename ImplB, typename ImplC>
void multiplySubmatrixSquare(
        MatFacade<ImplC, T, mm, pp>& c,
        const MatFacade<ImplA, T, mm, nn>& a,
        const MatFacade<ImplB, T, nn, pp>& b,
        typename StrassenWorkspace<T>::iterator workspace) {
    Index m = a.rows();
    Index n = a.cols();
    Index p = b.cols();
    assert(b.rows() == n);
    assert(c.rows() == m && c.cols() == p);

    // No need to subdivide when one of the dimensions is too small.
    if (m <= min_size || n <= min_size || p <= min_size) {
        multiplySubmatrixLeaf<op>(c, a, b);
        return;
    }
    
    Index m2 = (m + 1) / 2;
    Index n2 = (n + 1) / 2;
    Index p2 = (p + 1) / 2;
    
    auto [a_11, a_12, a_21, a_22] = splitSubmatrix(a);
    auto [b_11, b_12, b_21, b_22] = splitSubmatrix(b);
    auto [c_11, c_12, c_21, c_22] = splitSubmatrix(c);
    
    // M_1 \leftarrow (A_{11} + A_{22}) \times (B_{11} + B_{22})
    auto m_1 = workspace->m_1.submatrix(0, 0, m2, p2);
    auto a_11_plus_a_22 = workspace->a.submatrix(0, 0, m2, n2);
    auto b_11_plus_b_22 = workspace->b.submatrix(0, 0, n2, p2);
    addSubmatrix(a_11_plus_a_22, a_11, a_22, 1);
    addSubmatrix(b_11_plus_b_22, b_11, b_22, 1);
    multiplySubmatrixSquare<OVERWRITE, min_size>(
        m_1, a_11_plus_a_22, b_11_plus_b_22, workspace + 1);
    
    // M_4 \leftarrow A_{22} \times (B_{21} - B_{11})
    auto m_4 = workspace->m_4_m_6.submatrix(0, 0, m - m2, p2);
    auto b_21_minus_b_11 = workspace->b.submatrix(0, 0, n - n2, p2);
    addSubmatrix(b_21_minus_b_11, b_21, b_11, -1);
    multiplySubmatrixSquare<OVERWRITE, min_size>(
        m_4, a_22, b_21_minus_b_11, workspace + 1);
    
    // M_5 \leftarrow (A_{11} + A_{12}) \times B_{22}
    auto m_5 = workspace->m_5_m_2.submatrix(0, 0, m2, p - p2);
    auto a_11_plus_a_12 = workspace->a.submatrix(0, 0, m2, n - n2);
    addSubmatrix(a_11_plus_a_12, a_11, a_12, 1);
    multiplySubmatrixSquare<OVERWRITE, min_size>(
        m_5, a_11_plus_a_12, b_22, workspace + 1);
    
    // M_7 \leftarrow (A_{12} - A_{22}) \times (B_{21} + B_{22})
    auto m_7 = workspace->m_7_m_3.submatrix(0, 0, m2, p2);
    auto a_12_minus_a_22 = workspace->a.submatrix(0, 0, m2, n - n2);
    auto b_21_plus_b_22 = workspace->b.submatrix(0, 0, n - n2, p2);
    addSubmatrix(a_12_minus_a_22, a_12, a_22, -1);
    addSubmatrix(b_21_plus_b_22, b_21, b_22, +1);
    multiplySubmatrixSquare<OVERWRITE, min_size>(
        m_7, a_12_minus_a_22, b_21_plus_b_22, workspace + 1);
    
    // C_{11} \leftarrow M_1 + M_4 - M_5 + M_7
    addSubmatrix<op>(c_11, m_1, m_4, 1);
    addSubmatrix<ACCUMULATE>(c_11, m_7, m_5, -1);
    
    // M_3 \leftarrow A_{11} \times (B_{12} - B_{22})
    // M_7 is no longer needed, buffer is reused for M3
    auto m_3 = workspace->m_7_m_3.submatrix(0, 0, m2, p - p2);
    auto b_12_minus_b_22 = workspace->b.submatrix(0, 0, n2, p - p2);
    addSubmatrix(b_12_minus_b_22, b_12, b_22, -1);
    multiplySubmatrixSquare<OVERWRITE, min_size>(
        m_3, a_11, b_12_minus_b_22, workspace + 1);
    
    // C_{12} \leftarrow M_3 + M_5
    addSubmatrix<op>(c_12, m_3, m_5, 1);
    
    // M_2 \leftarrow (A_{21} + A_{22}) \times B_{11}
    // M_5 is no longer needed, buffer is reused;
    auto m_2 = workspace->m_5_m_2.submatrix(0, 0, m - m2, p2);
    auto a_21_plus_a_22 = workspace->a.submatrix(0, 0, m - m2, n2);
    addSubmatrix(a_21_plus_a_22, a_21, a_22, 1);
    multiplySubmatrixSquare<OVERWRITE, min_size>(
        m_2, a_21_plus_a_22, b_11, workspace + 1);
    
    // C_{21} \leftarrow M_2 + M_4
    addSubmatrix<op>(c_21, m_2, m_4, 1);
    
    // M_6 \leftarrow (A_{21} - A_{11}) \times (B_{11} + B_{12})
    // M_4 is no longer needed, buffer is reused
    auto m_6 = workspace->m_4_m_6.submatrix(0, 0, m - m2, p - p2);
    auto a_21_minus_a_11 = workspace->a.submatrix(0, 0, m - m2, n2);
    auto b_11_plus_b_12 = workspace->b.submatrix(0, 0, n2, p - p2);
    addSubmatrix(a_21_minus_a_11, a_21, a_11, -1);
    addSubmatrix(b_11_plus_b_12, b_11, b_12, 1);
    multiplySubmatrixSquare<OVERWRITE, min_size>(
        m_6, a_21_minus_a_11, b_11_plus_b_12, workspace + 1);
    
    // C_{22} \leftarrow M_1 - M_2 + M_3 + M_6
    addSubmatrix<op>(c_22, m_1, m_2, -1);
    addSubmatrix<ACCUMULATE>(c_22, m_3, m_6, 1);
}

// Performs the operation c <- a * b, a, b and c being submatrices.
// This top-level function breaks down the product into several products of
// smaller squarish matrices for which the use of the Strassen algorithm is
// efficient.
template <
    OpMode op=OVERWRITE,
    Index min_size=16,
    typename T,
    MatDim m_a, MatDim n_a, MatDim n_b, MatDim p_b,
    typename ImplA, typename ImplB, typename ImplC>
void multiplySubmatrix(
        MatFacade<ImplC, T, m_a, p_b>& c,
        const MatFacade<ImplA, T, m_a, n_a>& a,
        const MatFacade<ImplB, T, n_b, p_b>& b) {
    static_assert(n_a == n_b);
    Index m = a.rows();
    Index n = a.cols();
    Index p = b.cols();
    assert(b.rows() == n);
    assert(c.rows() == m && c.cols() == p);
    
    const Index kMaxDimRatio = 2;
    
    if (m <= min_size || n <= min_size || p <= min_size) {
        // A and B are skinny: directly use the explicit multiplication.
        multiplySubmatrixLeaf<op>(c, a, b);
    } else if (m >= kMaxDimRatio * n) {
        // A is too tall: break down into smaller blocks.
        Index num_outer_blocks = (m - 1) / n + 1;
        for (Index i = 0; i < num_outer_blocks; ++i) {
            Index h = min(n, m - i * n);
            auto sub_c = c.submatrix(i * n, 0, h, p);
            auto sub_a = a.submatrix(i * n, 0, h, n);
            multiplySubmatrix<op, min_size>(sub_c, sub_a, b);
        }
    } else if (n >= kMaxDimRatio * m) {
        // A is too wide: break down into smaller blocks.
        Index num_inner_blocks = (n - 1) / m + 1;
        for (Index i = 0; i < num_inner_blocks; ++i) {
            Index w = min(m, n - i * m);
            auto sub_a = a.submatrix(0, i * m, m, w);
            auto sub_b = b.submatrix(i * m, 0, w, p);
            multiplySubmatrix<ACCUMULATE, min_size>(c, sub_a, sub_b);
        }
    } else if (n >= kMaxDimRatio * p) {
        // B is too tall: break down into smaller blocks.
        Index num_inner_blocks = (n - 1) / p + 1;
        for (Index i = 0; i < num_inner_blocks; ++i) {
            Index h = min(p, n - i * p);
            auto sub_a = a.submatrix(0, i * p, m, h);
            auto sub_b = b.submatrix(i * p, 0, h, p);
            multiplySubmatrix<ACCUMULATE, min_size>(c, sub_a, sub_b);
        }
    } else if (p >= kMaxDimRatio * m) {
        // B is too wide: break down into smaller blocks.
        Index num_outer_blocks = (p - 1) / n + 1;
        for (Index i = 0; i  < num_outer_blocks; ++i) {
            Index w = min(n, p - i * n);
            auto sub_c = c.submatrix(0, i * n, m, w);
            auto sub_b = b.submatrix(0, i * n, n, w);
            multiplySubmatrix<op, min_size>(sub_c, a, sub_b);
        }
    } else {
        // A and B have a squarish size.
        auto workspace = makeStrassenWorkspace<T>(m, n, p, min_size);
        multiplySubmatrixSquare<op, min_size>(c, a, b, workspace.begin());
    }
}

template <
    typename T,
    typename ImplA, MatDim m, MatDim n_a,
    typename ImplB, MatDim n_b, MatDim p,
    Index min_size=16>
Mat<T, m, p> multiplyStrassen(
        const MatFacade<ImplA, T, m, n_a>& a,
        const MatFacade<ImplB, T, n_b, p>& b) {
    static_assert( n_a == n_b );
    Mat<T, m, p> c( a.rows(), b.cols() );
    multiplySubmatrix<OVERWRITE, min_size>(c, a, b);
    return c;
}

}