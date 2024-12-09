// MU4IN901 - MODEL - Implementation project
//
// Authors:
// Andrei Aralov (andrei.aralov@etu.sorbonne-universite.fr)
// Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
// Ruchi  Thareja (ruchi.thareja@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Implementation of the Strassen algorithm for matrix multiplication.
//
// This implementation does not use any particular data structure or
// representation, and directly operates on Mat objects.

#pragma once

#include <Mat.hpp>
#include <MatView.hpp>

#include <algorithm>
#include <tuple>
#include <vector>
#include <cstring>

namespace M
{

typedef std::tuple<Index, Index, Index> MatrixProductSize;

// These functors describe how the matrices will be divided in 4 sub-blocks at
// each iteration of the Strassen algorithms: given a dimension of a block n,
// the functor returns the point k at which the "cut" will be made.
//
// |------ n -------|
// |--- k ---|
//
// +---------+------+
// |         |      |
// |   A_11  | A_12 |
// |         |      |
// |         |      |
// +---------+------+
// |         |      |
// |   A_21  | A_22 |
// |         |      |
// +---------+------+
//
//
//
// NaiveSplitPolicy: cut in half, ie k = ceil(n / 2).
//
// PowerOfTwoSplitPolicy: cut so that one of the pieces has a dimension that
// is a power of 2. Cut in half if the split is too unbalanced.

struct NaiveSplitPolicy {

Index operator()(Index n, Index min_size) {
    return (n + 1) / 2;
}

};

struct PowerOfTwoSplitPolicy {

Index operator()(Index n, Index min_size) {
    // The ideal split point is in the half.
    Index half = n / 2;

    // Search for the nearest min_size x 2^k close to the half.
    Index next_pow = min_size;
    while (next_pow < half) {
        next_pow <<= 1;
    }
    if ((next_pow - half) > (half - next_pow / 2)) {
        next_pow >>= 1;
    }

    // If the split is too uneven, use the middle point.
    Index split_point = abs(next_pow - half) <= 8 ? next_pow : half;
    return  std::max(split_point, n - split_point);
}

};

// Addition and multiplication routines are instantiated in two versions:
// - the OVERWRITE version which performs C <- op
// - the ACCUMULATE version which performs C <- C + op
enum OpMode {
    OVERWRITE,
    ACCUMULATE
};

#define EXPLICIT_VECTORIZATION

#ifdef EXPLICIT_VECTORIZATION

// Perform the operation C <- A x B or C <- C + A x B on a submatrix (view).
// This will be called when the Strassen iterations has divided the
// multiplication problem into small enough matrices.
template<
    OpMode op=OVERWRITE,
    typename SubmatrixC,
    typename SubmatrixA,
    typename SubmatrixB>
void multiplySubmatrixLeaf(
        SubmatrixC& c,
        const SubmatrixA& a,
        const SubmatrixB& b) {
    Index m = a.rows();
    Index n = a.cols();
    Index p = b.cols();

    assert(b.rows() == n);
    assert(c.rows() == m && c.cols() == p);

    using ElemT = typename SubmatrixC::ElemT;

    // Define a 256-bit SIMD type
    typedef ElemT block __attribute__ (( vector_size(32) ));

    // Number of elements that can fit in a SIMD register
    constexpr Index N = (256 / 8) / sizeof(ElemT);

    // Copy A and B (transposed) to SIMD registers
    Index num_blocks = (n + N - 1) / N;
    alignas(32) block a_blocks[m * num_blocks];
    alignas(32) block b_t_blocks[p * num_blocks];
    std::memset(a_blocks, 0, sizeof(a_blocks));
    std::memset(b_t_blocks, 0, sizeof(b_t_blocks));

    for (Index i = 0; i < m; ++i) {
        for (Index j = 0; j < n; ++j) {
            a_blocks[i * num_blocks + (j / N)][j % N] = a(i, j);
        }
    }

    for (Index i = 0; i < n; ++i) {
        for (Index j = 0; j < p; ++j) {
            b_t_blocks[j * num_blocks + (i / N)][i % N] = b(i, j);
        }
    }

    for (Index i = 0; i < m; ++i) {
        for (Index j = 0; j < p; ++j) {
            block s{};
            for (Index k = 0; k < num_blocks; ++k) {
                s += a_blocks[i * num_blocks + k] * b_t_blocks[j * num_blocks + k];
            }
            typename SubmatrixC::ElemT sum = 0;
            for (Index k = 0; k < N; ++k) {
                sum += s[k];
            }
            if (op == ACCUMULATE) {
                c(i, j) += sum;
            } else {
                c(i, j) = sum;
            }
        }
    }
}

#else

// Perform the operation C <- A x B or C <- C + A x B on a submatrix (view).
// This will be called when the Strassen iterations has divided the
// multiplication problem into small enough matrices.
template<
    OpMode op=OVERWRITE,
    typename SubmatrixC,
    typename SubmatrixA,
    typename SubmatrixB>
void multiplySubmatrixLeaf(
        SubmatrixC& c,
        const SubmatrixA& a,
        const SubmatrixB& b) {
    Index m = a.rows();
    Index n = a.cols();
    Index p = b.cols();

    assert(b.rows() == n);
    assert(c.rows() == m && c.cols() == p);

    // Copy on the stack the transpose of B to speed up
    // the subsequent multiplication.
    typename SubmatrixC::ElemT b_transpose[p * n];
    for (Index i = 0; i < n; ++i) {
        for (Index j = 0; j < p; ++j) {
            b_transpose[j * n + i] = b(i, j);
        }
    }

    for (Index i = 0; i < m; ++i) {
        for (Index j = 0; j < p; ++j) {
            typename SubmatrixC::ElemT s = 0;
            for (Index k = 0; k < n; ++k) {
                s += a(i, k) * b_transpose[j * n + k];
                // Uncomment this to double the computation time.
                // s += a(i, k) * b(k, j);
            }
            if (op == ACCUMULATE) {
                c(i, j) += s;
            } else {
                c(i, j) = s;
            }
        }
    }
}

#endif // EXPLICIT_VECTORIZATION

// Perform the operation C <- A ± B or C <- C + A ± B on a
// submatrix (view) without copy. Handle mismatches in sizes of C
// A and B by padding with zeros.
template<
    int b_sign=1,
    OpMode op=OVERWRITE,
    typename SubmatrixC,
    typename SubmatrixA,
    typename SubmatrixB>
void accumulateSubmatrix(
        SubmatrixC& c,
        const SubmatrixA& a,
        const SubmatrixB& b) {
    Index a_m = a.rows();
    Index a_n = a.cols();
    Index b_m = b.rows();
    Index b_n = b.cols();
    Index c_m = c.rows();
    Index c_n = c.cols();

    Index m = std::min(std::min(a_m, b_m), c_m);
    Index n = std::min(std::min(a_n, b_n), c_n);

    // Handle the common domain: no boundary checks.
    constexpr typename SubmatrixC::ElemT b_coefficient = b_sign;
    for (Index i = 0; i < m; ++i) {
        for (Index j = 0; j < n; ++j) {
            if (op == ACCUMULATE) {
                c(i, j) += a(i, j) + b(i, j) * b_coefficient;
            } else {
                c(i, j) = a(i, j) + b(i, j) * b_coefficient;
            }
        }
    }

    if (a_m == b_m && b_m == c_m && a_n == b_n && b_n == c_n) {
        return;
    }

    // Handle the fringe, with boundary checks and zero-padding.
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

// The Strassen algorithm requires many sub-matrices to be computed. Allocating
// space in memory for them accounted for about 30% of the compute time in our
// tests. Thus, all the memory needed is pre-allocated, and the resulting
// workspace is passed throughout iterations.
//
// For a given recursion level, the workspace consists of 6 matrices:
// - 1 matrix to store the differences of elements from A (size m x n)
// - 1 matrix to store the differences of elements from B (size n x p)
// - 4 matrices storing sub-products M1 to M7, with reuse (size m x p)
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
using StrassenWorkspace = std::vector<StrassenIterationWorkspace<T> >;

// Create the workspace required to multiply an m * n matrix by an n * p matrix
// using the Strassen algorithm and naive multiplication for sizes <= min_size.
template<typename T, typename SplitPolicy>
StrassenWorkspace<T> makeStrassenWorkspace(
    Index m, Index n, Index p, Index min_size) {
    StrassenWorkspace<T> workspace;
    SplitPolicy split;
    while (m > min_size && n > min_size && p > min_size) {
        m = split(m, min_size);
        n = split(n, min_size);
        p = split(p, min_size);
        workspace.push_back( {
            Mat<T>(m, n),  // a
            Mat<T>(n, p),  // b
            Mat<T>(m, p),  // m1
            Mat<T>(m, p),  // m4, m6
            Mat<T>(m, p),  // m5, m2
            Mat<T>(m, p)   // m7, m3
        } );
    }
    return workspace;
}



// Split a submatrix into 4 smaller blocks, such that the upper-left block is of
// size m2 x n2
template<typename T>
auto splitSubmatrix(T& a, Index m2, Index n2) {
    Index m = a.rows();
    Index n = a.cols();
    return std::tuple{
        a.submatrix(0, 0, m2, n2),
        a.submatrix(0, n2, m2, n - n2),
        a.submatrix(m2, 0, m - m2, n2),
        a.submatrix(m2, n2, m - m2, n - n2)
    };
}

// Multiply two submatrices using the Strassen algorithm.
template<
    OpMode op=OVERWRITE,
    Index min_size=16,
    typename SplitPolicy=NaiveSplitPolicy,
    typename T,
    MatDim m_a, MatDim n_a, MatDim n_b, MatDim p_b,
    typename ImplA, typename ImplB, typename ImplC>
void multiplySubmatrixSquare(
        MatFacade<ImplC, T, m_a, p_b>& c,
        const MatFacade<ImplA, T, m_a, n_a>& a,
        const MatFacade<ImplB, T, n_b, p_b>& b,
        typename StrassenWorkspace<T>::iterator workspace) {
    static_assert(n_a == n_b);

    Index m = a.rows();
    Index n = a.cols();
    Index p = b.cols();

    assert(b.rows() == n);
    assert(c.rows() == m && c.cols() == p);

    // Stop the recursion when one of the dimensions is small enough.
    if (m <= min_size || n <= min_size || p <= min_size) {
        if (m == min_size && n == min_size && p == min_size) {
            // Use a specialization of multiplySubmatrixLeaf when the
            // dimensions are exactly those of the ideal leaf size.
            auto a_ = a.template submatrix<min_size, min_size>(0, 0);
            auto b_ = b.template submatrix<min_size, min_size>(0, 0);
            auto c_ = c.template submatrix<min_size, min_size>(0, 0);
            multiplySubmatrixLeaf<op>(c_, a_, b_);
        } else {
            multiplySubmatrixLeaf<op>(c, a, b);
        }
        return;
    }

    // Use the provided split policy to decide where to split A and B.
    SplitPolicy split;
    Index m2 = split(m, min_size);
    Index n2 = split(n, min_size);
    Index p2 = split(p, min_size);

    // Split A, B and C.
    auto [a_11, a_12, a_21, a_22] = splitSubmatrix(a, m2, n2);
    auto [b_11, b_12, b_21, b_22] = splitSubmatrix(b, n2, p2);
    auto [c_11, c_12, c_21, c_22] = splitSubmatrix(c, m2, p2);

    // M_1 \leftarrow (A_{11} + A_{22}) \times (B_{11} + B_{22})
    auto m_1 = workspace->m_1.submatrix(0, 0, m2, p2);
    auto a_11_plus_a_22 = workspace->a.submatrix(0, 0, m2, n2);
    auto b_11_plus_b_22 = workspace->b.submatrix(0, 0, n2, p2);
    accumulateSubmatrix<+1>(a_11_plus_a_22, a_11, a_22);
    accumulateSubmatrix<+1>(b_11_plus_b_22, b_11, b_22);
    multiplySubmatrixSquare<OVERWRITE, min_size, SplitPolicy>(
        m_1, a_11_plus_a_22, b_11_plus_b_22, workspace + 1);

    // M_4 \leftarrow A_{22} \times (B_{21} - B_{11})
    auto m_4 = workspace->m_4_m_6.submatrix(0, 0, m - m2, p2);
    auto b_21_minus_b_11 = workspace->b.submatrix(0, 0, n - n2, p2);
    accumulateSubmatrix<-1>(b_21_minus_b_11, b_21, b_11);
    multiplySubmatrixSquare<OVERWRITE, min_size, SplitPolicy>(
        m_4, a_22, b_21_minus_b_11, workspace + 1);

    // M_5 \leftarrow (A_{11} + A_{12}) \times B_{22}
    auto m_5 = workspace->m_5_m_2.submatrix(0, 0, m2, p - p2);
    auto a_11_plus_a_12 = workspace->a.submatrix(0, 0, m2, n - n2);
    accumulateSubmatrix<+1>(a_11_plus_a_12, a_11, a_12);
    multiplySubmatrixSquare<OVERWRITE, min_size, SplitPolicy>(
        m_5, a_11_plus_a_12, b_22, workspace + 1);

    // M_7 \leftarrow (A_{12} - A_{22}) \times (B_{21} + B_{22})
    auto m_7 = workspace->m_7_m_3.submatrix(0, 0, m2, p2);
    auto a_12_minus_a_22 = workspace->a.submatrix(0, 0, m2, n - n2);
    auto b_21_plus_b_22 = workspace->b.submatrix(0, 0, n - n2, p2);

    accumulateSubmatrix<-1>(a_12_minus_a_22, a_12, a_22);
    accumulateSubmatrix<+1>(b_21_plus_b_22, b_21, b_22);
    multiplySubmatrixSquare<OVERWRITE, min_size, SplitPolicy>(
        m_7, a_12_minus_a_22, b_21_plus_b_22, workspace + 1);

    // C_{11} \leftarrow M_1 + M_4 - M_5 + M_7
    accumulateSubmatrix<+1, op>(c_11, m_1, m_4);
    accumulateSubmatrix<-1, ACCUMULATE>(c_11, m_7, m_5);

    // M_3 \leftarrow A_{11} \times (B_{12} - B_{22})
    // M_7 is no longer needed, buffer is reused for M3
    auto m_3 = workspace->m_7_m_3.submatrix(0, 0, m2, p - p2);
    auto b_12_minus_b_22 = workspace->b.submatrix(0, 0, n2, p - p2);
    accumulateSubmatrix<-1>(b_12_minus_b_22, b_12, b_22);
    multiplySubmatrixSquare<OVERWRITE, min_size, SplitPolicy>(
        m_3, a_11, b_12_minus_b_22, workspace + 1);

    // C_{12} \leftarrow M_3 + M_5
    accumulateSubmatrix<+1, op>(c_12, m_3, m_5);

    // M_2 \leftarrow (A_{21} + A_{22}) \times B_{11}
    // M_5 is no longer needed, buffer is reused;
    auto m_2 = workspace->m_5_m_2.submatrix(0, 0, m - m2, p2);
    auto a_21_plus_a_22 = workspace->a.submatrix(0, 0, m - m2, n2);
    accumulateSubmatrix<+1>(a_21_plus_a_22, a_21, a_22);
    multiplySubmatrixSquare<OVERWRITE, min_size, SplitPolicy>(
        m_2, a_21_plus_a_22, b_11, workspace + 1);

    // C_{21} \leftarrow M_2 + M_4
    accumulateSubmatrix<1, op>(c_21, m_2, m_4);

    // M_6 \leftarrow (A_{21} - A_{11}) \times (B_{11} + B_{12})
    // M_4 is no longer needed, buffer is reused
    auto m_6 = workspace->m_4_m_6.submatrix(0, 0, m - m2, p - p2);
    auto a_21_minus_a_11 = workspace->a.submatrix(0, 0, m - m2, n2);
    auto b_11_plus_b_12 = workspace->b.submatrix(0, 0, n2, p - p2);
    accumulateSubmatrix<-1>(a_21_minus_a_11, a_21, a_11);
    accumulateSubmatrix<+1>(b_11_plus_b_12, b_11, b_12);
    multiplySubmatrixSquare<OVERWRITE, min_size, SplitPolicy>(
        m_6, a_21_minus_a_11, b_11_plus_b_12, workspace + 1);

    // C_{22} \leftarrow M_1 - M_2 + M_3 + M_6
    accumulateSubmatrix<-1, op>(c_22, m_1, m_2);
    accumulateSubmatrix<+1, ACCUMULATE>(c_22, m_3, m_6);
}

// Performs the operation C <- A * B or C <- C + A * B,
// A, B and C being submatrices.
// This top-level function breaks down the product into several products of
// smaller squarish matrices for which the Strassen algorithm is used.
template <
    OpMode op=OVERWRITE,
    Index min_size=64,
    typename SplitPolicy=PowerOfTwoSplitPolicy,
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

    const Index kOuterSplitRatio = 8;
    const Index kInnerSplitRatio = 2;
    const Index kMaxSizeWithoutSplit = 4096;

    if ((m <= min_size || n <= min_size || p <= min_size) && \
        (m <= kMaxSizeWithoutSplit) && \
        (n <= kMaxSizeWithoutSplit) && \
        (p <= kMaxSizeWithoutSplit)) {
        // A and B are small or skinny: directly use the
        // explicit multiplication.
        multiplySubmatrixLeaf<op>(c, a, b);
    } else if (m >= kOuterSplitRatio * n) {
        // A is too tall: break down into smaller blocks.
        Index num_outer_blocks = (m - 1) / n + 1;
        for (Index i = 0; i < num_outer_blocks; ++i) {
            Index h = std::min(n, m - i * n);
            auto sub_c = c.submatrix(i * n, 0, h, p);
            auto sub_a = a.submatrix(i * n, 0, h, n);
            multiplySubmatrix<op, min_size, SplitPolicy>(sub_c, sub_a, b);
        }
    } else if (n >= kInnerSplitRatio * m) {
        // A is too wide: break down into smaller blocks.
        Index num_inner_blocks = (n - 1) / m + 1;
        for (Index i = 0; i < num_inner_blocks; ++i) {
            Index w = std::min(m, n - i * m);
            auto sub_a = a.submatrix(0, i * m, m, w);
            auto sub_b = b.submatrix(i * m, 0, w, p);
            multiplySubmatrix<ACCUMULATE, min_size, SplitPolicy>(c, sub_a, sub_b);
        }
    } else if (n >= kInnerSplitRatio * p) {
        // B is too tall: break down into smaller blocks.
        Index num_inner_blocks = (n - 1) / p + 1;
        for (Index i = 0; i < num_inner_blocks; ++i) {
            Index h = std::min(p, n - i * p);
            auto sub_a = a.submatrix(0, i * p, m, h);
            auto sub_b = b.submatrix(i * p, 0, h, p);
            multiplySubmatrix<ACCUMULATE, min_size, SplitPolicy>(c, sub_a, sub_b);
        }
    } else if (p >= kOuterSplitRatio * m) {
        // B is too wide: break down into smaller blocks.
        Index num_outer_blocks = (p - 1) / n + 1;
        for (Index i = 0; i  < num_outer_blocks; ++i) {
            Index w = std::min(n, p - i * n);
            auto sub_c = c.submatrix(0, i * n, m, w);
            auto sub_b = b.submatrix(0, i * n, n, w);
            multiplySubmatrix<op, min_size, SplitPolicy>(sub_c, a, sub_b);
        }
    } else {
        // A and B have a squarish size: the Strassen algorithm can be used.
        auto workspace = makeStrassenWorkspace<T, SplitPolicy>(m, n, p, min_size);
        multiplySubmatrixSquare<op, min_size, SplitPolicy>(c, a, b, workspace.begin());
    }
}

// Top-level function to implement Strassen multiplication on two matrix-like
// objects.
template <
    Index min_size=64,
    typename SplitPolicy=PowerOfTwoSplitPolicy,
    typename T,
    typename ImplA, MatDim m, MatDim n_a,
    typename ImplB, MatDim n_b, MatDim p>
Mat<T, m, p> multiplyStrassen(
        const MatFacade<ImplA, T, m, n_a>& a,
        const MatFacade<ImplB, T, n_b, p>& b) {
    static_assert(n_a == n_b);
    Mat<T, m, p> c(a.rows(), b.cols());
    multiplySubmatrix<OVERWRITE, min_size, SplitPolicy>(c, a, b);
    return c;
}

// For comparison purposes, a non-Strassen, but efficient, multiplication
// algorithm. It breaks down the product into macro_size x macro_size products,
// then each of these into micro_size x micro_size products, and use the
// optimized routine multiplySubmatrixLeaf on the tiles.
template<
    OpMode op=OVERWRITE,
    Index micro_size=32,
    Index macro_size=256,
    typename T,
    MatDim m_a, MatDim n_a, MatDim n_b, MatDim p_b,
    typename ImplA, typename ImplB, typename ImplC>
void multiplyTiled(
        MatFacade<ImplC, T, m_a, p_b>& c,
        const MatFacade<ImplA, T, m_a, n_a>& a,
        const MatFacade<ImplB, T, n_b, p_b>& b) {
    static_assert(n_a == n_b);

    Index m = a.rows();
    Index n = a.cols();
    Index p = b.cols();

    assert(b.rows() == n);
    assert(c.rows() == m && c.cols() == p);

    if (op == OVERWRITE) {
        for (Index i = 0; i < m; ++i) {
            for (Index j = 0; j < p; ++j) {
                c(i, j) = 0;
            }
        }
    }

    for (Index j_macro = 0; j_macro < p; j_macro += macro_size) {
        for (Index k_macro = 0; k_macro < n; k_macro += macro_size) {
            for (Index i_macro = 0; i_macro < m; i_macro += macro_size) {
                for (Index j = j_macro; j < std::min(j_macro + macro_size, p); j += micro_size) {
                    for (Index k = k_macro; k < std::min(k_macro + macro_size, n); k += micro_size) {
                        for (Index i = i_macro; i < std::min(i_macro + macro_size, m); i += micro_size) {
                            Index mm = std::min(micro_size, m - i);
                            Index nn = std::min(micro_size, n - k);
                            Index pp = std::min(micro_size, p - j);
                            auto c_ij = c.submatrix(i, j, mm, pp);
                            auto a_ik = a.submatrix(i, k, mm, nn);
                            auto b_kj = b.submatrix(k, j, nn, pp);
                            if (mm == micro_size && nn == micro_size && pp == micro_size) {
                                // Force the instantiation of a
                                // specialized kernel when the size is
                                // exactly micro_size x micro_size x micro_size
                                auto a_ = a_ik.template submatrix<micro_size, micro_size>(0, 0);
                                auto b_ = b_kj.template submatrix<micro_size, micro_size>(0, 0);
                                auto c_ = c_ij.template submatrix<micro_size, micro_size>(0, 0);
                                multiplySubmatrixLeaf<ACCUMULATE>(c_, a_, b_);
                            } else {
                                multiplySubmatrixLeaf<ACCUMULATE>(c_ij, a_ik, b_kj);
                            }
                        }
                    }
                }
            }
        }
    }
}

template<
    Index micro_size=32,
    Index macro_size=256,
    typename T,
    typename ImplA, MatDim m, MatDim n_a,
    typename ImplB, MatDim n_b, MatDim p>
Mat<T, m, p> multiplyTiled(
        const MatFacade<ImplA, T, m, n_a>& a,
        const MatFacade<ImplB, T, n_b, p>& b) {
    static_assert(n_a == n_b);
    Mat<T, m, p> c(a.rows(), b.cols());
    multiplyTiled<ACCUMULATE, micro_size, macro_size>(c, a, b);
    return c;
}

}