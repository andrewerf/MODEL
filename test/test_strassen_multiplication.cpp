// S1 2024, MODEL
//
// Author: Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Unit tests for Strassen Multiplication.

// #include <BLASMultiplication.hpp>
#include <random.hpp>
#include <StrassenMultiplication.hpp>

// #include "MiniBenchmark.hpp"
#include "test_utils.hpp"

using namespace M;
using namespace std;

namespace {

class TestStrassenMultiplication :
    public ::testing::TestWithParam<MatrixProductSize> {
};

class BenchmarkStrassenMultiplication :
    public ::testing::TestWithParam<MatrixProductSize> {
};

TEST_P(TestStrassenMultiplication, CompareWithNaiveMultiplication) {
    auto [m, n, p] = GetParam();
    auto a = generateRandomMatrix(m, n);
    auto b = generateRandomMatrix(n, p);
    auto c = multiplyStrassen(a, b);
    EXPECT_TRUE(isMatrixNear(c, a * b));
}

TEST_P(BenchmarkStrassenMultiplication, Benchmark) {
    auto [m, n, p] = GetParam();
    auto a = generateRandomMatrix(m, n);
    auto b = generateRandomMatrix(n, p);
    auto c = multiplyStrassen(a, b);
}

const MatrixProductSize kTestSizes[] = {
    { 3, 4, 5 },
    { 30, 40, 50 },
    { 100, 100, 100 },
    { 257, 65, 33 },
    { 65, 33, 257 },
    { 64, 129, 65 },
    { 129, 65, 129 },
    { 129, 222, 128 },
};

const MatrixProductSize kBenchmarkSizes[] = {
    // 1.2s [BLAS: 68ms]
    { 2222, 2222, 2222 },
    // 3.9s [BLAS: 184ms]
    { 2222, 3333, 4444 },

    // The following test cases were used during development to validate the
    // splitting of wide or tall matrices into squarish blocks
    // in multiplySubmatrix and to set the constants kInnerSplitRatio
    // and kOuterSplitRatio

    // 1.1s (1.9s without split) [BLAS: 45ms]
    { 500, 4000, 4000 }, 
    // 1.1s (1.9s without split) [BLAS: 44ms]
    { 4000, 4000, 500 },
    // 1.3s [BLAS: 71ms]
    { 4000, 800, 4000 },
    // 0.8s (1.3s without split) [BLAS: 27ms]
    { 800, 8000, 800 },
};

INSTANTIATE_TEST_SUITE_P(
    TestStrassenMultiplicationAllSizes,
    TestStrassenMultiplication,
    testing::ValuesIn(kTestSizes),
    ProductSizeToString());

INSTANTIATE_TEST_SUITE_P(
    BenchmarkStrassenMultiplicationAllSizes,
    BenchmarkStrassenMultiplication,
    testing::ValuesIn(kBenchmarkSizes),
    ProductSizeToString());

}  // namespace
