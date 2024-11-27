// S1 2024, MODEL
//
// Author: Emilie Gillet (emilie.gillet@etu.sorbonne-universite.fr)
//
// -----------------------------------------------------------------------------
//
// Unit tests for Strassen Multiplication.

#include <io.hpp>
#include <StrassenMultiplication.hpp>

#include "test_utils.hpp"

using namespace M;
using namespace std;

namespace {

typedef std::tuple<Index, Index, Index> MatrixProductSize;

class TestStrassenMultiplication : 
    public ::testing::TestWithParam<MatrixProductSize> {
};

TEST_P(TestStrassenMultiplication, CompareWithNaiveMultiplication) {
    auto [m, n, p] = GetParam();
    auto a = generateRandomMatrix(m, n);
    auto b = generateRandomMatrix(n, p);
    auto c = multiplyStrassen(a, b);

    EXPECT_TRUE(isMatrixNear(c, a * b));
}

const MatrixProductSize kTestSizes[] = {
    { 3, 4, 5 },
    { 30, 40, 50 },
    { 100, 100, 100 },
    { 257, 65, 33 },
    { 65, 33, 257 },
    { 65, 129, 65 },
    { 129, 65, 129 },
};

INSTANTIATE_TEST_SUITE_P(
    TestStrassenMultiplicationAllSizes,
    TestStrassenMultiplication,
    testing::ValuesIn(kTestSizes));

TEST(BenchmarkStrassenMultiplication, SquareSize) {
    // 1.5s
    Index n = 2222;
    auto a = generateRandomMatrix(n, n);
    auto b = generateRandomMatrix(n, n);
    auto c = multiplyStrassen(a, b);
}

// The following test cases were used during development to validate the
// splitting of wide or tall matrices into squarish blocks in multiplySubmatrix
// and to set the constants kInnerSplitRatio and kOuterSplitRatio

TEST(BenchmarkStrassenMultiplicationSplitStrategy, WideA) {
    // 1.9s without split, 1.1s with split in multiplySubmatrix
    auto a = generateRandomMatrix(500, 4000);
    auto b = generateRandomMatrix(4000, 4000);
    auto c = multiplyStrassen(a, b);
}

TEST(BenchmarkStrassenMultiplicationSplitStrategy, TallB) {
    // 1.9s without split, 1.1s with split in multiplySubmatrix
    auto a = generateRandomMatrix(4000, 4000);
    auto b = generateRandomMatrix(4000, 500);
    auto c = multiplyStrassen(a, b);
}

TEST(BenchmarkStrassenMultiplicationSplitStrategy, TallAWideB) {
    // Split optimization has no impact - 1.5s in both cases.
    auto a = generateRandomMatrix(4000, 800);
    auto b = generateRandomMatrix(800, 4000);
    auto c = multiplyStrassen(a, b);
}

TEST(BenchmarkStrassenMultiplicationSplitStrategy, TallBWideA) {
    // 1.2s without split, 0.8s with split.
    auto a = generateRandomMatrix(800, 8000);
    auto b = generateRandomMatrix(8000, 800);
    auto c = multiplyStrassen(a, b);
}

}  // namespace
