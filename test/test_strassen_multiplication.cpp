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

TEST_P(TestStrassenMultiplication, CompareWithNaiveProduct) {
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
    { 256, 64, 32 },
    { 64, 32, 256 },
    { 64, 128, 64 },
    { 128, 64, 128 },
};

INSTANTIATE_TEST_SUITE_P(
    TestStrassenMultiplicationAllSizes,
    TestStrassenMultiplication,
    testing::ValuesIn(kTestSizes));

TEST(BenchMarkStrassenMultiplication, FixedSize) {
    Index n = 2300;
    auto a = generateRandomMatrix(n, n);
    auto b = generateRandomMatrix(n, n);
    auto c = multiplyStrassen(a, b);
}

/*TEST(BenchMarkStrassenMultiplication, SkinnyMatrix) {
    auto a = generateRandomMatrix(600, 10000);
    auto b = generateRandomMatrix(10000, 600);
    auto c = multiplyStrassen(a, b);
}*/

}  // namespace
