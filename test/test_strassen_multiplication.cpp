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

INSTANTIATE_TEST_SUITE_P(
    TestStrassenMultiplicationAllSizes,
    TestStrassenMultiplication,
    testing::ValuesIn(kTestSizes),
    ProductSizeToString());

}  // namespace
