#!/usr/bin/env bash

set -e

cmake -B ./build -DCMAKE_BUILD_TYPE=Release
cmake --build ./build --config Release --parallel 2
./build/test/test_all
./build/benchmark/mat_benchmark --benchmark_format=csv --benchmark_time_unit=ms > build/benchmark.csv

cat ./build/benchmark.csv
python3 ./benchmark/plot.py -f ./build/benchmark.csv --xfieldvar=-1 --ylabel="Run time (ms)" --xlabel="Input size" --output ./benchmark_plot.png
python3 ./benchmark/plot.py -f ./build/benchmark.csv --filter=".*Mult.*" --xfieldvar=-1 --ylabel="Run time (ms)" --xlabel="Input size" --logx --logy --output ./benchmark_multiplication_plot.pdf
python3 ./benchmark/plot.py -f ./build/benchmark.csv --filter=".*Inv.*" --xfieldvar=-1 --ylabel="Run time (ms)" --xlabel="Input size" --logx --logy --output ./benchmark_inversion_plot.pdf
python3 ./benchmark/plot.py -f ./build/benchmark.csv --filter=".*Mult.*" --xfieldvar=-1 --ylabel="GFlops" --xlabel="Input size" --flops="2*input_size**3" --output ./benchmark_multiplication_gflops.pdf