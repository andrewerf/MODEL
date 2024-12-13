#!/usr/bin/env bash

set -e

cmake -B ./build -DCMAKE_BUILD_TYPE=Release
cmake --build ./build --config Release --parallel 2
#./build/test/test_all
printf "\n\n *** The benchmark has started and may take some time to finish *** \n\n"
#./build/benchmark/mat_benchmark --benchmark_format=csv --benchmark_time_unit=ms | tee build/benchmark.csv

cat ./build/benchmark.csv

# Plot all results on a single graph
python3 ./benchmark/plot.py -f ./build/benchmark.csv --xfieldvar=-1 --ylabel="Run time (ms)" --clean_labels --xlabel="Input size" --output ./benchmark_plot.pdf

# All multiplication algorithms on a log-log plot
python3 ./benchmark/plot.py -f ./build/benchmark.csv --filter=".*Mult.*" --xfieldvar=-1 --ylabel="Run time (ms)" --clean_labels --xlabel="Input size" --logx --logy --output ./benchmark_multiplication_plot.pdf

# All inversion algorithms on a log-log plot
python3 ./benchmark/plot.py -f ./build/benchmark.csv --filter=".*Inv.*" --xfieldvar=-1 --ylabel="Run time (ms)" --clean_labels --xlabel="Input size" --logx --logy --output ./benchmark_inversion_plot.pdf

# GEMM GFlops
python3 ./benchmark/plot.py -f ./build/benchmark.csv --filter=".*Mult.*" --xfieldvar=-1 --ylabel="Effective GFlops" --clean_labels --xlabel="Input size" --flops="2*input_size**3" --output ./benchmark_multiplication_gflops.pdf