#!/usr/bin/env bash

set -e

cmake -B ./build -DCMAKE_BUILD_TYPE=Release
cmake --build ./build --config Release --parallel 2
./build/test/test_all
./build/benchmark/mat_benchmark --benchmark_format=csv --benchmark_time_unit=ms > build/benchmark.csv

cat ./build/benchmark.csv
python3 ./benchmark/plot.py -f ./build/benchmark.csv --xfieldvar=-1 --output ./benchmark_plot.png
