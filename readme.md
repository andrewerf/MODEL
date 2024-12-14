# MU4IN901 - Modèles de calcul - Implementation Project

**Group #12**: Emilie Gillet, Ruchi Thareja, Andrei Aralov

![benchmark_plot.png](resource%2Fbenchmark_plot.png)

## Build instructions

The library is header-only. Tests and benchmarks are compiled to separate executables.

### Requirements

* CMake
* GoogleTest and Google Benchmark (these libraries are automatically installed by CMake)
* Optionally, python, matplotlib and pandas for plotting the graphs

### End to end script

The most convenient way to get the tests and benchmarks results is to launch all-in-one script `./scripts/build_test_benchmark.sh` from the root of the repo.
The script outputs results of testing and benchmarks (see picture). Additionally, if python is installed, it plots the results of the benchmarks.
The list of generated files:
- `build/*` - build files of the project
- `build/benchmark.csv` - benchmark results in CSV format
- `benchmark_plot.pdf`, `benchmark_multiplication_plot.pdf`, `benchmark_inversion_plot.pdf`, `benchmark_multiplication_gflops.pdf` -
    different plots of the benchmarks results (generated automatically from the CSV file)

### Manual build
You can as well build the project manually using CMake. For reference see the script mentioned above, but in a nutshell, just launch from the project-root:
```shell
cmake -B ./build -DCMAKE_BUILD_TYPE=Release
cmake --build ./build --config Release --parallel 2
```
This will build the project into directory `./build`. Test and benchmark executables can be found in `./build/test/test_all` and `./build/benchmark/mat_benchmark` correspondingly. 
