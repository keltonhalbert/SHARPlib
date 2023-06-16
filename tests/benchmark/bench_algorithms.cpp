/**
 *
 *
 *
 */


#include <SHARPlib/algorithms.h>
#include <benchmark/benchmark.h>

#include <memory>

static void i32_addition(benchmark::State &state) {
    int32_t a, b, c;
	a = std::rand();
	b = std::rand();
    for (auto _ : state)
		benchmark::DoNotOptimize(c = (++a) + (++b));
}

BENCHMARK(i32_addition);
BENCHMARK_MAIN();
