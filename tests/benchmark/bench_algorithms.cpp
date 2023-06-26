/**
 *
 *
 *
 */

#include <SHARPlib/algorithms.h>
#include <benchmark/benchmark.h>

#include <memory>
#include <algorithm>
#include <iostream>


auto array_from_range = [](const float bottom, const float top,
				   const std::ptrdiff_t size) {
    auto arr = std::make_unique<float[]>(size);
	float delta = (top - bottom) / static_cast<float>(size);
    for (int k = 0; k < size; ++k) {
        arr[k] = bottom + static_cast<float>(k) * delta;
    }
    return arr;
};

static float random_lev(const float bot=0, const float top=15000.0) {
	// generate a random level and return it
	float lev = bot + static_cast <float>(std::rand())/(static_cast <float>(RAND_MAX/(top-bot)));
	return lev;
}

/**
 * \brief Benchmark the STL lower_bound for comparison against the SHARP implementation
 */
static void bench_std_lower_bound(benchmark::State& state) {
	// set up and fill our pressure array -- we use make_unique
	// so that memory is managed auromatically instead of manually
	ptrdiff_t N = state.range(0);
	auto height = array_from_range(0.0f, 15000.0f, N);
	float lev = random_lev(0.0f, 15000.0f);

	// This is the section that will be timed
	for (auto _ : state) {
		auto lower = std::lower_bound(&height[0], &height[N-1], lev);
		benchmark::DoNotOptimize(lower);
		benchmark::ClobberMemory();
	}

	state.SetItemsProcessed(state.iterations() * state.range(0));
	state.SetComplexityN(state.range(0));
}

/**
 * \brief Benchmark the STL upper_bound for comparison against the SHARP implementation
 */
static void bench_std_upper_bound(benchmark::State& state) {
	// set up and fill our pressure array -- we use make_unique
	// so that memory is managed auromatically instead of manually
	ptrdiff_t N = state.range(0);
	auto height = array_from_range(0.0f, 15000.0f, N);
	float lev = random_lev(0.0f, 15000.0f);

	// This is the section that will be timed
	for (auto _ : state) {
		auto lower = std::upper_bound(&height[0], &height[N-1], lev);
		benchmark::DoNotOptimize(lower);
		benchmark::ClobberMemory();
	}

	state.SetItemsProcessed(state.iterations() * state.range(0));
	state.SetComplexityN(state.range(0));
}

/**
 * \brief Benchmark the SHARP implementation of lower_bound on pressure data 
 */
static void bench_sharp_pressure_lower_bound(benchmark::State& state) {
	// set up and fill our pressure array -- we use make_unique
	// so that memory is managed auromatically instead of manually
	ptrdiff_t N = state.range(0);
	auto pressure = array_from_range(100000.0f, 500.0f, N);
	float lev = random_lev(100000.0f, 500.0f);
	constexpr auto cmp = std::greater<float>();

	// This is the section that will be timed
	for (auto _ : state) {
		// generate a random float between pbot and ptop
		auto lower = sharp::lower_bound(pressure.get(), N, lev, cmp);
		benchmark::DoNotOptimize(lower);
		benchmark::ClobberMemory();
	}

	state.SetItemsProcessed(state.iterations() * state.range(0));
	state.SetComplexityN(state.range(0));
}

/**
 * \brief Benchmark the SHARP implementation of upper_bound on pressure data 
 */
static void bench_sharp_pressure_upper_bound(benchmark::State& state) {
	// set up and fill our pressure array -- we use make_unique
	// so that memory is managed auromatically instead of manually
	ptrdiff_t N = state.range(0);
	auto pressure = array_from_range(100000.0f, 500.0f, N);
	float lev = random_lev(100000.0f, 500.0f);
	constexpr auto cmp = std::greater<float>();

	// This is the section that will be timed
	for (auto _ : state) {
		// generate a random float between pbot and ptop
		auto lower = sharp::upper_bound(pressure.get(), N, lev, cmp);
		benchmark::DoNotOptimize(lower);
		benchmark::ClobberMemory();
	}

	state.SetItemsProcessed(state.iterations() * state.range(0));
	state.SetComplexityN(state.range(0));
}

/**
 * \brief Benchmark the SHARP implementation of lower_bound on height data 
 */
static void bench_sharp_height_lower_bound(benchmark::State& state) {
	// set up and fill our pressure array -- we use make_unique
	// so that memory is managed auromatically instead of manually
	ptrdiff_t N = state.range(0);
	auto height = array_from_range(0.0f, 15000.0f, N);
	float lev = random_lev(0.0f, 15000.0f);
	constexpr auto cmp = std::less<float>();

	// This is the section that will be timed
	for (auto _ : state) {
		// generate a random float between pbot and ptop
		auto lower = sharp::lower_bound(height.get(), N, lev, cmp);
		benchmark::DoNotOptimize(lower);
		benchmark::ClobberMemory();
	}

	state.SetItemsProcessed(state.iterations() * state.range(0));
	state.SetComplexityN(state.range(0));
}

/**
 * \brief Benchmark the SHARP implementation of lower_bound on height data 
 */
static void bench_sharp_height_upper_bound(benchmark::State& state) {
	// set up and fill our pressure array -- we use make_unique
	// so that memory is managed auromatically instead of manually
	ptrdiff_t N = state.range(0);
	auto height = array_from_range(0.0f, 15000.0f, N);
	float lev = random_lev(0.0f, 15000.0f);
	constexpr auto cmp = std::less<float>();

	// This is the section that will be timed
	for (auto _ : state) {
		// generate a random float between pbot and ptop
		auto lower = sharp::upper_bound(height.get(), N, lev, cmp);
		benchmark::DoNotOptimize(lower);
		benchmark::ClobberMemory();
	}

	state.SetItemsProcessed(state.iterations() * state.range(0));
	state.SetComplexityN(state.range(0));
}

// Configure the benchmarks to run over a range of
// array sizes, and attempt to compute the time complexity
BENCHMARK(bench_std_lower_bound)
	->RangeMultiplier(2)
	->Range(4, 2 << 18)
    ->Complexity();	

BENCHMARK(bench_std_upper_bound)
	->RangeMultiplier(2)
	->Range(4, 2 << 18)
    ->Complexity();	

BENCHMARK(bench_sharp_pressure_lower_bound)
	->RangeMultiplier(2)
	->Range(4, 2 << 18)
    ->Complexity();	

BENCHMARK(bench_sharp_pressure_upper_bound)
	->RangeMultiplier(2)
	->Range(4, 2 << 18)
    ->Complexity();	

BENCHMARK(bench_sharp_height_lower_bound)
	->RangeMultiplier(2)
	->Range(4, 2 << 18)
    ->Complexity();	

BENCHMARK(bench_sharp_height_upper_bound)
	->RangeMultiplier(2)
	->Range(4, 2 << 18)
    ->Complexity();	

// main method macro
BENCHMARK_MAIN();
