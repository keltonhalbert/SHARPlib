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

/**
 * \brief Fill an array with values between hbot and htop and return a random value in that range
 *
 * Fill an array contained within a unique pointer of size N with values between hbot and htop,
 * and return a random value within that range to test on the array. 
 */
static float fill_array(std::unique_ptr<float[]>& array, const ptrdiff_t N, 
						const float bot=0, const float top=15000.0) {
	const float delta = (top - bot) / N;
	for (ptrdiff_t i = 0; i < N; ++i) {
		array[i] = bot + delta*static_cast<float>(i);
	}

	// generate a random pressure level and return it
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
	auto height = std::make_unique<float[]>(N);
	float lev = fill_array(height, N);

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
	auto height = std::make_unique<float[]>(N);
	float lev = fill_array(height, N);

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
	auto pressure = std::make_unique<float[]>(N);
	float lev = fill_array(pressure, N, 100000.0, 500.0);
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
	auto pressure = std::make_unique<float[]>(N);
	float lev = fill_array(pressure, N, 100000.0, 500.0);
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
	auto height = std::make_unique<float[]>(N);
	float lev = fill_array(height, N);
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
	auto height = std::make_unique<float[]>(N);
	float lev = fill_array(height, N);
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
