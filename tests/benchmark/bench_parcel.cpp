
#include <SHARPlib/constants.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/winds.h>
#include <benchmark/benchmark.h>

#include <memory>
#include <iostream>

static float random_lev(const float bot=0, const float top=15000.0) {
	// generate a random level and return it
	float lev = bot + static_cast <float>(std::rand())/(static_cast <float>(RAND_MAX/(top-bot)));
	return lev;
}

auto array_from_range = [](const float bottom, const float top,
                           const std::ptrdiff_t size) {
    auto arr = std::make_unique<float[]>(size);
    float delta = (top - bottom) / static_cast<float>(size);
    for (int k = 0; k < size; ++k) {
        arr[k] = bottom + static_cast<float>(k) * delta;
    }
    return arr;
};

static void lift_pcl_wobf(const float pres[], const std::ptrdiff_t N) {
	const float tmpk = random_lev(15.0f, 30.0f);
    static constexpr sharp::lifter_wobus lifter;

	for (std::ptrdiff_t k = 1; k < N; ++k) {
		float pcl_t = lifter(sharp::THETA_REF_PRESSURE, tmpk, pres[k]);
		benchmark::DoNotOptimize(pcl_t);
	}
}

static void bench_lift_pcl_wobf(benchmark::State& state) {
	// set up and fill our pressure array -- we use make_unique
	// so that memory is managed auromatically instead of manually
	std::ptrdiff_t N = state.range(0);
    auto pres = array_from_range(sharp::THETA_REF_PRESSURE, 500.0f, N);
	// This is the section that will be timed
	for (auto _ : state) {
		lift_pcl_wobf(pres.get(), N);
		benchmark::ClobberMemory();
	}

	state.SetItemsProcessed(state.iterations());
	state.SetComplexityN(state.range(0));
}

BENCHMARK(bench_lift_pcl_wobf)
	->DenseRange(50, 7000, 100)
    ->Complexity();	

// main method macro
BENCHMARK_MAIN();
