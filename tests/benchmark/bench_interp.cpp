
#include <SHARPlib/constants.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <benchmark/benchmark.h>

#include <memory>

static float random_lev(const float bot=0, const float top=15000.0) {
	// generate a random level and return it
	float lev = bot + static_cast <float>(std::rand())/(static_cast <float>(RAND_MAX/(top-bot)));
	return lev;
}

auto array_from_range = [](const float bottom, const float top,
				   const std::ptrdiff_t size) {
    auto arr = std::make_unique<float[]>(size);
	float delta = (top - bottom) / static_cast<float>(size);
    for (std::ptrdiff_t k = 0; k < size; ++k) {
        arr[k] = bottom + static_cast<float>(k) * delta;
    }
    return arr;
};

auto tmpk_snd = [](const float tmpk_sfc, const float delta_tmpk_cap,
                   const float delta_tmpk_troposphere, const float hght_pbl_top,
                   const float hght_trop_top, const float height[],
                   const std::ptrdiff_t size) {
    // get the indices corresponding to the pbl top
    // and the tropopause...
    sharp::HeightLayer free_troposphere =
        sharp::HeightLayer(hght_pbl_top, hght_trop_top);
    sharp::LayerIndex idx =
        sharp::get_layer_index(free_troposphere, height, size);

    const float dse_pbl = sharp::moist_static_energy(0.0f, tmpk_sfc, 0.0f);
    // use the temperature lapse rate to get a
    // dry static energy rate of change
    const float delta_dse_trop =
        sharp::GRAVITY - (sharp::CP_DRYAIR * delta_tmpk_troposphere);

    auto tmpk_arr = std::make_unique<float[]>(size);
    tmpk_arr[0] = tmpk_sfc;

    // calculate the temperature of the boundary layer
    // profile assuming a constant dry static energy
    for (int k = 1; k < idx.kbot; ++k) {
        tmpk_arr[k] = (dse_pbl - height[k] * sharp::GRAVITY) / sharp::CP_DRYAIR;
    }

    // add a capping inversion to the top of the boundary layer
    // and use that as the starting dry static energy for the
    // free troposphere, which will then be modified by the
    // tropospheric lapse rate.
    const float dse_troposphere_bottom = sharp::moist_static_energy(
        height[idx.kbot], tmpk_arr[idx.kbot - 1] + delta_tmpk_cap, 0.0f);

    // calculate the temperature profile of the free troposphere
    for (int k = idx.kbot; k < idx.ktop + 1; ++k) {
        float dse = dse_troposphere_bottom +
                    delta_dse_trop * (height[k] - hght_pbl_top);
        tmpk_arr[k] = (dse - sharp::GRAVITY * height[k]) / sharp::CP_DRYAIR;
    }

    // temperature above the tropopause is isothermal
    for (int k = idx.ktop + 1; k < size; ++k) {
        tmpk_arr[k] = tmpk_arr[idx.ktop];
    }

    return tmpk_arr;
};

auto pres_dry_snd = [](const float pres_sfc, const float height[],
                       const float tmpk[], const std::ptrdiff_t size) {
    auto pres_arr = std::make_unique<float[]>(size);
    pres_arr[0] = pres_sfc;

    for (std::ptrdiff_t k = 1; k < size; ++k) {
        float tmpk_mean = 0.5f * (tmpk[k] + tmpk[k - 1]);
        float delta_z = height[k] - height[k - 1];
        float delta_p =
            -(sharp::GRAVITY * delta_z) / (sharp::RDGAS * tmpk_mean);
        pres_arr[k] = pres_arr[k - 1] * std::exp(delta_p);
    }

    return pres_arr;
};


static void sharp_interp_hght(const float hght[], const float tmpk[], const std::ptrdiff_t N) {
	const float lev = random_lev(0.0f, 15000.0f);
	auto val = sharp::interp_height(lev, hght, tmpk, N);
	benchmark::DoNotOptimize(val);
}

static void sharp_interp_pres(const float pres[], const float tmpk[], const std::ptrdiff_t N) {
	const float lev = random_lev(100000.0f, 500.0f);
	auto val = sharp::interp_pressure(lev, pres, tmpk, N);
	benchmark::DoNotOptimize(val);
}

static void bench_sharp_interp_hght(benchmark::State& state) {
	// set up and fill our pressure array -- we use make_unique
	// so that memory is managed auromatically instead of manually
	ptrdiff_t N = state.range(0);
    //const float pres_sfc = 100000.0f;
    const float tmpk_sfc = 301.5f;
    const float hght_sfc = 0.0f;
    const float hght_top = 15000.0f;
    const float hght_tropopause = 12000.0f;
    const float delta_tmpk_cap = 1.0f;
    const float delta_tmpk_trop = 0.00725f;
    const float hght_pbl_top = 850.0f;

    // initialize the analytical sounding using
    // our input parameters
    auto hght = array_from_range(hght_sfc, hght_top, N);
    auto tmpk = tmpk_snd(tmpk_sfc, delta_tmpk_cap, delta_tmpk_trop,
                         hght_pbl_top, hght_tropopause, hght.get(), N);

	// This is the section that will be timed
	for (auto _ : state) {
		sharp_interp_hght(hght.get(), tmpk.get(), N);
		benchmark::ClobberMemory();
	}

	state.SetItemsProcessed(state.iterations());
	state.SetComplexityN(state.range(0));
}

static void bench_sharp_interp_pres(benchmark::State& state) {
	// set up and fill our pressure array -- we use make_unique
	// so that memory is managed auromatically instead of manually
	ptrdiff_t N = state.range(0);
    const float pres_sfc = 100000.0f;
    const float tmpk_sfc = 301.5f;
    const float hght_sfc = 0.0f;
    const float hght_top = 15000.0f;
    const float hght_tropopause = 12000.0f;
    const float delta_tmpk_cap = 1.0f;
    const float delta_tmpk_trop = 0.00725f;
    const float hght_pbl_top = 850.0f;

    // initialize the analytical sounding using
    // our input parameters
    auto hght = array_from_range(hght_sfc, hght_top, N);
    auto tmpk = tmpk_snd(tmpk_sfc, delta_tmpk_cap, delta_tmpk_trop,
                         hght_pbl_top, hght_tropopause, hght.get(), N);
    auto pres = pres_dry_snd(pres_sfc, hght.get(), tmpk.get(), N);

	// This is the section that will be timed
	for (auto _ : state) {
		sharp_interp_pres(pres.get(), tmpk.get(), N);
		benchmark::ClobberMemory();
	}

	state.SetItemsProcessed(state.iterations());
	state.SetComplexityN(state.range(0));
}

// Configure the benchmarks to run over a range of
// array sizes, and attempt to compute the time complexity
BENCHMARK(bench_sharp_interp_hght)
	->DenseRange(50, 7000, 100)
    ->Complexity();	

BENCHMARK(bench_sharp_interp_pres)
	->DenseRange(50, 7000, 100)
    ->Complexity();	


BENCHMARK_MAIN();
