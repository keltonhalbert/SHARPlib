
#include <SHARPlib/constants.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/winds.h>
#include <benchmark/benchmark.h>

#include <memory>

auto array_from_range = [](const float bottom, const float top,
                           const std::ptrdiff_t size) {
    auto arr = std::make_unique<float[]>(size);
    float delta = (top - bottom) / static_cast<float>(size);
    for (int k = 0; k < size; ++k) {
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

auto pres_moist_snd = [](const float pres_sfc, const float height[],
                         const float tmpk[], const float mixr[],
                         const std::ptrdiff_t size) {
    auto pres_arr = std::make_unique<float[]>(size);
    pres_arr[0] = pres_sfc;

    for (std::ptrdiff_t k = 1; k < size; ++k) {
        const float tmpk_mean = 0.5f * (tmpk[k] + tmpk[k - 1]);
        const float mixr_mean = 0.5f * (mixr[k] + mixr[k - 1]);
        const float tv_mean = sharp::virtual_temperature(tmpk_mean, mixr_mean);
        const float delta_z = height[k] - height[k - 1];
        const float delta_p =
            -(sharp::GRAVITY * delta_z) / (sharp::RDGAS * tv_mean);
        pres_arr[k] = pres_arr[k - 1] * std::exp(delta_p);
    }

    return pres_arr;
};

auto mixr_snd = [](const float pres_sfc, const float mixr_sfc,
                   const float relh_troposphere, const float height[],
                   const float tmpk[], const std::ptrdiff_t size) {
    auto pres_arr = pres_dry_snd(pres_sfc, height, tmpk, size);
    auto mixr_arr = std::make_unique<float[]>(size);
    mixr_arr[0] = mixr_sfc;

    for (std::ptrdiff_t k = 1; k < size; ++k) {
        const float sat_mixr = sharp::mixratio(pres_arr[k], tmpk[k]);
        mixr_arr[k] = relh_troposphere * sat_mixr;
    }

    return mixr_arr;
};

auto dwpk_snd = [](const float pres[], const float mixr[],
                   const std::ptrdiff_t size) {
    auto dwpk_arr = std::make_unique<float[]>(size);
    for (std::ptrdiff_t k = 0; k < size; ++k) {
        dwpk_arr[k] = sharp::temperature_at_mixratio(mixr[k], pres[k]);
    }
    return dwpk_arr;
};

auto vtmpk_snd = [](const float tmpk[], const float mixr[],
                    const std::ptrdiff_t size) {
    auto vtmpk_arr = std::make_unique<float[]>(size);
    for (std::ptrdiff_t k = 0; k < size; ++k) {
        vtmpk_arr[k] = sharp::virtual_temperature(tmpk[k], mixr[k]);
    }
    return vtmpk_arr;
};

static void compute_sbcape_wobf(const float pres[], const float hght[], 
							    const float tmpk[], const float vtmpk[],
								const float dwpk[], float buoy[], 
								const std::ptrdiff_t N) {
    sharp::Parcel sfc_pcl;
    sfc_pcl.pres = pres[0];
    sfc_pcl.tmpk = tmpk[0];
    sfc_pcl.dwpk = dwpk[0];

    static constexpr sharp::lifter_wobus lifter;
    // lift and integrate the surface parcel
    sharp::lift_parcel(lifter, pres, vtmpk, buoy, N,
                       &sfc_pcl);
	benchmark::DoNotOptimize(sfc_pcl);
	benchmark::DoNotOptimize(hght);
	benchmark::DoNotOptimize(buoy);
    //sharp::cape_cinh(pres, hght, buoy, N, &sfc_pcl);	
	//float cape = sfc_pcl.cape;
	//benchmark::DoNotOptimize(cape);
}

static void compute_sbcape_cm1_pseudo_ice(const float pres[], const float hght[], 
							    		  const float tmpk[], const float vtmpk[],
										  const float dwpk[], float buoy[], 
								          const std::ptrdiff_t N) {
    sharp::Parcel sfc_pcl;
    sfc_pcl.pres = pres[0];
    sfc_pcl.tmpk = tmpk[0];
    sfc_pcl.dwpk = dwpk[0];

    sharp::lifter_cm1 lifter;
	lifter.ma_type = sharp::adiabat::pseudo_ice;

    // lift and integrate the surface parcel
    sharp::lift_parcel(lifter, pres, vtmpk, buoy, N,
                       &sfc_pcl);
	benchmark::DoNotOptimize(sfc_pcl);
	benchmark::DoNotOptimize(hght);
	benchmark::DoNotOptimize(buoy);
    //sharp::cape_cinh(pres, hght, buoy, N, &sfc_pcl);	
	//float cape = sfc_pcl.cape;
	//benchmark::DoNotOptimize(cape);
}

static void bench_sbcape_wobf(benchmark::State& state) {
	// set up and fill our pressure array -- we use make_unique
	// so that memory is managed auromatically instead of manually
	std::ptrdiff_t N = state.range(0);
    const float pres_sfc = 100000.0f;
    const float tmpk_sfc = 301.5f;
    const float mixr_sfc = 0.0157f;
    const float relh_troposphere = 0.85f;
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
    auto mixr = mixr_snd(pres_sfc, mixr_sfc, relh_troposphere, hght.get(),
                         tmpk.get(), N);
    auto pres = pres_moist_snd(pres_sfc, hght.get(), tmpk.get(), mixr.get(), N);
    auto dwpk = dwpk_snd(pres.get(), mixr.get(), N);
    auto vtmpk = vtmpk_snd(tmpk.get(), mixr.get(), N);

    // allocate an array for our buoyancy data
    auto buoy = std::make_unique<float[]>(N);

	// This is the section that will be timed
	for (auto _ : state) {
		compute_sbcape_wobf(pres.get(), hght.get(), tmpk.get(), vtmpk.get(), dwpk.get(), buoy.get(), N);
		benchmark::ClobberMemory();
	}

	state.SetItemsProcessed(state.iterations());
	state.SetComplexityN(state.range(0));
}

static void bench_sbcape_cm1_pseudo_ice(benchmark::State& state) {
	// set up and fill our pressure array -- we use make_unique
	// so that memory is managed auromatically instead of manually
	std::ptrdiff_t N = state.range(0);
    const float pres_sfc = 100000.0f;
    const float tmpk_sfc = 301.5f;
    const float mixr_sfc = 0.0157f;
    const float relh_troposphere = 0.85f;
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
    auto mixr = mixr_snd(pres_sfc, mixr_sfc, relh_troposphere, hght.get(),
                         tmpk.get(), N);
    auto pres = pres_moist_snd(pres_sfc, hght.get(), tmpk.get(), mixr.get(), N);
    auto dwpk = dwpk_snd(pres.get(), mixr.get(), N);
    auto vtmpk = vtmpk_snd(tmpk.get(), mixr.get(), N);

    // allocate an array for our buoyancy data
    auto buoy = std::make_unique<float[]>(N);

	// This is the section that will be timed
	for (auto _ : state) {
		compute_sbcape_cm1_pseudo_ice(pres.get(), hght.get(), tmpk.get(), vtmpk.get(), dwpk.get(), buoy.get(), N);
		benchmark::ClobberMemory();
	}

	state.SetItemsProcessed(state.iterations());
	state.SetComplexityN(state.range(0));
}


BENCHMARK(bench_sbcape_wobf)
//	->RangeMultiplier(2)
	->DenseRange(50, 7000, 100)
    ->Complexity();	

BENCHMARK(bench_sbcape_cm1_pseudo_ice)
//	->RangeMultiplier(2)
	->DenseRange(50, 7000, 100)
    ->Complexity();	

// main method macro
BENCHMARK_MAIN();
