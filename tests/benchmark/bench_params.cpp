#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/params/winter.h>
#include <SHARPlib/thermo.h>
#include <benchmark/benchmark.h>

#include <cmath>
#include <iterator>
#include <memory>

// Precipitation-type API (modified Bourgouin method) next to the per-level
// work a caller already pays for: the wet-bulb precompute, and relative
// humidity as the cheap reference.

auto array_from_range = [](const float bottom, const float top,
                           const std::ptrdiff_t size) {
    auto arr = std::make_unique<float[]>(size);
    float delta = (top - bottom) / static_cast<float>(size);
    for (std::ptrdiff_t k = 0; k < size; ++k) {
        arr[k] = bottom + static_cast<float>(k) * delta;
    }
    return arr;
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

// Both soundings share one freezing-rain temperature profile, on N evenly
// spaced levels from the surface to 30 km like a 1 Hz radiosonde: -3 C at
// the surface, a +3 C warm nose at 1.5 km, 6.5 K/km above it, and
// isothermal above an 11 km tropopause. The wet-bulb temperature crosses
// 0 C near 0.8 and 1.9 km, so bourgouin_energy finds a near-surface cold
// layer under a warm layer aloft (RE and ME_aloft about 50 J/kg each, 100 %
// freezing rain) and walks up to the 250 hPa cap near 10.2 km, about a
// third of the levels. The two soundings differ only in moisture, which
// only the relative-humidity walk of precipitation_generation_layer reads.
enum class Moisture {
    // Cloud from the surface to about 4.1 km (the generation layer), a dry
    // layer 4 km deep above it, and a thin moist layer near 9 km. The
    // relative-humidity walk stops where that moist layer starts, near
    // 8.3 km and about 28 % of the levels, because the dry layer below it
    // eliminates everything above: the early exit that dry air over a cloud
    // gives a real sounding.
    cloud_under_dry_layer,
    // A dewpoint depression of 1 K at every level, so no dry layer stops
    // the relative-humidity walk. It reads all N levels, and the generation
    // layer, which layer_min then scans, is the whole column. The worst
    // case for the full-column call.
    moist_column,
};

struct WinterSounding {
    std::unique_ptr<float[]> pres;     // Pa
    std::unique_ptr<float[]> hght;     // m
    std::unique_ptr<float[]> tmpk;     // K
    std::unique_ptr<float[]> dwpk;     // K
    std::unique_ptr<float[]> wetbulb;  // K, Wobus lifter
};

static WinterSounding winter_snd(const std::ptrdiff_t N,
                                 const Moisture moisture) {
    // Piecewise linear in height between these nodes.
    constexpr float tmpc_hght[] = {0.0f, 1500.0f, 11000.0f, 30000.0f};
    constexpr float tmpc_node[] = {-3.0f, 3.0f, -58.75f, -58.75f};
    constexpr float dry_hght[] = {0.0f,    4000.0f, 4500.0f, 8000.0f,
                                  8500.0f, 9500.0f, 10000.0f, 30000.0f};
    constexpr float dry_depression[] = {1.0f, 1.0f, 15.0f, 15.0f,
                                        2.0f, 2.0f, 15.0f, 15.0f};
    constexpr float moist_hght[] = {0.0f, 30000.0f};
    constexpr float moist_depression[] = {1.0f, 1.0f};

    const bool dry_aloft = (moisture == Moisture::cloud_under_dry_layer);
    const float* depr_hght = (dry_aloft) ? dry_hght : moist_hght;
    const float* depr_node = (dry_aloft) ? dry_depression : moist_depression;
    const std::ptrdiff_t depr_N =
        (dry_aloft) ? std::size(dry_hght) : std::size(moist_hght);

    WinterSounding snd;
    snd.hght = array_from_range(0.0f, 30000.0f, N);
    snd.tmpk = std::make_unique<float[]>(N);
    snd.dwpk = std::make_unique<float[]>(N);
    for (std::ptrdiff_t k = 0; k < N; ++k) {
        snd.tmpk[k] =
            sharp::ZEROCNK + sharp::interp_height(snd.hght[k], tmpc_hght,
                                                  tmpc_node,
                                                  std::size(tmpc_hght));
        snd.dwpk[k] = snd.tmpk[k] - sharp::interp_height(snd.hght[k],
                                                         depr_hght, depr_node,
                                                         depr_N);
    }
    snd.pres = pres_dry_snd(100000.0f, snd.hght.get(), snd.tmpk.get(), N);

    snd.wetbulb = std::make_unique<float[]>(N);
    sharp::lifter_wobus lifter;
    for (std::ptrdiff_t k = 0; k < N; ++k) {
        snd.wetbulb[k] =
            sharp::wetbulb(lifter, snd.pres[k], snd.tmpk[k], snd.dwpk[k]);
    }
    return snd;
}

// Skips the benchmark if the sounding lets the API exit early, so every
// timing below includes the crossings, a generation layer, and both passes.
static bool skip_if_early_exit(benchmark::State& state,
                               const WinterSounding& snd,
                               const std::ptrdiff_t N) {
    const sharp::HeightLayer gen_layer = sharp::precipitation_generation_layer(
        snd.pres.get(), snd.hght.get(), snd.tmpk.get(), snd.dwpk.get(), N);
    const sharp::BourgouinEnergy energy = sharp::bourgouin_energy(
        snd.pres.get(), snd.hght.get(), snd.wetbulb.get(), N);
    if ((gen_layer.bottom == sharp::MISSING) ||
        !(energy.melting_energy_aloft > 0.0f) ||
        !(energy.refreezing_energy > 0.0f)) {
        state.SkipWithError("sounding has no generation layer or no warm "
                            "layer over a cold one");
        return true;
    }
    return false;
}

// The per-level precompute a caller does before calling the API, written to
// an output array the way the Python bindings do it.
template <typename Lft>
static void wetbulb_column(Lft& lifter, const float pres[], const float tmpk[],
                           const float dwpk[], float out[],
                           const std::ptrdiff_t N) {
    for (std::ptrdiff_t k = 0; k < N; ++k) {
        out[k] = sharp::wetbulb(lifter, pres[k], tmpk[k], dwpk[k]);
    }
    benchmark::DoNotOptimize(out);
}

static void relh_column(const float pres[], const float tmpk[],
                        const float dwpk[], float out[],
                        const std::ptrdiff_t N) {
    for (std::ptrdiff_t k = 0; k < N; ++k) {
        out[k] = sharp::relative_humidity(pres[k], tmpk[k], dwpk[k]);
    }
    benchmark::DoNotOptimize(out);
}

static void bench_relative_humidity(benchmark::State& state) {
    std::ptrdiff_t N = state.range(0);
    auto snd = winter_snd(N, Moisture::cloud_under_dry_layer);
    auto relh = std::make_unique<float[]>(N);

    for (auto _ : state) {
        relh_column(snd.pres.get(), snd.tmpk.get(), snd.dwpk.get(),
                    relh.get(), N);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
    state.SetComplexityN(state.range(0));
}

static void bench_wetbulb_wobus(benchmark::State& state,
                                const Moisture moisture) {
    std::ptrdiff_t N = state.range(0);
    auto snd = winter_snd(N, moisture);
    auto wetbulb = std::make_unique<float[]>(N);
    sharp::lifter_wobus lifter;

    for (auto _ : state) {
        wetbulb_column(lifter, snd.pres.get(), snd.tmpk.get(), snd.dwpk.get(),
                       wetbulb.get(), N);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
    state.SetComplexityN(state.range(0));
}

static void bench_wetbulb_cm1(benchmark::State& state) {
    std::ptrdiff_t N = state.range(0);
    auto snd = winter_snd(N, Moisture::cloud_under_dry_layer);
    auto wetbulb = std::make_unique<float[]>(N);
    sharp::lifter_cm1 lifter;

    for (auto _ : state) {
        wetbulb_column(lifter, snd.pres.get(), snd.tmpk.get(), snd.dwpk.get(),
                       wetbulb.get(), N);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
    state.SetComplexityN(state.range(0));
}

static void bench_bourgouin_energy(benchmark::State& state) {
    std::ptrdiff_t N = state.range(0);
    auto snd = winter_snd(N, Moisture::cloud_under_dry_layer);
    if (skip_if_early_exit(state, snd, N)) return;

    for (auto _ : state) {
        auto energy = sharp::bourgouin_energy(snd.pres.get(), snd.hght.get(),
                                              snd.wetbulb.get(), N);
        benchmark::DoNotOptimize(energy);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
    state.SetComplexityN(state.range(0));
}

static void bench_precipitation_generation_layer(benchmark::State& state,
                                                 const Moisture moisture) {
    std::ptrdiff_t N = state.range(0);
    auto snd = winter_snd(N, moisture);
    if (skip_if_early_exit(state, snd, N)) return;

    for (auto _ : state) {
        auto gen_layer = sharp::precipitation_generation_layer(
            snd.pres.get(), snd.hght.get(), snd.tmpk.get(), snd.dwpk.get(), N);
        benchmark::DoNotOptimize(gen_layer);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
    state.SetComplexityN(state.range(0));
}

// O(1): N only picks the sounding whose energies, probability of ice, and
// surface wet-bulb temperature feed it.
static void bench_modified_bourgouin_scalar(benchmark::State& state) {
    std::ptrdiff_t N = state.range(0);
    auto snd = winter_snd(N, Moisture::cloud_under_dry_layer);
    if (skip_if_early_exit(state, snd, N)) return;

    const sharp::HeightLayer gen_layer = sharp::precipitation_generation_layer(
        snd.pres.get(), snd.hght.get(), snd.tmpk.get(), snd.dwpk.get(), N);
    const float prob_ice = sharp::probability_of_ice(
        sharp::layer_min(gen_layer, snd.hght.get(), snd.tmpk.get(), N));
    const sharp::BourgouinEnergy energy = sharp::bourgouin_energy(
        snd.pres.get(), snd.hght.get(), snd.wetbulb.get(), N);
    const float surface_wetbulb = snd.wetbulb[0];

    for (auto _ : state) {
        auto probs =
            sharp::modified_bourgouin(energy, prob_ice, surface_wetbulb);
        benchmark::DoNotOptimize(probs);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
    state.SetComplexityN(state.range(0));
}

static void bench_modified_bourgouin_column(benchmark::State& state,
                                            const Moisture moisture) {
    std::ptrdiff_t N = state.range(0);
    auto snd = winter_snd(N, moisture);
    if (skip_if_early_exit(state, snd, N)) return;

    for (auto _ : state) {
        auto probs = sharp::modified_bourgouin(
            snd.pres.get(), snd.hght.get(), snd.tmpk.get(), snd.dwpk.get(),
            snd.wetbulb.get(), N);
        benchmark::DoNotOptimize(probs);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
    state.SetComplexityN(state.range(0));
}

// Gridded columns (50, 100) up to 1 Hz radiosondes (1,000, 6,000), with a
// complexity fit over the four sizes.
static void precip_type_sizes(benchmark::internal::Benchmark* bench) {
    bench->Arg(50)->Arg(100)->Arg(1000)->Arg(6000)->Complexity();
}

BENCHMARK(bench_relative_humidity)->Apply(precip_type_sizes);
BENCHMARK_CAPTURE(bench_wetbulb_wobus, cloud_under_dry_layer,
                  Moisture::cloud_under_dry_layer)
    ->Apply(precip_type_sizes);
BENCHMARK_CAPTURE(bench_wetbulb_wobus, moist_column, Moisture::moist_column)
    ->Apply(precip_type_sizes);
BENCHMARK(bench_wetbulb_cm1)->Apply(precip_type_sizes);
BENCHMARK(bench_bourgouin_energy)->Apply(precip_type_sizes);
BENCHMARK_CAPTURE(bench_precipitation_generation_layer, cloud_under_dry_layer,
                  Moisture::cloud_under_dry_layer)
    ->Apply(precip_type_sizes);
BENCHMARK_CAPTURE(bench_precipitation_generation_layer, moist_column,
                  Moisture::moist_column)
    ->Apply(precip_type_sizes);
BENCHMARK(bench_modified_bourgouin_scalar)->Apply(precip_type_sizes);
BENCHMARK_CAPTURE(bench_modified_bourgouin_column, cloud_under_dry_layer,
                  Moisture::cloud_under_dry_layer)
    ->Apply(precip_type_sizes);
BENCHMARK_CAPTURE(bench_modified_bourgouin_column, moist_column,
                  Moisture::moist_column)
    ->Apply(precip_type_sizes);

BENCHMARK_MAIN();
