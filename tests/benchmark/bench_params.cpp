#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/params/winter.h>
#include <SHARPlib/thermo.h>
#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
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

// Spectral bin classifier (SBC). The sample sounding of its Python reference
// (sample_data.csv; case 1000 of data/sbc_reference): HRRR 2022-02-02 23Z
// f00 at 34.94 N, 97.18 W, 39 levels from the surface (2 m) up, with the
// reference's wet-bulb temperature. The cloud top is at 9496 m, Tw is above
// 0 C only at 1370 and 1617 m, and the core gives ice pellets with each
// drop-size distribution below.
constexpr float SBC_SAMPLE_PRES[] = {
    98210.0f, 97500.0f, 95000.0f, 92500.0f, 90000.0f, 87500.0f, 85000.0f,
    82500.0f, 80000.0f, 77500.0f, 75000.0f, 72500.0f, 70000.0f, 67500.0f,
    65000.0f, 62500.0f, 60000.0f, 57500.0f, 55000.0f, 52500.0f, 50000.0f,
    47500.0f, 45000.0f, 42500.0f, 40000.0f, 37500.0f, 35000.0f, 32500.0f,
    30000.0f, 27500.0f, 25000.0f, 22500.0f, 20000.0f, 17500.0f, 15000.0f,
    12500.0f, 10000.0f, 7500.0f,  5000.0f};
constexpr float SBC_SAMPLE_HGHT[] = {
    0.0f,       56.442383f, 260.70694f, 469.10275f, 681.9126f,  901.26404f,
    1130.5469f, 1369.6487f, 1616.9025f, 1871.5979f, 2133.8792f, 2403.5466f,
    2680.673f,  2966.7834f, 3261.1917f, 3565.882f,  3881.2122f, 4207.1167f,
    4545.4614f, 4895.5923f, 5261.2954f, 5641.041f,  6039.1777f, 6457.2085f,
    6895.511f,  7356.904f,  7843.8374f, 8359.449f,  8908.406f,  9495.977f,
    10127.756f, 10809.0205f, 11555.859f, 12411.309f, 13401.112f, 14553.105f,
    15933.92f,  17685.596f, 20181.354f};
constexpr float SBC_SAMPLE_TMPK[] = {
    269.5766f,  269.07043f, 267.5222f,  265.7832f,  265.22726f, 267.7831f,
    271.326f,   273.35565f, 273.753f,   272.99402f, 271.7636f,  270.32312f,
    268.6287f,  267.1676f,  265.66608f, 264.2079f,  262.53174f, 260.54184f,
    258.4745f,  256.37653f, 254.39236f, 252.43417f, 250.50314f, 248.35762f,
    245.80988f, 242.88681f, 239.50732f, 236.00835f, 232.58592f, 228.92673f,
    224.08499f, 218.03888f, 216.89758f, 220.19148f, 217.91052f, 214.21368f,
    208.7692f,  209.2101f,  211.9147f};
constexpr float SBC_SAMPLE_DWPK[] = {
    267.84027f, 267.2251f,  267.17688f, 265.34628f, 264.62677f, 267.3271f,
    270.72638f, 272.98447f, 273.4998f,  272.6873f,  271.4373f,  270.1248f,
    268.5623f,  267.1248f,  265.3748f,  263.4373f,  261.3748f,  259.1873f,
    256.9998f,  254.56229f, 252.43729f, 250.18729f, 248.12479f, 245.81229f,
    243.06229f, 239.81229f, 236.18729f, 232.49979f, 228.81229f, 224.62479f,
    217.56229f, 203.68729f, 199.31229f, 192.49979f, 192.12479f, 192.12479f,
    192.12479f, 192.12479f, 192.12479f};
constexpr float SBC_SAMPLE_RELH[] = {
    0.874f,      0.8661301f,   0.9670676f,   0.9688999f,   0.9583809f,
    0.9959251f,  0.99636406f,  0.9963686f,   0.9968148f,   0.9961528f,
    0.9861289f,  0.9834164f,   0.9911664f,   0.98461264f,  0.9653304f,
    0.93880683f, 0.9168299f,   0.8952456f,   0.8777837f,   0.8567047f,
    0.8419585f,  0.8228494f,   0.8083167f,   0.7945054f,   0.775455f,
    0.75370675f, 0.7282403f,   0.7048327f,   0.6836666f,   0.6608954f,
    0.6338591f,  0.6003386f,   0.4446618f,   0.036687344f, 0.023816185f,
    0.022249507f, 0.029207146f, 0.021175709f, 0.010973916f};
constexpr float SBC_SAMPLE_WETBULB[] = {
    268.94424f, 268.40955f, 267.4046f,  265.6432f,  265.03622f, 267.6174f,
    271.07037f, 273.18256f, 273.63135f, 272.84793f, 271.612f,   270.23395f,
    268.6001f,  267.1498f,  265.54956f, 263.91086f, 262.10638f, 260.0743f,
    257.99957f, 255.83368f, 253.84691f, 251.84987f, 249.92656f, 247.79037f,
    245.26207f, 242.35379f, 239.02518f, 235.58836f, 232.21382f, 228.5855f,
    223.70995f, 217.5039f,  216.24294f, 218.64377f, 216.56673f, 213.25244f,
    208.24313f, 208.49329f, 210.38112f};
constexpr std::ptrdiff_t SBC_SAMPLE_N = std::size(SBC_SAMPLE_PRES);

// Columns on the levels of the sample. The classifier reads no wet-bulb
// temperature above the cloud top, so where the columns below dry the air
// aloft, they keep the sample's wet-bulb temperature.
enum class SBCColumn {
    // The sample: snow from a frozen cloud top melts in the warm nose and
    // refreezes below it, so the core runs every branch of a frozen top.
    ice_pellets,
    // The sample 1 K colder (T, Td, and Tw), so no level is above 0 C. Snow
    // from the pre-classifier, after one pass up to the cloud top.
    all_subfreezing,
    // The sample with dry air above its warm nose (T - Td of 15 K and
    // relative humidity of 0.3 above 1617 m), so the cloud top has Tw of
    // +0.5 C over a subfreezing surface. Freezing rain from the
    // pre-classifier, before any pass over the column.
    warm_top_cold_surface,
    // The sample 10 K warmer, with dry air above 3261 m, where the cloud top
    // has Tw of +2.4 C. Rain from the pre-classifier, after one pass up to
    // the cloud top.
    all_warm,
};

struct SBCSounding {
    std::array<float, SBC_SAMPLE_N> pres;     // Pa
    std::array<float, SBC_SAMPLE_N> hght;     // m AGL
    std::array<float, SBC_SAMPLE_N> tmpk;     // K
    std::array<float, SBC_SAMPLE_N> dwpk;     // K
    std::array<float, SBC_SAMPLE_N> relh;     // fraction
    std::array<float, SBC_SAMPLE_N> wetbulb;  // K
    float cloud_top;                          // m AGL
    sharp::PrecipType precip_type;
};

static SBCSounding sbc_snd(const SBCColumn column) {
    SBCSounding snd;
    std::copy_n(SBC_SAMPLE_PRES, SBC_SAMPLE_N, snd.pres.begin());
    std::copy_n(SBC_SAMPLE_HGHT, SBC_SAMPLE_N, snd.hght.begin());
    std::copy_n(SBC_SAMPLE_TMPK, SBC_SAMPLE_N, snd.tmpk.begin());
    std::copy_n(SBC_SAMPLE_DWPK, SBC_SAMPLE_N, snd.dwpk.begin());
    std::copy_n(SBC_SAMPLE_RELH, SBC_SAMPLE_N, snd.relh.begin());
    std::copy_n(SBC_SAMPLE_WETBULB, SBC_SAMPLE_N, snd.wetbulb.begin());

    float shift = 0.0f;
    std::ptrdiff_t dry_from = SBC_SAMPLE_N;
    switch (column) {
        case SBCColumn::ice_pellets:
            snd.precip_type = sharp::PrecipType::ice_pellets;
            break;
        case SBCColumn::all_subfreezing:
            shift = -1.0f;
            snd.precip_type = sharp::PrecipType::snow;
            break;
        case SBCColumn::warm_top_cold_surface:
            dry_from = 9;
            snd.precip_type = sharp::PrecipType::freezing_rain;
            break;
        case SBCColumn::all_warm:
            shift = 10.0f;
            dry_from = 15;
            snd.precip_type = sharp::PrecipType::rain;
            break;
    }
    for (std::ptrdiff_t k = 0; k < SBC_SAMPLE_N; ++k) {
        snd.tmpk[k] += shift;
        snd.dwpk[k] += shift;
        snd.wetbulb[k] += shift;
        if (k >= dry_from) {
            snd.dwpk[k] = snd.tmpk[k] - 15.0f;
            snd.relh[k] = 0.3f;
        }
    }
    snd.cloud_top = sharp::spectral_bin_cloud_top(
        snd.pres.data(), snd.hght.data(), snd.tmpk.data(), snd.dwpk.data(),
        snd.relh.data(), SBC_SAMPLE_N);
    return snd;
}

// The default drop-size distribution, and bins of DSD25 (Reeves et al.
// 2016) every 0.1 mm (19 bins, as run_sbc.py builds them for deld = 0.1)
// and every 0.033 mm (64 bins, up to 2.129 mm).
enum class SBCBins { default_4, dsd25_19, dsd25_64 };

struct SBCBinArrays {
    std::array<float, sharp::SBC_MAX_BINS> diameter{};       // mm
    std::array<float, sharp::SBC_MAX_BINS> concentration{};  // per bin
    std::ptrdiff_t nbins = 0;
};

static SBCBinArrays sbc_bins(const SBCBins bins) {
    SBCBinArrays arrays;
    if (bins == SBCBins::default_4) {
        const sharp::SpectralBinDSD dsd = sharp::spectral_bin_dsd_default();
        arrays.diameter = dsd.diameter();
        arrays.concentration = dsd.concentration();
        arrays.nbins = dsd.nbins();
        return arrays;
    }
    // DSD25 as run_sbc.py tables it, 0.05 to 1.65 mm every 0.1 mm. Like
    // run_sbc.py, bins every step mm from 0.05 mm interpolate it linearly,
    // and bins past its end take its last value.
    constexpr float dsd25_diameter[] = {0.05f, 0.15f, 0.25f, 0.35f, 0.45f,
                                        0.55f, 0.65f, 0.75f, 0.85f, 0.95f,
                                        1.05f, 1.15f, 1.25f, 1.35f, 1.45f,
                                        1.55f, 1.65f};
    constexpr float dsd25_concentration[] = {
        55.1843f, 66.0695f, 130.272f, 154.556f, 203.649f, 171.814f,
        206.606f, 146.647f, 94.9404f, 79.4013f, 61.0083f, 35.6567f,
        25.4924f, 16.2522f, 11.6891f, 7.49152f, 3.60886f};
    constexpr std::ptrdiff_t dsd25_N = std::size(dsd25_diameter);
    const float step = (bins == SBCBins::dsd25_19) ? 0.1f : 0.033f;
    arrays.nbins = (bins == SBCBins::dsd25_19) ? 19 : 64;
    for (std::ptrdiff_t j = 0; j < arrays.nbins; ++j) {
        const float D = 0.05f + static_cast<float>(j) * step;
        arrays.diameter[j] = D;
        arrays.concentration[j] = sharp::interp_height(
            std::min(D, dsd25_diameter[dsd25_N - 1]), dsd25_diameter,
            dsd25_concentration, dsd25_N);
    }
    return arrays;
}

static sharp::SpectralBinDSD sbc_dsd(const SBCBins bins) {
    const SBCBinArrays arrays = sbc_bins(bins);
    return sharp::spectral_bin_dsd(arrays.diameter.data(),
                                   arrays.concentration.data(), arrays.nbins);
}

// Skips the benchmark unless both overloads give the column its category.
static bool skip_unless_category(benchmark::State& state,
                                 const SBCSounding& snd,
                                 const sharp::SpectralBinDSD& dsd) {
    const sharp::SpectralBinResult given_top = sharp::spectral_bin_classifier(
        snd.pres.data(), snd.hght.data(), snd.tmpk.data(), snd.dwpk.data(),
        snd.relh.data(), snd.wetbulb.data(), SBC_SAMPLE_N, snd.cloud_top, dsd);
    const sharp::SpectralBinResult full = sharp::spectral_bin_classifier(
        snd.pres.data(), snd.hght.data(), snd.tmpk.data(), snd.dwpk.data(),
        snd.relh.data(), snd.wetbulb.data(), SBC_SAMPLE_N, dsd);
    if ((given_top.precip_type != snd.precip_type) ||
        (full.precip_type != snd.precip_type)) {
        state.SkipWithError("column does not give its precipitation type");
        return true;
    }
    return false;
}

static void bench_spectral_bin_dsd(benchmark::State& state,
                                   const SBCBins bins) {
    const SBCBinArrays arrays = sbc_bins(bins);

    for (auto _ : state) {
        auto dsd = sharp::spectral_bin_dsd(arrays.diameter.data(),
                                           arrays.concentration.data(),
                                           arrays.nbins);
        benchmark::DoNotOptimize(dsd);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
}

static void bench_spectral_bin_cloud_top(benchmark::State& state) {
    const SBCSounding snd = sbc_snd(SBCColumn::ice_pellets);

    for (auto _ : state) {
        auto cloud_top = sharp::spectral_bin_cloud_top(
            snd.pres.data(), snd.hght.data(), snd.tmpk.data(),
            snd.dwpk.data(), snd.relh.data(), SBC_SAMPLE_N);
        benchmark::DoNotOptimize(cloud_top);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
}

// The overload that takes the cloud top.
static void bench_spectral_bin_classifier(benchmark::State& state,
                                          const SBCColumn column,
                                          const SBCBins bins) {
    const SBCSounding snd = sbc_snd(column);
    const sharp::SpectralBinDSD dsd = sbc_dsd(bins);
    if (skip_unless_category(state, snd, dsd)) return;

    for (auto _ : state) {
        auto result = sharp::spectral_bin_classifier(
            snd.pres.data(), snd.hght.data(), snd.tmpk.data(),
            snd.dwpk.data(), snd.relh.data(), snd.wetbulb.data(),
            SBC_SAMPLE_N, snd.cloud_top, dsd);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
}

// The overload that finds the cloud top.
static void bench_spectral_bin_classifier_full(benchmark::State& state,
                                               const SBCColumn column,
                                               const SBCBins bins) {
    const SBCSounding snd = sbc_snd(column);
    const sharp::SpectralBinDSD dsd = sbc_dsd(bins);
    if (skip_unless_category(state, snd, dsd)) return;

    for (auto _ : state) {
        auto result = sharp::spectral_bin_classifier(
            snd.pres.data(), snd.hght.data(), snd.tmpk.data(),
            snd.dwpk.data(), snd.relh.data(), snd.wetbulb.data(),
            SBC_SAMPLE_N, dsd);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
}

static void bench_modified_bourgouin_sbc_column(benchmark::State& state,
                                                const SBCColumn column) {
    const SBCSounding snd = sbc_snd(column);

    for (auto _ : state) {
        auto probs = sharp::modified_bourgouin(
            snd.pres.data(), snd.hght.data(), snd.tmpk.data(),
            snd.dwpk.data(), snd.wetbulb.data(), SBC_SAMPLE_N);
        benchmark::DoNotOptimize(probs);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
}

BENCHMARK_CAPTURE(bench_spectral_bin_dsd, default_4, SBCBins::default_4);
BENCHMARK_CAPTURE(bench_spectral_bin_dsd, dsd25_19, SBCBins::dsd25_19);
BENCHMARK_CAPTURE(bench_spectral_bin_dsd, dsd25_64, SBCBins::dsd25_64);
BENCHMARK(bench_spectral_bin_cloud_top);
BENCHMARK_CAPTURE(bench_spectral_bin_classifier, all_warm, SBCColumn::all_warm,
                  SBCBins::default_4);
BENCHMARK_CAPTURE(bench_spectral_bin_classifier, all_subfreezing,
                  SBCColumn::all_subfreezing, SBCBins::default_4);
BENCHMARK_CAPTURE(bench_spectral_bin_classifier, warm_top_cold_surface,
                  SBCColumn::warm_top_cold_surface, SBCBins::default_4);
BENCHMARK_CAPTURE(bench_spectral_bin_classifier, ice_pellets_default_4,
                  SBCColumn::ice_pellets, SBCBins::default_4);
BENCHMARK_CAPTURE(bench_spectral_bin_classifier, ice_pellets_dsd25_19,
                  SBCColumn::ice_pellets, SBCBins::dsd25_19);
BENCHMARK_CAPTURE(bench_spectral_bin_classifier, ice_pellets_dsd25_64,
                  SBCColumn::ice_pellets, SBCBins::dsd25_64);
BENCHMARK_CAPTURE(bench_spectral_bin_classifier_full, all_warm,
                  SBCColumn::all_warm, SBCBins::default_4);
BENCHMARK_CAPTURE(bench_spectral_bin_classifier_full, all_subfreezing,
                  SBCColumn::all_subfreezing, SBCBins::default_4);
BENCHMARK_CAPTURE(bench_spectral_bin_classifier_full, warm_top_cold_surface,
                  SBCColumn::warm_top_cold_surface, SBCBins::default_4);
BENCHMARK_CAPTURE(bench_spectral_bin_classifier_full, ice_pellets_default_4,
                  SBCColumn::ice_pellets, SBCBins::default_4);
BENCHMARK_CAPTURE(bench_spectral_bin_classifier_full, ice_pellets_dsd25_19,
                  SBCColumn::ice_pellets, SBCBins::dsd25_19);
BENCHMARK_CAPTURE(bench_spectral_bin_classifier_full, ice_pellets_dsd25_64,
                  SBCColumn::ice_pellets, SBCBins::dsd25_64);
BENCHMARK_CAPTURE(bench_modified_bourgouin_sbc_column, all_warm,
                  SBCColumn::all_warm);
BENCHMARK_CAPTURE(bench_modified_bourgouin_sbc_column, all_subfreezing,
                  SBCColumn::all_subfreezing);
BENCHMARK_CAPTURE(bench_modified_bourgouin_sbc_column, warm_top_cold_surface,
                  SBCColumn::warm_top_cold_surface);
BENCHMARK_CAPTURE(bench_modified_bourgouin_sbc_column, ice_pellets,
                  SBCColumn::ice_pellets);

BENCHMARK_MAIN();
