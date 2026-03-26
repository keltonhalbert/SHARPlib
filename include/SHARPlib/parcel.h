/**
 * \file
 * \brief Routines used for parcel lifting and integration
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2022-11-09
 *
 * Written for the NWS Storm Prediction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */
#ifndef SHARP_PARCEL_H
#define SHARP_PARCEL_H

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/thermo.h>
#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace sharp {

////////////    FUNCTORS    ///////////
//

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief A functor that calls the Wobus Wetlift funtion
 *
 * This functor is used to wrap the Wobus Wetlift function for parcel
 * lifting routines. Functors - classes with their operator()
 * overloaded - are used so that functions can be
 * passed to templates in a way that the compiler can still
 * optimize, rather than using function pointers or lambdas.
 *
 * Specifically, this functor is designed to be passed as a template
 * argument to sharp::Parcel::lift_parcel, so that the method of computing
 * moist adiabats can be changed without changing the overall parcel
 * lifting code. The reason this is awesome is that the compiler
 * can still optimize and inline this code, while the user can configure
 * the parcel lifting algorithm to their specifications.
 */
struct lifter_wobus {
    static constexpr bool lift_from_lcl = true;
    static constexpr adiabat ma_type = adiabat::pseudo_liq;

    /**
     * \brief The iterative convergence criteria (K)
     */
    float converge = 0.001f;

    /**
     * \brief Perform the setup step for the parcel lifter.
     *
     * Some parcel lifters require setup in order to handle
     * adiabatic ascent and tracking of vapor, liquid, and
     * ice mixing ratios. The Wobus lifter does not require
     * this, however, so this function does nothing.
     */
    inline void setup([[maybe_unused]] const float lcl_pres,
                      [[maybe_unused]] const float lcl_tmpk) const {}

    /**
     * \brief Overloads operator() to call sharp::wetlift.
     * \param   pres        Parcel pressure (Pa)
     * \param   tmpk        Parcel temperature (K)
     * \param   new_pres    Final level of parcel after lift (Pa)
     *
     * \return  The temperature of the lifted parcel
     */
    [[nodiscard]] inline float operator()(const float pres, const float tmpk,
                                          const float new_pres) const {
        return wetlift(pres, tmpk, new_pres, converge);
    }

    /**
     * \brief Computes the virtual temperature of the parcel
     * \param   pres        Parcel pressure (Pa)
     * \param   tmpk        Parcel temperature (K)
     *
     * \return  The virtual temperature of the parcel
     */
    [[nodiscard]] inline float parcel_virtual_temperature(
        const float pres, const float tmpk) const {
        return virtual_temperature(tmpk, mixratio(pres, tmpk));
    }
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief A functor that calls the CM1 moist adiabat solver

 * This functor is used to wrap the CM1 moist adiabat function
 * for parcel lifting routines. Functors - classes with their
 * operator() overloaded - are used so that functions can be
 * passed to templates in a way that the compiler can still
 * optimize, rather than using function pointers or lambdas.
 *
 * Specifically, this functor is designed to be passed as a template
 * argument to sharp::Parcel::lift_parcel, so that the method of computing
 * moist adiabats can be changed without changing the overall parcel
 * lifting code. The reason this is awesome is that the compiler
 * can still optimize and inline this code, while the user can configure
 * the parcel lifting algorithm to their specifications.
 */
struct lifter_cm1 {
    static constexpr bool lift_from_lcl = false;

    /**
     * \brief The type of moist adiabat to use, as defined by sharp::adiabat
     */
    adiabat ma_type = adiabat::pseudo_liq;

    /**
     * \brief The pressure increment (Pa) to use for the iterative solver
     */
    float pressure_incr = 500.0f;

    /**
     * \brief The iterative convergence criteria (K)
     */
    float converge = 0.001f;

    /**
     * \brief Used to keep track of mixing ratio for conserved/adiabatic lifting
     */
    float rv_total = MISSING;

    /**
     * \brief Water vapor mixing ratio variable updated during parcel lifts
     */
    float rv = MISSING;

    /**
     * \brief Liquid water mixing ratio variable updated during parcel lifts
     */
    float rl = MISSING;

    /**
     * \brief Ice water mixing ratio variable updated during parcel lifts
     */
    float ri = MISSING;

    /**
     * \brief perform the necessary setup for parcel ascent.
     *
     * This function sets the total water mixing ratio for
     * adiabatic parcel ascent, and zeroes out the vapor,
     * liquid, and ice mixing ratios from previous parcel
     * ascents.
     */
    inline void setup(const float lcl_pres, const float lcl_tmpk) {
        this->rv_total = mixratio(lcl_pres, lcl_tmpk);
        this->rv = 0.0f;
        this->rl = 0.0f;
        this->ri = 0.0f;
    }

    /**
     * \brief Overloads operator() to call sharp::moist_adiabat_cm1
     *
     * \param   pres        Parcel pressure (Pa)
     * \param   tmpk        Parcel temperature (K)
     * \param   new_pres    Final level of parcel after lift (Pa)
     *
     * \return  The temperature of the lifted parcel
     */
    [[nodiscard]] inline float operator()(const float pres, const float tmpk,
                                          const float new_pres) {
        return moist_adiabat_cm1(pres, tmpk, new_pres, this->rv_total, this->rv,
                                 this->rl, this->ri, this->pressure_incr,
                                 this->converge, this->ma_type);
    }

    /**
     * \brief Computes the virtual temperature of the parcel
     * \param   pres        Parcel pressure (Pa)
     * \param   tmpk        Parcel temperature (K)
     *
     * \return  The virtual temperature of the parcel
     */
    [[nodiscard]] inline float parcel_virtual_temperature(
        [[maybe_unused]] const float pres, const float tmpk) const {
        return virtual_temperature(tmpk, this->rv, this->rl, this->ri);
    }
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Pseudoadiabatic ascent lookup table (LUT)
 *
 * This data structure contains everything necessary to construct,
 * store, and interrogate a pseudoadiabatic ascent lookup table (LUT).
 * The LUT is generated using the passed lifter (e.g. sharp::lifter_cm1).
 * Due to the complex nature of reversible adiabatic ascent, it is
 * not supported in this lookup table. Any sharp::adiabat::adiab_liq
 * or sharp::adiabat::adiab_ice configured lifters will throw an error.
 *
 * This is intended to be used by sharp::lifter_lut for actual parcel
 * lifting operations. While this data structure can be accessed by
 * multiple threads, each thread must have its own sharp::lifter_lut.
 */
template <typename Lft>
struct lut_data {
    /**
     * \brief The moist adiabatic ascent solver to use
     */
    Lft m_lifter;

    /**
     * \brief The minimum pressure of the lookup table
     */
    const float pres_min;

    /**
     * \brief The maximum pressure of the lookup table
     */
    const float pres_max;

    /**
     * \brief The minimum thetae of the lookup table
     */
    const float thetae_min;

    /**
     * \brief The maximum thetae of the lookup table
     */
    const float thetae_max;

    /**
     * \brief The number of pressure levels (in log space) for the LUT
     */
    const std::size_t num_logp;

    /**
     * \brief The number of thetae levels for the LUT
     */
    const std::size_t num_thetae;

    /**
     * \brief Storage for the maximum logp
     */
    const float m_logp_max;

    /**
     * \brief Storage for inverse delta theta
     */
    const float m_delta_thetae_inv;

    /**
     * \brief Storage for inverse delta logarithmic pressure
     */
    const float m_delta_logp_inv;

    /**
     * \brief The parcel temperature LUT
     */
    std::vector<float> m_LUT_pcl_tmpk;

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Create a pseudoadiabatic lookup table (LUT)
     *
     * Construct the pseudoadiabatic lookup table (LUT) using
     * the provided lifter and LUT configuration options. This
     * lookup data is safe to share with many threads, but each
     * thread must use its own sharp::lifter_lut.
     *
     * \param   lifter      Parcel lifting function/functor
     * \param   pmin        Minimum Pressure (Pa)
     * \param   pmax        Maximum pressure (Pa)
     * \param   thte_min    Minimum theta-e (K)
     * \param   thte_max    Maximum theta-e (K)
     * \param   n_logp      Number of logarithmic pressure levels
     * \param   n_thetae    Number of theta-e levels
     */
    explicit lut_data(Lft lifter, float pmin = 5000.0f, float pmax = 110000.0f,
                      float thte_min = 210.0f, float thte_max = 430.0f,
                      std::size_t n_logp = 201, std::size_t n_thetae = 221)
        : m_lifter(std::move(lifter)),
          pres_min(pmin),
          pres_max(pmax),
          thetae_min(thte_min),
          thetae_max(thte_max),
          num_logp(n_logp),
          num_thetae(n_thetae),
          m_logp_max(std::log(pmax)),
          m_delta_thetae_inv(static_cast<float>(n_thetae - 1) /
                             (thte_max - thte_min)),
          m_delta_logp_inv(static_cast<float>(n_logp - 1) /
                           (std::log(pmin) - std::log(pmax))) {
        if (m_lifter.ma_type == sharp::adiabat::adiab_liq ||
            m_lifter.ma_type == sharp::adiabat::adiab_ice) {
            throw std::invalid_argument(
                "CRITICAL ERROR: sharp::lut_data does not support reversible "
                "adiabatic ascent (adiab_liq or adiab_ice). The 2D lookup "
                "table cannot track conserved total water mass. Please "
                "configure the lifter to use pseudo_liq or pseudo_ice "
                "instead.");
        }
        generate_table();
    }

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Get the bounding vertical indices for a given pressure
     *
     * Returns the bounding vertical indices for the LUT given a
     * pressure value, along with the interpolation weight.
     *
     * \param   pres    Pa
     *
     * \return  tuple of (k0, k1, weight_k): bottom index, top index,
     *          and interpolation weight
     */
    [[nodiscard]] inline std::tuple<std::size_t, std::size_t, float>
    get_logp_indices(const float pres) const {
        const float target_logp = std::log(pres);

        const float frac_k =
            std::max((target_logp - m_logp_max) * m_delta_logp_inv, 0.0f);
        const std::size_t k0 =
            std::min(static_cast<std::size_t>(frac_k), num_logp - 2);

        const std::size_t k1 = k0 + 1;
        const float weight_k =
            std::clamp(frac_k - static_cast<float>(k0), 0.0f, 1.0f);

        return {k0, k1, weight_k};
    }

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Get the fractional thetae coordinate index
     *
     * Returns the fractional index for the thetae coordinate
     * of the LUT given the lcl temperature and the bounding
     * vertical coordinate indices, along with the vertical
     * interpolation weight.
     *
     * \param   lcl_tmpk    lifted condensation level temperature (K)
     * \param   k0          The bottom bounding index of pressure
     * \param   k1          The top bounding index of pressure
     * \param   weight_k    The vertical interpolation weight
     *
     * \return  fractional index of thetae axis for LUT
     */
    [[nodiscard]] inline float find_thetae_index(const float lcl_tmpk,
                                                 const std::size_t k0,
                                                 const std::size_t k1,
                                                 const float weight_k) const {
        float frac_i_lo = 0.0f;
        float frac_i_hi = static_cast<float>(num_thetae - 1);
        float frac_i_mid = 0.0f;

        constexpr int max_iter = 20;
        constexpr float tol = 0.005f;  // Kelvin
        const float* lut = m_LUT_pcl_tmpk.data();
        for (int step = 0; step < max_iter; ++step) {
            frac_i_mid = 0.5f * (frac_i_lo + frac_i_hi);

            const std::size_t i0 =
                std::min(static_cast<std::size_t>(frac_i_mid), num_thetae - 2);
            const std::size_t i1 = i0 + 1;
            const float weight_i =
                std::clamp(frac_i_mid - static_cast<float>(i0), 0.0f, 1.0f);

            const std::size_t idx1 = i0 * num_logp + k0;
            const std::size_t idx2 = i0 * num_logp + k1;
            const std::size_t idx3 = i1 * num_logp + k0;
            const std::size_t idx4 = i1 * num_logp + k1;

            const float t0 = sharp::lerp(lut[idx1], lut[idx2], weight_k);
            const float t1 = sharp::lerp(lut[idx3], lut[idx4], weight_k);
            const float t_interp = sharp::lerp(t0, t1, weight_i);

            const float residual = t_interp - lcl_tmpk;
            if (std::abs(residual) < tol) break;

            if (residual > 0.0f) {
                frac_i_hi = frac_i_mid;
            } else {
                frac_i_lo = frac_i_mid;
            }
        }
        return frac_i_mid;
    }

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Perform bilinear interpolation of LUT
     *
     * Performs the bilinear interpolation of the LUT to
     * return a parcel temperature along a pseudoadiabat.
     *
     * \param   i0  thetae left bounding index
     * \param   i1  thetae right bounding index
     * \param   wi  thetae interpolation weight
     * \param   k0  pressure bottom bounding index
     * \param   k1  pressure top bounding index
     * \param   wk  pressure interpolation weight
     *
     * \return  parcel temperature (K)
     */
    [[nodiscard]] inline float bilinear_interp(
        const std::size_t i0, const std::size_t i1, const float wi,
        const std::size_t k0, const std::size_t k1, const float wk) const {
        const float* lut = m_LUT_pcl_tmpk.data();
        const float val00 = lut[i0 * num_logp + k0];
        const float val01 = lut[i0 * num_logp + k1];
        const float val10 = lut[i1 * num_logp + k0];
        const float val11 = lut[i1 * num_logp + k1];

        const float val0 = sharp::lerp(val00, val01, wk);
        const float val1 = sharp::lerp(val10, val11, wk);
        return sharp::lerp(val0, val1, wi);
    }

   private:
    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Constructs pseudoadiabatic LUT
     *
     * Called by the constructor, this creates and fills
     * the pseudoadiabatic parcel ascent LUT.
     */
    inline void generate_table() {
        std::vector<float> logp_coord(num_logp);
        std::vector<float> thetae_coord(num_thetae);

        const float delta_thetae = (thetae_max - thetae_min) / (num_thetae - 1);
        for (std::size_t i = 0; i < num_thetae; ++i) {
            thetae_coord[i] = thetae_min + i * delta_thetae;
        }

        const float logp_min = std::log(pres_min);
        const float delta_logp = (logp_min - m_logp_max) / (num_logp - 1);
        for (std::size_t i = 0; i < num_logp; ++i) {
            logp_coord[i] = m_logp_max + i * delta_logp;
        }

        const std::size_t size_2D = num_logp * num_thetae;
        m_LUT_pcl_tmpk.resize(size_2D);

        for (std::size_t i = 0; i < num_thetae; ++i) {
            const float thetae_target = thetae_coord[i];
            const float tmpk = solve_tmpk_for_thetae(pres_max, thetae_target);

            m_lifter.setup(pres_max, tmpk);
            float pres_curr = pres_max;
            float tmpk_curr = tmpk;

            for (std::size_t k = 0; k < num_logp; ++k) {
                const float pres = std::exp(logp_coord[k]);
                if (k > 0) {
                    tmpk_curr = m_lifter(pres_curr, tmpk_curr, pres);
                }

                const std::size_t idx_2D = i * num_logp + k;
                m_LUT_pcl_tmpk[idx_2D] = tmpk_curr;

                pres_curr = pres;
            }
        }
    }

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Finds and returns a temperature for a given pressure and thetae.
     *
     * Iteratively inverts theta-e to find the saturated temperature at a
     * single pressure level. Used for finding the starting parcel
     * temperature during LUT construction.
     *
     * \param   pressure            (Pa)
     * \param   thetae_target       (K)
     *
     * \return  lcl_temperature     (K)
     */
    [[nodiscard]] float solve_tmpk_for_thetae(float pressure,
                                              float thetae_target) const {
        float tmpk_lo = 200.0f;
        float tmpk_hi = 330.0f;
        float tmpk_mid = ZEROCNK;

        const float tol = 0.001f;
        const std::size_t max_iters = 50;

        for (std::size_t i = 0; i < max_iters; ++i) {
            tmpk_mid = 0.5f * (tmpk_lo + tmpk_hi);
            const float curr_thetae =
                sharp::thetae(pressure, tmpk_mid, tmpk_mid);

            if (std::abs(curr_thetae - thetae_target) < tol) {
                return tmpk_mid;
            }

            if (curr_thetae < thetae_target) {
                tmpk_lo = tmpk_mid;
            } else {
                tmpk_hi = tmpk_mid;
            }
        }
        fmt::println("Failed to converge on a LCL temperature for thetae.");
        return tmpk_mid;
    }
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Parcel lifter functor for LUT based calculations.
 *
 * This functor is used to wrap sharp::lut_data, which creates
 * a pseudoadiabatic ascent lookup table based on the type of
 * parcel lifter used to construct the table. Instead of
 * directly solving the moist ascent ODEs, this lifter uses
 * the LUT to perform bilinear interpolation to get the
 * parcel temperature.
 *
 * LUT based parcel ascent only works for pseudoadiabats.
 * Constructing the LUT with a reversible adiabat type will
 * result in an error being thrown.
 */
template <typename Lft>
struct lifter_lut {
    static constexpr bool lift_from_lcl = Lft::lift_from_lcl;

    /**
     * \brief Shared pointer to a previously created LUT
     */
    std::shared_ptr<const lut_data<Lft>> m_data;

    /**
     * Used to signal a fallback to the direct solver
     */
    bool m_use_lifter = false;

    /**
     * Fractional index used to determine the pseudoadiabat
     */
    float m_fi = 0.0f;

    /**
     * An instantiation of a parcel lifter for fallback computation
     */
    Lft m_lifter;

    explicit lifter_lut(std::shared_ptr<const lut_data<Lft>> data)
        : m_data(std::move(data)), m_lifter(m_data->m_lifter) {}

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Performs a setup step based on the LCL attributes
     *
     * Computes the fractional index needed to select the
     * correct pseudoadiabat for lookup. If the LCL is outside
     * the table bounds, it falls back to the direct solver.
     */
    inline void setup(const float lcl_pres, const float lcl_tmpk) {
        const auto [k0, k1, wk] = m_data->get_logp_indices(lcl_pres);
        const std::size_t i_max = m_data->num_thetae - 1;

        const float t_min_bound =
            sharp::lerp(m_data->m_LUT_pcl_tmpk[0 * m_data->num_logp + k0],
                        m_data->m_LUT_pcl_tmpk[0 * m_data->num_logp + k1], wk);

        const float t_max_bound = sharp::lerp(
            m_data->m_LUT_pcl_tmpk[i_max * m_data->num_logp + k0],
            m_data->m_LUT_pcl_tmpk[i_max * m_data->num_logp + k1], wk);

        if (lcl_tmpk < t_min_bound || lcl_tmpk > t_max_bound ||
            lcl_pres < m_data->pres_min || lcl_pres > m_data->pres_max) {
            m_use_lifter = true;
            m_lifter.setup(lcl_pres, lcl_tmpk);
            return;
        }

        m_use_lifter = false;
        m_fi = m_data->find_thetae_index(lcl_tmpk, k0, k1, wk);
    }

    /**
     * \brief Overloads operator() to perform LUT interpolation
     *
     * \param   pres        Parcel pressure (Pa)
     * \param   tmpk        Parcel temperature (K)
     * \param   new_pres    Final level of parcel after lift (Pa)
     *
     * \return  The temperature of the lifted parcel
     */
    [[nodiscard]] inline float operator()(const float pres, const float tmpk,
                                          const float new_pres) {
        if (m_use_lifter) {
            return m_lifter(pres, tmpk, new_pres);
        }

        if ((new_pres < m_data->pres_min) || (new_pres > m_data->pres_max))
            return MISSING;

        const std::size_t i0 =
            std::min(static_cast<std::size_t>(m_fi), m_data->num_thetae - 2);
        const std::size_t i1 = i0 + 1;
        const float weight_i =
            std::clamp(m_fi - static_cast<float>(i0), 0.0f, 1.0f);
        const auto [k0, k1, weight_k] = m_data->get_logp_indices(new_pres);

        return m_data->bilinear_interp(i0, i1, weight_i, k0, k1, weight_k);
    }

    /**
     * \brief Computes the virtual temperature of the parcel
     * \param   pres        Parcel pressure (Pa)
     * \param   tmpk        Parcel temperature (K)
     *
     * \return  The virtual temperature of the parcel
     */
    [[nodiscard]] inline float parcel_virtual_temperature(
        const float pres, const float tmpk) const {
        if (m_use_lifter) {
            return m_lifter.parcel_virtual_temperature(pres, tmpk);
        } else {
            return sharp::virtual_temperature(tmpk,
                                              sharp::mixratio(pres, tmpk));
        }
    }
};

//
////////////  END FUNCTORS   ///////////

/**
 * \brief Enum that defines the lifted parcel level (LPL) of origin.
 *
 * The SFC parcel is a surface-based parcel, where the parcel initial attributes
 * are set to the surface pressure, temperature, and dewpoint.
 *
 * The FCST parcel is a forecast-surface-based parcel, in which the afternoon
 * surface temperature and dewpoint are estimated and set as the parcel starting
 * values.
 *
 * The MU parcel is the most unstable parcel, in which the parcel attributes are
 * set to the pressure, temperature, and dewpoint of the maximum Theta-E level
 * within the bottom 400 hPa of the profile.
 *
 * The ML parcel is the mixed-layer parcel, in which the mean theta and water
 * vapor mixing ratio within the lowest 100 hPa are used to estimate a boundary
 * layer mean, and lifted from the surface.
 *
 * The USR parcel means that the parcel initial lifting attributes have already
 * been set by the programmer or user, and there is no need for them to be
 * set or modified.
 */
enum class LPL : int {
    /**
     * \brief Surface Based Parcel
     */
    SFC = 1,

    /**
     * \brief Forecast Surface Parcel
     */
    FCST = 2,

    /**
     * \brief Most Unstable Parcel
     */
    MU = 3,  // most unstable

    /**
     * \brief Mixed Layer Parcel
     */
    ML = 4,  // Mixed layer

    /**
     * \brief User-defined Parcel
     */
    USR = 5,  // user-defined
    END,
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Data that defines a Parcel, its attributes, and derived quantities.
 *
 * Contains information about a Parcel's starting level and
 * thermodynamic attributes, as well as paramaters computed
 * using the parcel.
 */
struct Parcel {
    /**
     * \brief Parcel starting pressure (Pa)
     */
    float pres = MISSING;

    /**
     * \brief Parcel starting temperature (K)
     */
    float tmpk = MISSING;

    /**
     * \brief Parcel starting dewpoint (K)
     */
    float dwpk = MISSING;

    /**
     * \brief Pressure at the Lifted Condensation Level (Pa)
     */
    float lcl_pressure = MISSING;

    /**
     * \brief Pressure at the Level of Free Convection (Pa)
     */
    float lfc_pressure = MISSING;

    /**
     * \brief Pressure at the parcel Equilibrium Level (Pa)
     */
    float eql_pressure = MISSING;

    /**
     * \brief Pressure at the Maximum Parcel Level (MPL) (Pa)
     */
    float mpl_pressure = MISSING;

    /**
     * \brief Parcel Convective Available Potential Energy (J/kg) between the
     * LFC and EL
     */
    float cape = 0.0;

    /**
     * \brief Parcel Convective Inhibition (J/kg) between the LFC and EL
     */
    float cinh = std::nanf("");

    /**
     * \brief The type of parcel this is
     */
    LPL source;

    Parcel();
    Parcel(const float pressure, const float temperature, const float dewpoint,
           const LPL lpl);

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Lifts a sharp::Parcel to compute its virtual temperature
     *
     * Lifts a sharp::Parcel dry adiabatically from its sharp::LPL to its
     * LCL dry adiabatically, and then moist adiabatically from the
     * LCL to the top of the profile. The moist adiabat used is determined
     * bu the type of lifting functor passed to the function (i.e.
     * sharp::lifter_wobus or sharp::lifter_cm1).
     *
     * \param   liftpcl         Parcel lifting function/functor
     * \param   pressure_arr    Array of env pressure (Pa)
     * \param   pcl_vtmpk_arr   The array to fill with virtual temp (K)
     * \param   N               The length of the arrays
     */
    template <typename Lft>
    void lift_parcel(Lft& liftpcl, const float pressure_arr[],
                     float pcl_vtmpk_arr[], const std::ptrdiff_t N) {
        // Lift the parcel from the LPL to the LCL
        float pres_lcl;
        float tmpk_lcl;
        drylift(this->pres, this->tmpk, this->dwpk, pres_lcl, tmpk_lcl);
        // If we are lifting elevated parcel (i.e. EIL), we need to make
        // sure our LCL isnt above the top of our data.
        if (pres_lcl < pressure_arr[N - 1]) return;

        this->lcl_pressure = pres_lcl;

        const float qv_lcl = mixratio(pres_lcl, tmpk_lcl);
        const float thetav_lcl =
            theta(pres_lcl, virtual_temperature(tmpk_lcl, qv_lcl),
                  THETA_REF_PRESSURE);

        // Define the dry and saturated lift layers
        PressureLayer dry_lyr = {this->pres, this->lcl_pressure};
        PressureLayer sat_lyr = {this->lcl_pressure, pressure_arr[N - 1]};

        // The LayerIndex excludes the top and bottom for interpolation reasons
        const LayerIndex dry_idx = get_layer_index(dry_lyr, pressure_arr, N);
        const LayerIndex sat_idx = get_layer_index(sat_lyr, pressure_arr, N);

        for (std::ptrdiff_t k = 0; k < dry_idx.kbot; ++k) {
            pcl_vtmpk_arr[k] = MISSING;
        }

        // Virtual potential temperature (Theta-V)
        // is conserved for a parcels dry ascent to the LCL
        for (std::ptrdiff_t k = dry_idx.kbot; k < dry_idx.ktop + 1; ++k) {
            const float pcl_vtmp =
                theta(THETA_REF_PRESSURE, thetav_lcl, pressure_arr[k]);
            pcl_vtmpk_arr[k] = pcl_vtmp;
        }

        float pres_bot = pres_lcl;
        float tmpk_bot = tmpk_lcl;
        liftpcl.setup(pres_lcl, tmpk_lcl);

        // fill the array with the moist parcel buoyancy
        for (std::ptrdiff_t k = sat_idx.kbot; k < N; ++k) {
            // compute above-lcl buoyancy here
            const float pcl_pres = pressure_arr[k];
            const float pcl_tmpk = liftpcl(pres_bot, tmpk_bot, pcl_pres);

            if constexpr (!Lft::lift_from_lcl) {
                pres_bot = pcl_pres;
                tmpk_bot = pcl_tmpk;
            }

            const float pcl_vtmpk =
                liftpcl.parcel_virtual_temperature(pcl_pres, pcl_tmpk);
            pcl_vtmpk_arr[k] = pcl_vtmpk;
        }
    }

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Find the LFC and EL that bounds the layer with the maximum CAPE
     *
     * Searches the buoyancy array for the LFC and EL combination that
     * results in the most CAPE in the given profile. The buoyancy array is
     * typically computed by calling sharp::Parcel::lift_parcel. Once the LFC
     * and EL are found, the values are set in sharp::Parcel::lfc_pres and
     * sharp::Parcel::eql_pres.
     *
     * The value of sharp::Parcel::eql_pres is sharp::MISSING if there there
     * is no qualifying level found within the data bounds (e.g. incomplete
     * data, or an EL above the available data). Any calls to
     * sharp::Parcel::cape_cinh will still compute CAPE without the presence of
     * an EL, using the best-available data.
     *
     * \param   pres_arr    The pressure coordinate array (Pa)
     * \param   hght_arr    The height coordinate array (meters)
     * \param   buoy_arr    The profile buoyancy array (m/s^2)
     * \param   N           The length of the arrays
     */
    void find_lfc_el(const float pres_arr[], const float hght_arr[],
                     const float buoy_arr[], const std::ptrdiff_t N);

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Find the pressure of the Maximum Parcel Level (MPL).
     *
     * The Maximum Parcel Level (MPL) is the level a parcel would reach
     * if it expended all of its integrated positive buoyancy past the
     * Eqilibrium Level. It is found by integrating negatively buoyant
     * area above the Equilibrium Level until the integrated negative
     * buoyancy is equal in magnitude to the Convective Available
     * Potential Energy betwee the Level of Free Convection and the
     * Equilibrium Level.
     *
     * For valid calculations, sharp::Parcel::cape_cinh must be
     * called first, or sharp::Parcel::cape and
     * sharp::Parcel::eql_pressure must be set.
     *
     * A value of sharp::MISSING is returned if:
     * - CAPE is 0
     * - sharp::Parcel::eql_pressure is sharp::MISSING
     * - No valid MPL candidate is found within the profile.
     *     - In this scenario, it likely exceeds the top of the available data.
     *
     * In addition to being returned, the result is stored inside of
     * sharp::Parcel::mpl_pressure.
     *
     * \param   pres_arr    Array of pressure   (Pa)
     * \param   hght_arr    Array of height     (meters)
     * \param   buoy_arr    Array of buoyancy   (m/s^2)
     * \param   N           Lengh of arrays
     *
     * \return  Maximum Parcel Level (MPL) Pressure (Pa)
     */
    float maximum_parcel_level(const float pres_arr[], const float hght_arr[],
                               const float buoy_arr[], const std::ptrdiff_t N);

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Compute CAPE and CINH for a previously lifted sharp::Parcel.
     *
     * Assuming that sharp::Parcel::lift_parcel has been called, cape_cinh
     * will integrate the area between the LFC and EL to compute CAPE,
     * and integrate the area between the LPL and LCL to compute CINH.
     *
     * If sharp::Parcel::eql_pressure is sharp::MISSING, but the
     * sharp::Parcel::lfc_pressure is defined, the routine will
     * compute CAPE with the available data despite the lack of
     * a defined equilibrium level. This is useful for incomplete
     * profile data, or pressure-level data where the EL is above
     * the top pressure value.
     *
     * The results are set in pcl->cape and pcl->cinh.
     *
     * \param   pres_arr    Array of pressure (Pa)
     * \param   hght_arr    Array of height (meters)
     * \param   buoy_arr    Array of buoyancy (m/s^2)
     * \param   N           Length of arrays
     */
    void cape_cinh(const float pres_arr[], const float hght_arr[],
                   const float buoy_arr[], const std::ptrdiff_t N);

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Computes the Lifted Index for the sharp::Parcel
     *
     * Assuming that sharp::Parcel::lift_parcel has been called,
     * this routine will compute the lifted index at the requested
     * pressure level.
     *
     * The lifted index is the difference between the environment
     * virtual temperature and the parcel virtual temperature.
     *
     * \param   pres_lev        Pressure level for LI (Pa)
     * \param   pres_arr        Array of pressure (Pa)
     * \param   vtmpk_arr       Array of environment virtual temperature (K)
     * \param   pcl_vtmpk_arr   Array of parcel virtual temperature (K)
     * \param   N               Length of arrays
     *
     * \return  Lifted Index (K)
     */
    float lifted_index(const float pres_lev, const float pres_arr[],
                       const float vtmpk_arr[], const float pcl_vtmpk_arr[],
                       const std::ptrdiff_t N);

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Construct and return a surface-based Parcel
     *
     * Given input values of surface pressure, temperature,
     * and dewpoint temperature, construct and return a Surface
     * Based Parcel.
     *
     * \param    pressure        Surface pressure (Pa)
     * \param    temperature     Surface temperature (K)
     * \param    dewpoint        Surface dewpoint (K)
     * \return   sharp::Parcel with surface values
     */
    static inline Parcel surface_parcel(const float pressure,
                                        const float temperature,
                                        const float dewpoint) noexcept {
        return Parcel(pressure, temperature, dewpoint, LPL::SFC);
    }

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Construct and return a mixed-layer Parcel
     *
     * Given input arrays of pressure, height, potential
     * temperature, and water vapor mixing ratio, as well
     * as a defined sharp::PressureLayer or
     * sharp::HeightLayer, compute and return a mixed-layer
     * Parcel.
     *
     * \param   mix_layer           sharp::PressureLayer or sharp::HeightLayer
     * \param   pressure            Array of pressure (Pa)
     * \param   height              Array of height (meters)
     * \param   pot_temperature     Array of potential temperature (K)
     * \param   wv_mixratio         Array of water vapor mixing ratio (unitless)
     * \param   N                   Length of arrays
     * \return sharp::Parcel with mixed-layer values
     */
    template <typename Lyr>
    static Parcel mixed_layer_parcel(Lyr& mix_layer, const float pressure[],
                                     [[maybe_unused]] const float height[],
                                     const float pot_temperature[],
                                     const float wv_mixratio[],
                                     const std::ptrdiff_t N) {
        float mean_mixr, mean_thta, pcl_pres;
        if constexpr (Lyr::coord == LayerCoordinate::pressure) {
            mean_mixr = layer_mean(mix_layer, pressure, wv_mixratio, N);
            mean_thta = layer_mean(mix_layer, pressure, pot_temperature, N);
            pcl_pres = mix_layer.bottom;

        } else {
            mean_mixr = layer_mean(mix_layer, height, pressure, wv_mixratio, N);
            mean_thta =
                layer_mean(mix_layer, height, pressure, pot_temperature, N);
            pcl_pres =
                sharp::interp_height(mix_layer.bottom, height, pressure, N);
        }
        const float tmpk = theta(THETA_REF_PRESSURE, mean_thta, pcl_pres);
        const float dwpk = temperature_at_mixratio(mean_mixr, pcl_pres);

        return Parcel(pcl_pres, tmpk, dwpk, LPL::ML);
    }

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Construct and return a most-unstable Parcel
     *
     * Given input arrays of pressure, height, temperature,
     * virtual temperature, and dewpoint, a scratch/working
     * array to store values of buoyancy, and a defined
     * sharp::PressureLayer or sharp::HeightLayer to search
     * over, find and return the most-unstable parcel.
     *
     * \param   pressure            Array of pressure (Pa)
     * \param   height              Array of height (meters)
     * \param   temperature         Array of temperature (K)
     * \param   virtemp             Array of virtual temperature (K)
     * \param   dewpoint            Array of dewpoint temperature (K)
     * \param   pcl_virtemp         Writeable array for parcel lifting calcs (K)
     * \param   buoy_arr            Writeable array for buoyancy calcs (m/s^2)
     * \param   N                   Length of arrays
     * \param   search_layer        sharp::PressureLayer or sharp::HeightLay
     * \param   lifter              The parcel moist adiabatic ascent function
     * \return The most unstable sharp::Parcel within the search layer
     */
    template <typename Lyr, typename Lft>
    static Parcel most_unstable_parcel(
        Lyr& search_layer, Lft& lifter, const float pressure[],
        const float height[], const float temperature[], const float virtemp[],
        const float dewpoint[], float pcl_virtemp[], float buoy_arr[],
        const std::ptrdiff_t N) {
        LayerIndex lyr_idx;
        if constexpr (Lyr::coord == LayerCoordinate::pressure) {
            lyr_idx = get_layer_index(search_layer, pressure, N);

        } else {
            lyr_idx = get_layer_index(search_layer, height, N);
        }

        Parcel max_parcel;
        for (std::ptrdiff_t idx = lyr_idx.kbot; idx <= lyr_idx.ktop; ++idx) {
            const float pres = pressure[idx];
            const float tmpk = temperature[idx];
            const float dwpk = dewpoint[idx];
            Parcel pcl(pres, tmpk, dwpk, LPL::MU);

            pcl.lift_parcel(lifter, pressure, pcl_virtemp, N);
            buoyancy(pcl_virtemp, virtemp, buoy_arr, N);
            pcl.cape_cinh(pressure, height, buoy_arr, N);
            if (pcl.cape > max_parcel.cape) max_parcel = pcl;
        }

        return max_parcel;
    }
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Data that defines a DowndraftParcel, its attributes, and derived
 * quantities.
 *
 * Contains information about a DowndraftParcel's starting level and
 * thermodynamic attributes, as well as paramaters computed
 * using the parcel.
 */
struct DowndraftParcel {
    /**
     * \brief DowndraftParcel starting pressure (Pa)
     */
    float pres = MISSING;

    /**
     * \brief DowndraftParcel starting temperature (K)
     */
    float tmpk = MISSING;

    /**
     * \brief DowndraftParcel starting dewpoint (K)
     */
    float dwpk = MISSING;

    /**
     * \brief DowndraftParcel Convective Available Potential Energy (J/kg)
     * between the LFC and EL
     */
    float cape = 0.0;

    /**
     * \brief DowndraftParcel Convective Inhibition (J/kg) between the LFC and
     * EL
     */
    float cinh = std::nanf("");

    DowndraftParcel();
    DowndraftParcel(const float pressure, const float temperature,
                    const float dewpoint);

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Lowers a saturated sharp::Parcel
     *
     * Lowers a saturated sharp::Parcel moist adiabatically from its
     * sharp::LPL to the surface. The moist adiabat used is determined
     * bu the type of lifting functor passed to the function (i.e.
     * sharp::lifter_wobus or sharp::lifter_cm1).
     *
     * Unlike sharp::lift_parcel, the virtual temperature correction
     * is not used for downdraft parcels.
     *
     * \param   liftpcl         Parcel lifting function/functor
     * \param   pressure_arr    Array of env pressure (Pa)
     * \param   pcl_tmpk_arr    The array to fill with parcel temperature
     * (K)
     * \param   N               The length of the arrays
     */
    template <typename Lft>
    void lower_parcel(Lft& liftpcl, const float pressure_arr[],
                      float pcl_tmpk_arr[], const std::ptrdiff_t N) {
        PressureLayer downdraft_layer = {pressure_arr[0], this->pres};
        LayerIndex downdraft_idx =
            get_layer_index(downdraft_layer, pressure_arr, N);

        // temperature is MISSING outside of the downdraft layer
        for (std::ptrdiff_t k = N - 1; k > downdraft_idx.ktop; --k) {
            pcl_tmpk_arr[k] = MISSING;
        }

        float wbt_pcl = wetbulb(liftpcl, this->pres, this->tmpk, this->dwpk);
        liftpcl.setup(this->pres, wbt_pcl);

        pcl_tmpk_arr[downdraft_idx.ktop] =
            liftpcl(this->pres, wbt_pcl, pressure_arr[downdraft_idx.ktop]);

        for (std::ptrdiff_t k = downdraft_idx.ktop - 1; k >= downdraft_idx.kbot;
             --k) {
            pcl_tmpk_arr[k] = liftpcl(pressure_arr[k + 1], pcl_tmpk_arr[k + 1],
                                      pressure_arr[k]);
        }
    }

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Compute CAPE and CINH for a previously lowered
     * sharp::DowndraftParcel.
     *
     * Assuming that sharp::DowndraftParcel::lower_parcel has been called,
     * cape_cinh will integrate the area between the LPL and the surface to
     * compute downdraft CAPE and downdraft CINH.
     *
     * The results are set in sharp::DowndraftParcel::cape and
     * sharp::DowndraftParcel::cinh.
     *
     * \param   pres_arr    Array of pressure   (Pa)
     * \param   hght_arr    Array of height     (meters)
     * \param   buoy_arr    Array of buoyancy   (m/s^2)
     * \param   N           Length of arrays
     */
    void cape_cinh(const float pres_arr[], const float hght_arr[],
                   const float buoy_arr[], const ptrdiff_t N);

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Define a downdraft parcel.
     *
     * Defines a downdraft parcel within a given search layer.
     * The downdraft parcel is defined as the mimumim layer-mean
     * equivalent potential temperature (Theta-E) within the
     * search layer. Typical values are to search within the lowest
     * 400 hPa of the profile, and a mean depth of 100 hPa.
     *
     * \param   search_layer    sharp::PressureLayer to search
     * \param   pressure        (Pa)
     * \param   temperature     (K)
     * \param   dewpoint        (K)
     * \param   thetae          (K)
     * \param   N               (length of arrays)
     * \param   mean_depth      (Pa)
     *
     * \return  the sharp::DowndraftParcel defining a downdraft parcel.
     */
    static DowndraftParcel min_thetae(
        sharp::PressureLayer& search_layer, const float pressure[],
        const float temperature[], const float dewpoint[], const float thetae[],
        const std::ptrdiff_t N, const float mean_depth = 10000.0f) {
        const sharp::LayerIndex lyr_idx =
            sharp::get_layer_index(search_layer, pressure, N);

        float min_thetae = 99999.0;
        float pres_of_min = sharp::MISSING;
        float half_depth = mean_depth / 2.0f;

        for (std::ptrdiff_t k = lyr_idx.ktop; k >= lyr_idx.kbot; --k) {
            sharp::PressureLayer mn_lyr = {pressure[k] + half_depth,
                                           pressure[k] - half_depth};
            float mean_thetae = sharp::layer_mean(mn_lyr, pressure, thetae, N);
            if (mean_thetae < min_thetae) {
                min_thetae = mean_thetae;
                pres_of_min = pressure[k];
            }
        }

        float pcl_t =
            sharp::interp_pressure(pres_of_min, pressure, temperature, N);
        float pcl_td =
            sharp::interp_pressure(pres_of_min, pressure, dewpoint, N);

        return DowndraftParcel(pres_of_min, pcl_t, pcl_td);
    }
};

/// @cond DOXYGEN_IGNORE

extern template Parcel
Parcel::most_unstable_parcel<PressureLayer, lifter_wobus>(
    PressureLayer& search_layer, lifter_wobus& lifter, const float pressure[],
    const float height[], const float temperature[], const float virtemp[],
    const float dewpoint[], float pcl_virtemp[], float buoy_arr[],
    const std::ptrdiff_t N);

extern template Parcel Parcel::most_unstable_parcel<HeightLayer, lifter_wobus>(
    HeightLayer& search_layer, lifter_wobus& lifter, const float pressure[],
    const float height[], const float temperature[], const float virtemp[],
    const float dewpoint[], float pcl_virtemp[], float buoy_arr[],
    const std::ptrdiff_t N);

extern template Parcel Parcel::most_unstable_parcel<PressureLayer, lifter_cm1>(
    PressureLayer& search_layer, lifter_cm1& lifter, const float pressure[],
    const float height[], const float temperature[], const float virtemp[],
    const float dewpoint[], float pcl_virtemp[], float buoy_arr[],
    const std::ptrdiff_t N);

extern template Parcel Parcel::most_unstable_parcel<HeightLayer, lifter_cm1>(
    HeightLayer& search_layer, lifter_cm1& lifter, const float pressure[],
    const float height[], const float temperature[], const float virtemp[],
    const float dewpoint[], float pcl_virtemp[], float buoy_arr[],
    const std::ptrdiff_t N);

extern template Parcel
Parcel::most_unstable_parcel<PressureLayer, lifter_lut<lifter_wobus>>(
    PressureLayer& search_layer, lifter_lut<lifter_wobus>& lifter,
    const float pressure[], const float height[], const float temperature[],
    const float virtemp[], const float dewpoint[], float pcl_virtemp[],
    float buoy_arr[], const std::ptrdiff_t N);

extern template Parcel
Parcel::most_unstable_parcel<HeightLayer, lifter_lut<lifter_wobus>>(
    HeightLayer& search_layer, lifter_lut<lifter_wobus>& lifter,
    const float pressure[], const float height[], const float temperature[],
    const float virtemp[], const float dewpoint[], float pcl_virtemp[],
    float buoy_arr[], const std::ptrdiff_t N);

extern template Parcel
Parcel::most_unstable_parcel<PressureLayer, lifter_lut<lifter_cm1>>(
    PressureLayer& search_layer, lifter_lut<lifter_cm1>& lifter,
    const float pressure[], const float height[], const float temperature[],
    const float virtemp[], const float dewpoint[], float pcl_virtemp[],
    float buoy_arr[], const std::ptrdiff_t N);

extern template Parcel
Parcel::most_unstable_parcel<HeightLayer, lifter_lut<lifter_cm1>>(
    HeightLayer& search_layer, lifter_lut<lifter_cm1>& lifter,
    const float pressure[], const float height[], const float temperature[],
    const float virtemp[], const float dewpoint[], float pcl_virtemp[],
    float buoy_arr[], const std::ptrdiff_t N);

extern template Parcel Parcel::mixed_layer_parcel<PressureLayer>(
    PressureLayer& mix_layer, const float pressure[], const float height[],
    const float pot_temperature[], const float wv_mixratio[],
    const std::ptrdiff_t N);

extern template Parcel Parcel::mixed_layer_parcel<HeightLayer>(
    HeightLayer& mix_layer, const float pressure[], const float height[],
    const float pot_temperature[], const float wv_mixratio[],
    const std::ptrdiff_t N);

extern template void Parcel::lift_parcel<lifter_wobus>(
    lifter_wobus& liftpcl, const float pressure_arr[], float pcl_vtmpk_arr[],
    const std::ptrdiff_t N);

extern template void Parcel::lift_parcel<lifter_cm1>(lifter_cm1& liftpcl,
                                                     const float pressure_arr[],
                                                     float pcl_vtmpk_arr[],
                                                     const std::ptrdiff_t N);

extern template void Parcel::lift_parcel<lifter_lut<lifter_wobus>>(
    lifter_lut<lifter_wobus>& liftpcl, const float pressure_arr[],
    float pcl_vtmpk_arr[], const std::ptrdiff_t N);

extern template void Parcel::lift_parcel<lifter_lut<lifter_cm1>>(
    lifter_lut<lifter_cm1>& liftpcl, const float pressure_arr[],
    float pcl_vtmpk_arr[], const std::ptrdiff_t N);

extern template void DowndraftParcel::lower_parcel<lifter_wobus>(
    lifter_wobus& liftpcl, const float pressure_arr[], float pcl_tmpk_arr[],
    const std::ptrdiff_t N);

extern template void DowndraftParcel::lower_parcel<lifter_cm1>(
    lifter_cm1& liftpcl, const float pressure_arr[], float pcl_tmpk_arr[],
    const std::ptrdiff_t N);

extern template void DowndraftParcel::lower_parcel<lifter_lut<lifter_wobus>>(
    lifter_lut<lifter_wobus>& liftpcl, const float pressure_arr[],
    float pcl_tmpk_arr[], const std::ptrdiff_t N);

extern template void DowndraftParcel::lower_parcel<lifter_lut<lifter_cm1>>(
    lifter_lut<lifter_cm1>& liftpcl, const float pressure_arr[],
    float pcl_tmpk_arr[], const std::ptrdiff_t N);

/// @endcond

}  // end namespace sharp

#endif  // SHARP_PARCEL_H
