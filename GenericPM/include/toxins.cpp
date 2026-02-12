/**
 * @file toxins.cpp
 * @brief Toxins class implementation for DON (mycotoxin) production modeling
 * 
 * @author Thiago Berton Ferreira
 * @author Willingthon Pavan (wpavan.us@gmail.com)
 * @author Jose Mauricio Cunha Fernandes (jmauricio.fernandes@icloud.com)
 * 
 * @copyright Copyright (c) 2017–2025, DSSAT Foundation
 * @license BSD-3-Clause. See the LICENSE file in the root folder for details.
 */
#include "toxins.h"
#include "../../FlexibleIO/Data/FlexibleIO.hpp"
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>

Toxins::FHB_DON Toxins::calculateDailyDON(
    int YRDOY,
    int RUN,
    int YRSIM,
    double ZSTAGE,
    double total_biomass,
    double total_damaged_tissue,
    size_t cohort_count,
    double WSDD,
    double SDWT,
    double dIdt,
    int heading_yrdoy_global,
    bool heading_detected_global
) {
    // Static variables to persist between calls
    static double cumulative_DON_total_ug = 0.0;
    static int last_run_don = -1;
    static int last_yrsim_don = -1;
    static int last_run_file = -1;
    static int last_yrsim_file = -1;
    
    FlexibleIO *fio = FlexibleIO::getInstance();
    
    // Reset cumulative DON when new run starts
    if (last_run_don != RUN) {
        cumulative_DON_total_ug = 0.0;
        last_run_don = RUN;
        last_yrsim_don = YRSIM;
    }
    
    // ═══════════════════════════════════════════════════════════════════════
    // DAILY DON PRODUCTION MODEL (Temperature and Weather-Dependent)
    // ═══════════════════════════════════════════════════════════════════════
    
    // Get weather data for current day
    double TAVG = fio->getReal("PEST", "TAVG");
    double TMIN = fio->getRealYrdoy("WTH", std::to_string(YRDOY), "TMIN");
    double TMAX = fio->getRealYrdoy("WTH", std::to_string(YRDOY), "TMAX");
    double RAIN = fio->getRealYrdoy("WTH", std::to_string(YRDOY), "RAIN");
    
    // Calculate days after heading (using global tracking variable)
    int days_after_heading = heading_detected_global ? (YRDOY - heading_yrdoy_global) : 0;
    
    // ═══════════════════════════════════════════════════════════════════════
    // TEMPERATURE-DEPENDENT BASE RATE
    // ═══════════════════════════════════════════════════════════════════════
    // Double sigmoid bell curve: optimal 15-25°C, max rate 75 µg/mg
    // Centers at 10°C (left) and 30°C (right), 1% minimum baseline
    double base_rate = 75.0 * (0.01 + 0.99 / ((1.0 + std::exp(-0.5 * (TAVG - 10.0))) * 
                                                (1.0 + std::exp(0.5 * (TAVG - 30.0)))));
    
    // ═══════════════════════════════════════════════════════════════════════
    // RELATIVE HUMIDITY FACTOR
    // ═══════════════════════════════════════════════════════════════════════
    
    // Spore type parameters
    std::string spore_type = "combined";  // Options: "macroconidia", "ascospores", "combined"
    double macro_proportion = 0.6;        // Weight for macroconidia (60%) vs ascospores (40%)
    
    double rh_factor = 1.0;
    double relative_humidity_percent = -1.0; // Set to -1 to indicate RH not available
    
    // Optional: Estimate RH from TMIN/TMAX (uncomment if you want to use this)
    // relative_humidity_percent = 100.0 * std::exp(0.06 * (TMIN - TMAX));
    
    if (relative_humidity_percent < 0) {
        // If RH not provided, assume optimal conditions
        rh_factor = 1.0;
    } else {
        // Constrain RH to valid range
        double rh = std::max(0.0, std::min(100.0, relative_humidity_percent));
        
        if (spore_type == "macroconidia") {
            // Rain-splash dispersed: optimal >90%, center 85%, range 0.1-1.0
            rh_factor = 0.1 + 0.9 / (1.0 + std::exp(-0.15 * (rh - 85.0)));
            
        } else if (spore_type == "ascospores") {
            // Wind-dispersed: optimal >85%, center 80%, range 0.2-1.0
            rh_factor = 0.2 + 0.8 / (1.0 + std::exp(-0.12 * (rh - 80.0)));
            
        } else {  // combined (default)
            // Weighted average of both spore types
            double macro_factor = 0.1 + 0.9 / (1.0 + std::exp(-0.15 * (rh - 85.0)));
            double asco_factor = 0.2 + 0.8 / (1.0 + std::exp(-0.12 * (rh - 80.0)));
            rh_factor = macro_proportion * macro_factor + (1.0 - macro_proportion) * asco_factor;
        }
    }
    
    // ═══════════════════════════════════════════════════════════════════════
    // WEATHER MODULATION FACTOR
    // ═══════════════════════════════════════════════════════════════════════
    // Combined weather: rain (0.5-1.2), cold penalty (<10°C), hot penalty (>32°C)
    double weather_factor = (0.5 + 0.7 / (1.0 + std::exp(-1.2 * (RAIN - 3.0)))) * 
                           (1.0 - 0.7 / (1.0 + std::exp(0.8 * (TMIN - 10.0)))) * 
                           (1.0 - 0.6 / (1.0 + std::exp(-0.8 * (TMAX - 32.0))));
    
    // ═══════════════════════════════════════════════════════════════════════
    // TEMPORAL FACTOR (Days after heading)
    // ═══════════════════════════════════════════════════════════════════════
    // Gompertz decay: optimal 0-10 DAH, asymptotic minimum 0.1 at 40+ DAH
    double temporal_factor = 0.1 + 0.9 * std::exp(-std::exp(-2.0 + 0.12 * days_after_heading));
    
    // ═══════════════════════════════════════════════════════════════════════
    // CALCULATE DAILY DON PRODUCTION
    // ═══════════════════════════════════════════════════════════════════════
    // Convert fungal biomass g/m² to mg/m²
    double fungal_biomass_mg_m2 = total_biomass * 1000.0;
    
    // DON = base_rate(T) × biomass × weather × temporal
    // Daily DON production in total mass (µg)
    double daily_DON_total_ug = base_rate * fungal_biomass_mg_m2 * weather_factor * 
                                temporal_factor; //rh_factor;
    
    // ═══════════════════════════════════════════════════════════════════════
    // DON ACCUMULATION STRATEGY
    // ═══════════════════════════════════════════════════════════════════════
    // We accumulate total DON MASS (µg) over the season, NOT concentration.
    // This is critical because:
    //   1. DON concentration = cumulative_DON_mass / current_grain_weight
    //   2. As grain fills (SDWT increases), concentration naturally dilutes
    //   3. Adding concentrations directly (wrong approach) causes mathematical artifacts
    //   4. This matches laboratory measurement: total toxin per total grain sample
    
    // Accumulate total DON mass (µg) produced daily
    cumulative_DON_total_ug += daily_DON_total_ug;
    
    // Calculate current concentration for reporting (µg/kg)
    double grain_yield_g_m2 = SDWT;
    double cumulative_DON_ug_kg = (grain_yield_g_m2 > 0) ? 
                                   (cumulative_DON_total_ug / grain_yield_g_m2) : 0.0;
    double daily_DON_ug_kg = (grain_yield_g_m2 > 0) ? 
                             (daily_DON_total_ug / grain_yield_g_m2) : 0.0;

    printf("YRDOY: %i, ZSTAGE: %.4f, dIdt: %.4f, HSDWT: %.4f, total_damaged_tissue: %.4f, cohorts.size(): %zu, WSDD: %.6f, total_biomass: %.6f, daily_DON: %.2f ug/kg, cumulative_DON: %.2f ug/kg\n",
           YRDOY, ZSTAGE, dIdt, SDWT - total_damaged_tissue, total_damaged_tissue,
           cohort_count, WSDD, total_biomass, daily_DON_ug_kg, cumulative_DON_ug_kg);
    
    // ═══════════════════════════════════════════════════════════════════════
    // WRITE DAILY DON OUTPUT FILE (grouped by RUN and YRSIM like simulation_results)
    // ═══════════════════════════════════════════════════════════════════════
    
    std::string YRSIM_str = std::to_string(YRSIM);
    std::string year = YRSIM_str.substr(0, YRSIM_str.size() - 3);
    std::string don_fileName = "DON_results_RUN" + std::to_string(RUN) + "_" + year + ".csv";
    
    // Check if new run or new simulation year started - create new file with header
    if (RUN != last_run_file || YRSIM != last_yrsim_file) {
        std::remove(don_fileName.c_str());
        std::ofstream fdon(don_fileName);
        fdon << "YRDOY,ZSTAGE,Infection_Rate,Fungal_Biomass_g_m2,Damaged_Tissue_g_m2,Cohorts,"
             << "WSDD_g_m2,Daily_DON_ug_kg,Cumulative_DON_ug_kg\n";
        fdon.close();
        last_run_file = RUN;
        last_yrsim_file = YRSIM;
    }
    
    // Append daily DON data
    std::ofstream fdon(don_fileName, std::ios::app);
    fdon << YRDOY << ',' 
         << ZSTAGE << ',' 
         << dIdt << ',' 
         << total_biomass << ',' 
         << total_damaged_tissue << ',' 
         << cohort_count << ',' 
         << WSDD << ',' 
         << daily_DON_ug_kg << ',' 
         << cumulative_DON_ug_kg << '\n';
    fdon.close();
    
    // Return results
    FHB_DON result;
    result.daily_DON_total_ug = daily_DON_total_ug;
    result.cumulative_DON_total_ug = cumulative_DON_total_ug;
    result.daily_DON_ug_kg = daily_DON_ug_kg;
    result.cumulative_DON_ug_kg = cumulative_DON_ug_kg;
    
    return result;
}

