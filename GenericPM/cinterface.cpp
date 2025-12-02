/**
 * @file cinterface.cpp
 * 
 * @author Willingthon Pavan (wpavan.us@gmail.com)
 * @author Jose Mauricio Cunha Fernandes (jmauricio.fernandes@icloud.com)
 * 
 * @copyright Copyright (c) 2017–2025, DSSAT Foundation
 * @license BSD-3-Clause. See the LICENSE file in the root folder for details.
 */
#include "include/simulator.h"
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <tuple>
#include <utility>
#include "include/basic.h"
#include "../FlexibleIO/Data/FlexibleIO.hpp"
using namespace std;

extern "C" {
    // Coupling Functions
    int couplingInit(int *YRDOY, int *YRPLT);
    int couplingRate(int *YRDOY,
            float *AREALF, float *CLW, float *CSW, float *PCLMT, float *PCSTMD,
            float *PDLA, float *PLFAD, float *PLFMD, float *PSTMD, float *PVSTGD,
            float *SLA, float *SLDOT, float *SSDOT, float *STMWT, float *TDLA,
            float *VSTGD, float *WLFDOT, float *WSTMD, float *WTLF,
            float *TLFAD, float *TLFMD, float *VSTAGE, float *WLIDOT,
            float *CLAI, float *CLFM, float *CSTEM, float *DISLA, float *DISLAP,
            float *LAIDOT, float *WSIDOT, float *SDWT, float *WSDD, 
            float *PSDD, int *DAS, int *YRPLT);
    int couplingIntegration(int *YRDOY,
            float *AREALF, float *CLW, float *CSW, float *PCLMT, float *PCSTMD,
            float *PDLA, float *PLFAD, float *PLFMD, float *PSTMD, float *PVSTGD,
            float *SLA, float *SLDOT, float *SSDOT, float *STMWT, float *TDLA,
            float *VSTGD, float *WLFDOT, float *WSTMD, float *WTLF,
            float *TLFAD, float *TLFMD, float *VSTAGE, float *WLIDOT,
            float *CLAI, float *CLFM, float *CSTEM, float *DISLA, float *DISLAP,
            float *LAIDOT, float *WSIDOT, float *SDWT, float *WSDD, 
            float *PSDD, int *DAS);
    int couplingOutput(int *doy);

}

float CLWp, SLAp, SDWTp, SW, SL1, SLL1, SSAT1, SDUL1;

// NEW **************************************************

// Structure to hold treatment summary data for DON output
struct TreatmentSummary {
    int last_YRDOY;
    double cumulative_DON;
    double total_biomass;
    double total_damaged_tissue;
    size_t cohort_count;
};

// Global static variables for DON summary tracking
static std::map<std::pair<int,int>, TreatmentSummary> treatment_summaries;
static std::map<std::pair<int,int>, bool> summary_written;
static bool header_written_summary = false;
static int last_run_summary = -1;
static int last_yrsim_summary = -1;

double kr;
double M_thresh;
double kd;
double kg;
double alpha;
double kdw;
double k_inf;
double k_rec;
double rain_threshold;
double time_step_size;
double S , W , I ;
using Vec3 = std::array<double,3>;
Vec3 y = {S, W, I};
double dSdt, dWdt, dIdt;
double daily_dI;           // daily increment in infection (for cohort biomass calculation)
double daily_dW;           // daily increment in wheat (for cohort biomass calculation)

// Global variables for DON model (shared between couplingRate and couplingIntegration)
int heading_yrdoy_global = -1;
bool heading_detected_global = false;
int first_day_sus_global = 0;

// NEW **************************************************


// Coupling Functions Implementation 

int couplingInit(int *YRDOY, int *YRPLT) {
    // Set the start day for Disease Model
    // Get an instance of Simulator
    Simulator *s = Simulator::newInstance();
    // Set the start day for Disease Model
    s->setCurrentYearDoy(*YRDOY);
    // Set the sowing/planting date
    s->getCropInterface()->setPlantingDate(*YRPLT);

        kr = 300;
        M_thresh = 0.65;
        kd = 0.1;
        kg = 0.05;
        alpha = 0.7;
        kdw = 0.08;
        k_inf = 0.0025;
        k_rec = 0.0;
        rain_threshold = 2.0;
        time_step_size = 1;
        S = 100, W = 0.0, I = 0.0;
        dSdt = S;
        dWdt = W;
        dIdt = I;
        daily_dI = 0.0;
        daily_dW = 0.0;
        y = {S, W, I};
        
    // Initialize DON model global variables
    heading_yrdoy_global = -1;
    heading_detected_global = false;
    first_day_sus_global = 0;
    
    return (1);
}

// NEW **************************************************

double temperature_beta(double T, double tmin, double topt, double tmax) {
    // Check for invalid or non-finite inputs
    if (!std::isfinite(T) || !std::isfinite(tmin) || !std::isfinite(topt) || !std::isfinite(tmax))
        return 0.0;

    // Degenerate parameter case
    if (tmin >= topt || topt >= tmax)
        return 1.0;

    if (T <= tmin || T >= tmax)
        return 0.0;

    // Avoid division by zero
    double denom1 = std::max(topt - tmin, 1e-9);
    double denom2 = std::max(tmax - topt, 1e-9);

    double a = denom1 / denom2;
    double b = denom2 / denom1;

    double x1 = std::max((T - tmin) / denom1, 1e-12);
    double x2 = std::max((tmax - T) / denom2, 1e-12);

    double f = std::pow(x1, a) * std::pow(x2, b);

    return std::min(std::max(f, 0.0), 1.0);
}

Vec3 ode_rhs(Vec3 &state, double Rain_t, double day, double SW, int is_rainy_day, double TAVG, double inf_temp_max, double inf_temp_min, double inf_temp_opt) {

    
    double A = 0.022900;
    double B = 3.612468;
    double C = 0.464022;
    double anther_prop_t = A * std::pow(day, B) * std::exp(-C * day);
    anther_prop_t = std::max(0.0, std::min(1.0, anther_prop_t));

    double spore_release_rate = kr * std::max(0.0, SW - M_thresh);

    double S = state[0];
    double W = state[1];
    double I = state[2];

    Rain_t = std::max(0.0, Rain_t);
    
    double temp_factor = temperature_beta(TAVG, inf_temp_min, inf_temp_opt, inf_temp_max);

    dSdt = std::max(0.0,spore_release_rate - (kd * S) - (kg * S));
    dWdt = std::max(0.0, (alpha * kg * S) - (kdw * W));
    dIdt = std::max(0.0,(k_inf * W * is_rainy_day * anther_prop_t * temp_factor) - (k_rec * I));

    // infection floor (prevent negative across dt)
    double dt = time_step_size;
    if ((I + dIdt * dt) < 0.0) dIdt = -I / dt;
    
    //printf("DEPOIS S %f W %f I %f SW %f \n", S, W , I, SW);
    
    //std::ofstream fout("simulation_results.csv", std::ios::app);
    //fout << spore_release_rate << ',' << S << ',' << W << ',' << I <<'\n';
    //fout.close();    
    

    return Vec3{dSdt, dWdt, dIdt};
}

Vec3 rk4_step(Vec3 &y, double h, double Rain_t, double day, double SW, int is_rainy_day, 
    double TAVG, double inf_temp_max, double inf_temp_min, double inf_temp_opt) {
    Vec3 k1 = ode_rhs(y, Rain_t, day, SW, is_rainy_day, TAVG, inf_temp_max, inf_temp_min, inf_temp_opt);
    Vec3 y2;
    for (int i = 0; i < 3; ++i) y2[i] = y[i] + 0.5*h*k1[i];
    Vec3 k2 = ode_rhs(y2, Rain_t, day, SW, is_rainy_day, TAVG, inf_temp_max, inf_temp_min, inf_temp_opt);
    Vec3 y3;
    for (int i = 0; i < 3; ++i) y3[i] = y[i] + 0.5*h*k2[i];
    Vec3 k3 = ode_rhs(y3, Rain_t, day, SW, is_rainy_day, TAVG, inf_temp_max, inf_temp_min, inf_temp_opt);
    Vec3 y4;
    for (int i = 0; i < 3; ++i) y4[i] = y[i] + h*k3[i];
    Vec3 k4 = ode_rhs(y4, Rain_t, day, SW, is_rainy_day, TAVG, inf_temp_max, inf_temp_min, inf_temp_opt);
    Vec3 ynew;
    for (int i = 0; i < 3; ++i) ynew[i] = y[i] + (h/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    // enforce non-negatives on some states
    if (ynew[0] < 0.0) ynew[0] = 0.0; // S
    if (ynew[1] < 0.0) ynew[1] = 0.0; // W
    if (ynew[2] < 0.0) ynew[2] = 0.0; // I
        
    return ynew;
}
// NEW **************************************************


int couplingRate(int *YRDOY,
        float *AREALF, float *CLW, float *CSW, float *PCLMT, float *PCSTMD,
        float *PDLA, float *PLFAD, float *PLFMD, float *PSTMD, float *PVSTGD,
        float *SLA, float *SLDOT, float *SSDOT, float *STMWT, float *TDLA,
        float *VSTGD, float *WLFDOT, float *WSTMD, float *WTLF,
        float *TLFAD, float *TLFMD, float *VSTAGE, float *WLIDOT,
        float *CLAI, float *CLFM, float *CSTEM, float *DISLA, float *DISLAP,
        float *LAIDOT, float *WSIDOT, float *SDWT, 
        float *WSDD, float *PSDD, int *DAS, int *YRPLT) {
    // Temporary variable used for computations 
    float temp = 0, newOrgan = 0, SW = 0, SL1 = 0, SLL1 = 0, SSAT1 = 0, SDUL1 = 0, TAVG = 0;
    float CloudField = 0;
    FlexibleIO *fio = FlexibleIO::getInstance();
    // Get an instance of Simulator
    Simulator *s = Simulator::getInstance();
    newOrgan = s->getCropInterface()->getOrgansQtd()+1;

    // Set the sowing/planting date
    if(s->getCropInterface()->getPlantingDate() < 0) {
        s->getCropInterface()->setPlantingDate(*YRPLT);
    }    
    // Set the current YearDOY for next Disease step computation
    s->updateCurrentYearDoy(*YRDOY);

    TAVG = fio->getReal("PEST", "TAVG");
    double inf_temp_max = 32.0;//fio->getRealIndex("PST", "TFS", 1);
    double inf_temp_min = 10.0; //fio->getRealIndex("PST", "TFS", 2);
    double inf_temp_opt = 23.0;//fio->getRealIndex("PST", "TFS", 3);
  
    
    //printf("inf_temp_max %f inf_temp_min %f inf_temp_opt %f\n");
    
    
    double ZSTAGE = fio->getReal("PEST", "ZSTAGE");
    
    fio->setIntegerMemory("PEST", "YRDOY", *YRDOY);
    
    SL1 = fio->getReal("PEST", "SL1");
    SLL1 = fio->getReal("PEST", "SLL1");
    SDUL1 = fio->getReal("PEST", "SDUL1");
    SSAT1 = fio->getReal("PEST", "SSAT1");

    SW = std::min(100.0f, std::max(0.0f, (SL1-SLL1)/(SSAT1-SLL1)));
    
    //printf("YRDOY %i SL1 %f SLL1 %f SDUL1 %f SSAT1 %f SW %f\n", *YRDOY, SL1, SLL1, SDUL1, SSAT1, SW);
    
    //NEW CODE
    //*****************************************************************
    
    int day = 0;
    static int SUSTAGE = 0;

    
    double Rain_t = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "RAIN");
    int dap = fio->getReal("PEST", "YRPLT");
    
    
    int is_rainy_day = (Rain_t > rain_threshold) ? 1 : 0;
    
    if(ZSTAGE >= 51 && SUSTAGE == 0){
        first_day_sus_global = *YRDOY;
        heading_yrdoy_global = *YRDOY;
        heading_detected_global = true;
        SUSTAGE = 1;
    }
    
    day = *YRDOY - first_day_sus_global;
    

    double A = 0.022900;
    double B = 3.612468;
    double C = 0.464022;
    double anther_prop_t = A * std::pow(day, B) * std::exp(-C * day);
    anther_prop_t = std::max(0.0, std::min(1.0, anther_prop_t));

    //double S = fio->getReal("PST", "II"); // Initial number of spores in the air
    //double W = 0; // Initial number of spores on wheat spikes
    //double I = 0; // Initial level of infection
    
    //double S = 0, W = 0.0, I = 0.0;


    //printf("YRDOY %i kr %f SW %f SL1 %f M_thresh %f\n", *YRDOY, kr, SW, SL1, M_thresh);

    //double spore_release_rate = kr * std::max(0.0, SW - M_thresh);
    

    // **************************************************

    // ---- RK4 integration ----
        
    ////double h = dap;
    double h = 1; // daily step

    Vec3 out_dis = rk4_step(y, h, Rain_t, day, SW, is_rainy_day, TAVG, inf_temp_max, inf_temp_min, inf_temp_opt);
    
    // Store old infection level before integration
    double I_old = y[2];
    double W_old = y[1];
    
    dSdt = out_dis[0];
    dWdt = out_dis[1];
    dIdt = out_dis[2];
    
    y = {dSdt, dWdt, dIdt};
    
    daily_dI = dIdt - I_old;
    daily_dW = dWdt - W_old;

    
    int YRSIM = fio->getReal("PEST", "YRSIM");
    int RUN = fio->getInteger("PEST", "RUN");
    
    double SRAD_t = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "SRAD");
    double TMAX_t = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "TMAX");
    double TMIN_t = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "TMIN");
    
    static int last_run_output = -1;
    static int last_yrsim_output = -1;

    std::string YRSIM_str = std::to_string(YRSIM);
    std::string year = YRSIM_str.substr(0, YRSIM_str.size() - 3);
    std::string fileName = "simulation_results_RUN" + std::to_string(RUN) + "_" + year + ".csv";

    // Check if new run or new simulation year started
    if (RUN != last_run_output || YRSIM != last_yrsim_output) {
        // Remove any previous results for this run/year combination
        std::remove(fileName.c_str());
        
        // Initialize new results file with header
        std::ofstream fout(fileName);
        fout << "YRDOY,Spores_in_Air,Spores_on_Wheat,Infection,Phen_Stage,Soil_Water(%),SRAD,TMAX,TMIN,RAIN\n";
        fout.close();

        // Reset disease susceptibility stage tracker
        SUSTAGE = 0;
        last_run_output = RUN;
        last_yrsim_output = YRSIM;
    }
    
    //printf("YRSIM %i YRSIMp %i SUSTAGE %i \n", YRSIM, YRSIMp, SUSTAGE);

    // Append new row each day
    // Also print to terminal for debugging/traceability
    //printf("CSV YRDOY=%d, dSdt=%.6f, dWdt=%.6f, dIdt=%.6f, ZSTAGE=%.3f, SW=%.4f, SRAD=%.3f, TMAX=%.2f, TMIN=%.2f, RAIN=%.3f\n",
    //    *YRDOY, dSdt, dWdt, dIdt, ZSTAGE, SW, SRAD_t, TMAX_t, TMIN_t, Rain_t);
    std::ofstream fout(fileName, std::ios::app);
    fout << *YRDOY << ',' << dSdt << ',' << dWdt << ',' << dIdt << ',' << ZSTAGE << ',' << SW << ',' << SRAD_t << ',' << TMAX_t << ',' << TMIN_t << ',' << Rain_t <<'\n';
    fout.close();
    
    // NEW **************************************************
    return (1);
}

// Structure to represent a disease cohort
struct DiseaseCohort {
    int infection_day;        // YRDOY when cohort was created
    int age;                  // days since infection
    double biomass;           // fungal biomass in this cohort (g)
    double damaged_tissue;    // amount of grain tissue damaged by this cohort (g)
    std::string stage;        // "latent", "infectious", or "necrotic"
};

int couplingIntegration(int *YRDOY,
    float *AREALF, float *CLW, float *CSW, float *PCLMT, float *PCSTMD,
    float *PDLA, float *PLFAD, float *PLFMD, float *PSTMD, float *PVSTGD,
    float *SLA, float *SLDOT, float *SSDOT, float *STMWT, float *TDLA,
    float *VSTGD, float *WLFDOT, float *WSTMD, float *WTLF,
    float *TLFAD, float *TLFMD, float *VSTAGE, float *WLIDOT,
    float *CLAI, float *CLFM, float *CSTEM, float *DISLA, float *DISLAP,
    float *LAIDOT, float *WSIDOT, float *SDWT, 
    float *WSDD, float *PSDD, int *DAS) {
      
    // Fungal growth model parameters
    double Y = 0.4;                // yield: fraction of consumed substrate converted to biomass
    double r_max = 0.3;            // intrinsic fungal growth rate (1/day)
    double t_lag = 5;              // lag duration (days)
    double lag_slope = 1.0;        // controls sharpness of lag activation
    
    // Cohort life cycle parameters (days)
    // int latent_period = 5;         // days before lesion becomes infectious
    // int necrotic_start = 20;       // day when lesion becomes necrotic
    
    // Initial infection parameters
    double B0 = 0.001;             // initial fungal biomass per new infection (g)
    
    // Static variables to persist between calls
    static std::vector<DiseaseCohort> cohorts; //Cohorts -- check if this sort of implementation is ok
    static double total_damaged_tissue = 0.0;  // cumulative damaged grain tissue (g)
    static int last_YRDOY = -1;
    static int last_run = -1;
    
    // Reset cohorts when new run starts (critical for sequential runs in batch mode)
    FlexibleIO *fio = FlexibleIO::getInstance();
    int YRSIM = fio->getReal("PEST", "YRSIM");
    int RUN = fio->getInteger("PEST", "RUN");
    double ZSTAGE = fio->getReal("PEST", "ZSTAGE");
    
    if (last_run != RUN) {
        cohorts.clear();
        total_damaged_tissue = 0.0;
        last_run = RUN;
        last_YRDOY = -1;
    }
        
    // SDWT - current grain weight (substrate - g/m²)
    // HSDWT - healthy grain weight
    double HSDWT = std::max(0.0, *SDWT - total_damaged_tissue);
    
    // Only process if there's available substrate
    if (HSDWT > 0) {
        
        //printf("AAAAAABBCCCCCCC YRDOY=%d, dSdt=%.6f, dWdt=%.6f, dIdt=%.6f, ZSTAGE=%.3f, SW=%.4f\n",
        //*YRDOY, dSdt, dWdt, dIdt, ZSTAGE, SW);
        
        //printf("(dIdt*dWdt): %.6f\n", (dIdt*dWdt));
        //printf("(daily_dI*daily_dW): %.6f\n", (daily_dI*daily_dW));
        
        if (dIdt > 0.0) {
            DiseaseCohort new_cohort;
            new_cohort.infection_day = *YRDOY;
            new_cohort.age = 0;
            //new_cohort.biomass = B0; // this should be B0 * (dIdt*dWdt)?
            //new_cohort.biomass = B0 * dIdt; 
            //new_cohort.biomass = B0 * (dIdt*dWdt); 
            //new_cohort.biomass = B0 * (daily_dI*daily_dW);
            new_cohort.biomass = B0 * (daily_dI);
            new_cohort.damaged_tissue = 0.0;
            cohorts.push_back(new_cohort);
            
            //printf("  [NEW COHORT] Day %i - daily_dI=%.6f, initial_biomass=%.8f (Total cohorts: %zu)\n", 
            //       *YRDOY, daily_dI, new_cohort.biomass, cohorts.size());
        }
        
        double total_cohort_biomass = 0.0;
        double daily_new_damage = 0.0;
        
        for (auto& cohort : cohorts) {
            // Update cohort age
            cohort.age = (*YRDOY - cohort.infection_day);
            
            double activation = 1.0 / (1.0 + std::exp(-lag_slope * (cohort.age - t_lag)));
            double r_eff = r_max * activation;
            
            // Calculate total biomass across all cohorts for competition
            total_cohort_biomass = 0.0;
            for (const auto& c : cohorts) {
                total_cohort_biomass += c.biomass;
            }
            
            // Growth is proportional to available healthy tissue
            double growth_limit = 1.0 - (total_cohort_biomass / (Y * HSDWT));
            growth_limit = std::max(0.0, std::min(1.0, growth_limit));
            
            double dB = r_eff * cohort.biomass * growth_limit;
            cohort.biomass += dB;
            
            // Calculate damage (substrate consumption)
            // Damaged tissue = fungal biomass / yield coefficient
            double tissue_consumed = (1.0 / Y) * dB;
            cohort.damaged_tissue += tissue_consumed;
            daily_new_damage += tissue_consumed;
                
        }
        
        total_damaged_tissue += daily_new_damage;
        // *WSDD is the actual mass damaged (g/m²)
        *WSDD = total_damaged_tissue; 
        
        if (total_damaged_tissue > HSDWT) {
            total_damaged_tissue = HSDWT;
        }
        
    } 

    // Compute total fungal biomass across all cohorts for reporting
    double total_biomass = 0.0;
    for (const auto& c : cohorts) {
        total_biomass += c.biomass;
    }
    
    // ═══════════════════════════════════════════════════════════════════════
    // DAILY DON PRODUCTION MODEL (Temperature and Weather-Dependent)
    // ═══════════════════════════════════════════════════════════════════════
    
    // Get weather data for current day
    double TAVG = fio->getReal("PEST", "TAVG");
    double TMIN = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "TMIN");
    double TMAX = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "TMAX");
    double RAIN = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "RAIN");
    
    // Calculate days after heading (using global tracking variable)
    int days_after_heading = heading_detected_global ? (*YRDOY - heading_yrdoy_global) : 0;
    
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
                           (1.0 - 0.7 / (1.0 + std::exp(0.8 * (TMIN - 8.0)))) * 
                           (1.0 - 0.6 / (1.0 + std::exp(-0.8 * (TMAX - 34.0))));
    
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
    
    // DON = base_rate(T) × biomass × weather × temporal × RH
    double daily_DON_total_ug = base_rate * fungal_biomass_mg_m2 * weather_factor * 
                                temporal_factor * rh_factor;
    
    // Convert to concentration (µg/kg grain)
    double grain_yield_g_m2 = *SDWT;
    double daily_DON_ug_kg = (grain_yield_g_m2 > 0) ? 
                             (daily_DON_total_ug / grain_yield_g_m2) : 0.0;
    
    // Accumulate total DON over season
    static double cumulative_DON_ug_kg = 0.0;
    static int last_run_don = -1;
    static int last_yrsim_don = -1;
    
    // IMPORTANT: Save final values BEFORE resetting for new RUN
    static double saved_cumulative_DON = 0.0;
    static double saved_total_biomass = 0.0;
    static double saved_total_damaged_tissue = 0.0;
    static size_t saved_cohort_count = 0;
    
    if (last_run_don != -1 && last_run_don != RUN) {
        // Save the previous RUN's final values before reset
        saved_cumulative_DON = cumulative_DON_ug_kg;
        saved_total_biomass = total_biomass;
        saved_total_damaged_tissue = total_damaged_tissue;
        saved_cohort_count = cohorts.size();
        
        printf("\n>>> SAVED FINAL VALUES for RUN %d, YRSIM %d before reset:\n", last_run_don, last_yrsim_don);
        printf("    Final DON: %.2f µg/kg\n", saved_cumulative_DON);
        printf("    Final Biomass: %.4f g/m²\n", saved_total_biomass);
        printf("    Cohorts: %zu\n\n", saved_cohort_count);
    }
    
    // Reset cumulative DON when new run starts
    if (last_run_don != RUN) {
        cumulative_DON_ug_kg = 0.0;
        last_run_don = RUN;
        last_yrsim_don = YRSIM;
    }
    
    cumulative_DON_ug_kg += daily_DON_ug_kg;

    printf("YRDOY: %i, ZSTAGE: %.4f, dIdt: %.4f, HSDWT: %.4f, total_damaged_tissue: %.4f, cohorts.size(): %zu, WSDD: %.6f, total_biomass: %.6f, daily_DON: %.2f ug/kg, cumulative_DON: %.2f ug/kg\n",
           *YRDOY, ZSTAGE, dIdt, HSDWT, total_damaged_tissue,
           cohorts.size(), *WSDD, total_biomass, daily_DON_ug_kg, cumulative_DON_ug_kg);
    
    // ═══════════════════════════════════════════════════════════════════════
    // WRITE DAILY DON OUTPUT FILE (grouped by RUN and YRSIM like simulation_results)
    // ═══════════════════════════════════════════════════════════════════════
    static int last_run_file = -1;
    static int last_yrsim_file = -1;
    
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
    fdon << *YRDOY << ',' 
         << ZSTAGE << ',' 
         << dIdt << ',' 
         << total_biomass << ',' 
         << total_damaged_tissue << ',' 
         << cohorts.size() << ',' 
         << *WSDD << ',' 
         << daily_DON_ug_kg << ',' 
         << cumulative_DON_ug_kg << '\n';
    fdon.close();
    
    // ═══════════════════════════════════════════════════════════════════════
    // WRITE END-OF-SEASON DON SUMMARY (last value for each TRT)
    // ═══════════════════════════════════════════════════════════════════════
    
    std::string summary_fileName = "DON_concentration_summary.csv";
    
    // DEBUG: Print current RUN/YRSIM on first day of each treatment
    static int debug_last_run = -1;
    if (debug_last_run != RUN) {
        printf("\n>>> DEBUG: Treatment changed! RUN=%d, YRSIM=%d, YRDOY=%d\n", RUN, YRSIM, *YRDOY);
        printf("    Previous: last_run_summary=%d, last_yrsim_summary=%d\n", last_run_summary, last_yrsim_summary);
        printf("    cumulative_DON so far: %.2f µg/kg\n\n", cumulative_DON_ug_kg);
        debug_last_run = RUN;
    }
    
    // When treatment changes, write summary for previous treatment FIRST before updating the map
    if (last_run_summary != -1 && (RUN != last_run_summary || YRSIM != last_yrsim_summary)) {
        // FIRST: Update the map with the SAVED final values from before the reset
        std::pair<int,int> prev_key = std::make_pair(last_run_summary, last_yrsim_summary);
        if (treatment_summaries.find(prev_key) != treatment_summaries.end()) {
            treatment_summaries[prev_key].cumulative_DON = saved_cumulative_DON;
            treatment_summaries[prev_key].total_biomass = saved_total_biomass;
            treatment_summaries[prev_key].total_damaged_tissue = saved_total_damaged_tissue;
            treatment_summaries[prev_key].cohort_count = saved_cohort_count;
            
            TreatmentSummary& prev_summary = treatment_summaries[prev_key];
            
            printf("\n>>> UPDATED MAP with saved values for RUN %d:\n", last_run_summary);
            printf("    Saved DON: %.2f µg/kg\n", saved_cumulative_DON);
            printf("    Saved Biomass: %.4f g/m²\n\n", saved_total_biomass);
            
            // NOW write immediately to file (only if not already written)
            if (!summary_written[prev_key]) {
                std::string summary_fileName = "DON_concentration_summary.csv";
                std::ofstream fsum;
                if (!header_written_summary) {
                fsum.open(summary_fileName);
                fsum << "RUN,YRSIM,Last_YRDOY,Final_DON_Concentration_ug_kg\n";
                header_written_summary = true;
            } else {
                fsum.open(summary_fileName, std::ios::app);
            }
            
            fsum << last_run_summary << ','
                 << last_yrsim_summary << ','
                 << prev_summary.last_YRDOY << ','
                 << prev_summary.cumulative_DON << '\n';
            fsum.close();
            
            summary_written[prev_key] = true;
            
            printf("\n>>> TREATMENT SUMMARY WRITTEN (on change) - RUN %d, YRSIM %d <<<\n", 
                   last_run_summary, last_yrsim_summary);
            printf("    Final DON Concentration: %.2f µg/kg\n", prev_summary.cumulative_DON);
            printf("    Final Fungal Biomass: %.4f g/m²\n", prev_summary.total_biomass);
            printf("    Final YRDOY: %d\n", prev_summary.last_YRDOY);
            printf("    Summary written to DON_concentration_summary.csv\n\n");
            }  // Close the if (!summary_written[prev_key]) block
        }
    }
    
    // NOW update the map with current treatment's latest data EVERY day
    // This ensures we always have the most recent values, including for the last treatment
    std::pair<int,int> trt_key = std::make_pair(RUN, YRSIM);
    
    TreatmentSummary current_summary;
    current_summary.last_YRDOY = *YRDOY;
    current_summary.cumulative_DON = cumulative_DON_ug_kg;
    current_summary.total_biomass = total_biomass;
    current_summary.total_damaged_tissue = total_damaged_tissue;
    current_summary.cohort_count = cohorts.size();
    treatment_summaries[trt_key] = current_summary;
    
    // Write current treatment summary ONLY once at the end of the season
    // We detect end-of-season when ZSTAGE reaches near physiological maturity (>= 89)
    // This captures the final accumulated DON for all treatments including the last one
    if (ZSTAGE >= 89.0 && !summary_written[trt_key]) {
        std::string summary_fileName = "DON_concentration_summary.csv";
        std::ofstream fsum;
        if (!header_written_summary) {
            fsum.open(summary_fileName);
            fsum << "RUN,YRSIM,Last_YRDOY,Final_DON_Concentration_ug_kg\n";
            header_written_summary = true;
        } else {
            fsum.open(summary_fileName, std::ios::app);
        }
        
        fsum << RUN << ','
             << YRSIM << ','
             << current_summary.last_YRDOY << ','
             << current_summary.cumulative_DON << '\n';
        fsum.close();
        
        summary_written[trt_key] = true;
        
        printf("\n>>> TREATMENT SUMMARY WRITTEN (at maturity) - RUN %d, YRSIM %d <<<\n", RUN, YRSIM);
        printf("    Final DON Concentration: %.2f µg/kg\n", current_summary.cumulative_DON);
        printf("    Final Fungal Biomass: %.4f g/m²\n", current_summary.total_biomass);
        printf("    Final YRDOY: %d, ZSTAGE: %.2f\n", current_summary.last_YRDOY, ZSTAGE);
        printf("    Summary written to DON_concentration_summary.csv\n\n");
    }
    
    // Update tracking variables
    last_run_summary = RUN;
    last_yrsim_summary = YRSIM;
    last_YRDOY = *YRDOY;
    
    return 1;
}

int couplingOutput(int *doy) {
    
}