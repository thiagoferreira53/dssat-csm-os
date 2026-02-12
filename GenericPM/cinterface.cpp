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
#include "include/utilities.h"
#include "include/toxins.h"
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

double kr;
double M_thresh;
double kd;
double kg;
double alpha;
double kdw;
double k_inf;
double k_rec;
double rain_threshold;
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
        S = 100, W = 0.0, I = 0.0;
        y = {S, W, I};
        
    // Initialize DON model global variables
    heading_yrdoy_global = -1;
    heading_detected_global = false;
    first_day_sus_global = 0;
    
    return (1);
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
    
    // ═══════════════════════════════════════════════════════════════════════
    // DAILY INTEGRATION
    // ═══════════════════════════════════════════════════════════════════════
    
    // Anther maturation phenology (fraction of anthers mature)
    double A = 0.022900;
    double B = 3.612468;
    double C = 0.464022;
    double anther_prop_t = A * std::pow(day, B) * std::exp(-C * day);
    anther_prop_t = std::max(0.0, std::min(1.0, anther_prop_t));
    
    // Store old state values before integration
    double S_old = y[0];
    double W_old = y[1];
    double I_old = y[2];
    
    // Calculate rates of change for current day
    double spore_release_rate = kr * std::max(0.0, SW - M_thresh);
    double temp_fact = Utilities::temperatureFactor(TAVG);
    
    double rate_S = std::max(0.0, spore_release_rate - (kd * S_old) - (kg * S_old));
    double rate_W = std::max(0.0, (alpha * kg * S_old) - (kdw * W_old));
    double rate_I = std::max(0.0, (k_inf * W_old * is_rainy_day * anther_prop_t * temp_fact) - (k_rec * I_old));
    
    // Update state: state_new = state_old + rate * 1.0 day
    y[0] = std::max(0.0, S_old + rate_S);  // S (spores in air) - cloudF
    y[1] = std::max(0.0, W_old + rate_W);  // W (spores on wheat) - cloudP
    y[2] = std::max(0.0, I_old + rate_I);  // I (infection level) - cloudO
    
    // Store rates and daily changes for output and biomass calculation
    // (mimicking original code naming for compatibility)
    dSdt = y[0];  // new S
    dWdt = y[1];  // new W
    dIdt = y[2];  // new I
    
    daily_dI = y[2] - I_old;  // Daily change in infection
    daily_dW = y[1] - W_old;  // Daily change in wheat spores

    
    int YRSIM = fio->getReal("PEST", "YRSIM");
    int RUN = fio->getInteger("PEST", "RUN");
    
    double SRAD_t = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "SRAD");
    double TMAX_t = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "TMAX");
    double TMIN_t = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "TMIN");
    
    std::string YRSIM_str = std::to_string(YRSIM);
    std::string year = YRSIM_str.substr(0, YRSIM_str.size() - 3);
    std::string fileName = "simulation_results_RUN" + std::to_string(RUN) + "_" + year + ".csv";
    
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
    
    Toxins::FHB_DON don_result = Toxins::calculateDailyDON(
        *YRDOY, RUN, YRSIM, ZSTAGE,
        total_biomass, total_damaged_tissue, cohorts.size(),
        *WSDD, *SDWT, dIdt,
        heading_yrdoy_global, heading_detected_global
    );
    
    return 1;
}

int couplingOutput(int *doy) {
    return 1;
}