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
double time_step_size;
double S , W , I ;
using Vec3 = std::array<double,3>;
Vec3 y = {S, W, I};
double dSdt, dWdt, dIdt;
double daily_dI;           // daily increment in infection (for cohort biomass calculation)
double daily_dW;           // daily increment in wheat (for cohort biomass calculation)


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
    static int SUSTAGE = 0, first_day_sus = 0;

    
    double Rain_t = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "RAIN");
    int dap = fio->getReal("PEST", "YRPLT");
    
    
    int is_rainy_day = (Rain_t > rain_threshold) ? 1 : 0;
    
    if(ZSTAGE >= 51 && SUSTAGE == 0){
        first_day_sus = *YRDOY;
        SUSTAGE = 1;
    }
    
    day = *YRDOY - first_day_sus;
    

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
    
    double SRAD_t = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "SRAD");
    double TMAX_t = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "TMAX");
    double TMIN_t = fio->getRealYrdoy("WTH", std::to_string(*YRDOY), "TMIN");
    
    static int YRSIMp = -1;

    std::string YRSIM_str = std::to_string(YRSIM);
    std::string year = YRSIM_str.substr(0, YRSIM_str.size() - 3);
    std::string fileName = "simulation_results_" + year + ".csv";

    // Check if new simulation year started
    if (YRSIM != YRSIMp) {
        // Remove any previous results for this simulation year
        std::remove(fileName.c_str());
        
        // Initialize new results file with header
        std::ofstream fout(fileName);
        fout << "YRDOY,Spores_in_Air,Spores_on_Wheat,Infection,Phen_Stage,Soil_Water(%),SRAD,TMAX,TMIN,RAIN\n";
        fout.close();

        // Reset disease susceptibility stage tracker
        SUSTAGE = 0;
        YRSIMp = YRSIM;
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
    static int simulation_year = -1;
    
    // Reset cohorts when new simulation year starts
    FlexibleIO *fio = FlexibleIO::getInstance();
    int YRSIM = fio->getReal("PEST", "YRSIM");
    double ZSTAGE = fio->getReal("PEST", "ZSTAGE");
    if (simulation_year != YRSIM) {
        cohorts.clear();
        total_damaged_tissue = 0.0;
        simulation_year = YRSIM;
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

    printf("YRDOY: %i, ZSTAGE: %.4f, dIdt: %.4f, HSDWT: %.4f, total_damaged_tissue: %.4f, cohorts.size(): %zu, WSDD: %.6f, total_biomass: %.6f\n",
           *YRDOY, ZSTAGE, dIdt, HSDWT, total_damaged_tissue,
           cohorts.size(), *WSDD, total_biomass);
    
    last_YRDOY = *YRDOY;
    
    return 1;
}
int couplingOutput(int *doy) {
    
}