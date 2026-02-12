/**
 * @file toxins.h
 * 
 * @author Thiago Berton Ferreira
 * @author Willingthon Pavan (wpavan.us@gmail.com)
 * @author Jose Mauricio Cunha Fernandes (jmauricio.fernandes@icloud.com)
 * 
 * @copyright Copyright (c) 2017–2025, DSSAT Foundation
 * @license BSD-3-Clause. See the LICENSE file in the root folder for details.
 */
#ifndef TOXINS_H
#define TOXINS_H

#include <string>

class Toxins {
public:
    // Structure to hold DON calculation results
    struct FHB_DON {
        double daily_DON_total_ug;
        double cumulative_DON_total_ug;
        double daily_DON_ug_kg;
        double cumulative_DON_ug_kg;
    };
    
    // Calculate daily DON production
    static FHB_DON calculateDailyDON(
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
    );

};

#endif // TOXINS_H
