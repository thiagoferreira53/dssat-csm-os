#include "../../GenericPM/include/weather.h"
#include "simulatorS.h"
#include "../../FlexibleIO/Data/FlexibleIO.hpp"
#include<new>
#include<iostream>

Weather::Weather() {
    //update();
}

Weather* Weather::instance = nullptr;

Weather* Weather::getInstance() {
    if (instance == nullptr)
        instance = new Weather();
    return instance;
}

bool Weather::update() {
    try {
        yearDoy = SimulatorSpore::getInstance()->getCurrentYearDoy();
        year = yearDoy / 1000;
        doy = yearDoy - (year * 1000);
        

        sRad = FlexibleIO::getInstance()->getRealYrdoy("WTH", std::to_string(yearDoy), "SRAD");
        tMax = FlexibleIO::getInstance()->getRealYrdoy("WTH", std::to_string(yearDoy), "TMAX");
        tMin = FlexibleIO::getInstance()->getRealYrdoy("WTH", std::to_string(yearDoy), "TMIN");
        rain = FlexibleIO::getInstance()->getRealYrdoy("WTH", std::to_string(yearDoy), "RAIN");
        hRH90 = FlexibleIO::getInstance()->getRealYrdoy("WTH", std::to_string(yearDoy), "RH90");
        //rh = FlexibleIO::getInstance()->getRealYrdoy("WTH", std::to_string(yearDoy), "RH");
        tMean = (tMax + tMin) / 2;
        par = 0;
        wetDur = 0;
        hWetDur = 0;
        rhMax = 0;
        rhMin = 0;
        //        getRealYrdoy(&yearDoy,(char *) "RHUM=",&rhMean);  //?????
        rhMean = 0;
                //printf(" Weather::update(). yearDoy: %i - sRad: %f - tMax: %f - tMin: %f - rain: %f - rhMean: %f - hRH90: %i - tMean: %f\n", 
                //        yearDoy, sRad, tMax, tMin, rain, rhMean, hRH90, tMean);
    } catch (bool) {
        return false;
    }

    return true;
}
