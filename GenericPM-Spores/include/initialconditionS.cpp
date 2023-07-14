#include "initialconditionS.h"

#include<sstream>
#include<iostream>
#include<fstream>

int InitialConditionS::qtdS = 0;

/** Calculate the daily favorability based on temp * wetness favorability*/
void InitialConditionS::rate() {
    if (!stop) {
        dailyFavorability = Utilities::temperatureFavorability(Basic::getWeather()->getTMean(),
                                                               cloudf.getDisease()->getTemperatureFavorabilitySet()) 
                            *
                            Utilities::wetnessFavorability(Basic::getWeather()->getWetDur(), cloudf.getDisease()->getWetnessFunction()); //
    }
    cloudf.rate();
}

/** Accumulate the daily favorability. If the this value hits the pre-determinated value, stop the process */
void InitialConditionS::integration() {
    if (!stop) {
        acumulateFavorability += dailyFavorability;
        if (acumulateFavorability >= cloudf.getDisease()->getAcumulateFavorability()) {
            cloudf.setFirstSporeCloud(cloudf.getDisease()->getInitialInoculum());
            stop = true;
        }

        std::ostringstream convert;
        convert << Basic::getWeather()->getYearDoy() << "," << acumulateFavorability;
        Basic::output.push_back(convert.str());

    }
    cloudf.integration();

}

void InitialConditionS::output() {
    std::ostringstream convert;
    convert << "Cpp_InitialCondition_" << getID() << ".txt";
    Basic::getOutput(convert.str());

    // Speedup the model removing outputs
    //std::cout << "\nInitialCondition " << getID() << ":";
    //for(unsigned int i=0; i<Basic::output.size(); i++)
    //{
    //    std::cout << Basic::output[i] << std::endl;
    //}
    cloudf.output();
}
