#include "initialconditionS.h"

#include<sstream>
#include<iostream>
#include<fstream>

int InitialConditionS::qtdS = 0;

/** Calculate the daily favorability based on temp * wetness favorability*/
void InitialConditionS::rate() {
    if (!stop) {
        dailyFavorability = Utilities::temperatureFavorability(BasicS::getWeather()->getTMean(),
                                                               cloudf.getDisease()->getTemperatureFavorabilitySet()) 
                            *
                            Utilities::wetnessFavorability(BasicS::getWeather()->getWetDur(), cloudf.getDisease()->getWetnessFunction()); //
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
        convert << BasicS::getWeather()->getYearDoy() << "," << acumulateFavorability;
        BasicS::output.push_back(convert.str());

    }
    cloudf.integration();

}

void InitialConditionS::output() {
    std::ostringstream convert;
    convert << "Cpp_InitialCondition_" << getID() << ".txt";
    BasicS::getOutput(convert.str());

    // Speedup the model removing outputs
    //std::cout << "\nInitialCondition " << getID() << ":";
    //for(unsigned int i=0; i<BasicS::output.size(); i++)
    //{
    //    std::cout << BasicS::output[i] << std::endl;
    //}
    cloudf.output();
}
