#include "diseaseS.h"
#include<cmath>
#include<iostream>
#include "../../FlexibleIO/Data/FlexibleIO.hpp"
#include "../../GenericPM/include/utilities.h"

std::vector<DiseaseS*> DiseaseS::listDiseasesS;

double DiseaseS::getSporulationCrowdingFactorS(double proportionDiseaseArea) {
    double a = (1 / (sporulationCrowdingFactorsSet[0] + sporulationCrowdingFactorsSet[1] * pow(proportionDiseaseArea, sporulationCrowdingFactorsSet[2])));
    return (fmin(a,1));
}

int DiseaseS::newLesionsS(double cloudDensity, double healthyAreaProportion) {
    Utilities util;
    double newLesionsS = 0;
    double fitWetnessThreshold = getWetnessThreshold();

    if (healthyAreaProportion > 0 && BasicS::getWeather()->getWetDur() >= fitWetnessThreshold) {
        newLesionsS = (cloudDensity * healthyAreaProportion * getInfectionEfficiency() *
                util.temperatureFavorability(BasicS::getWeather()->getTMean(),
                                             getTemperatureFavorabilitySet()) *
                util.wetnessFavorability(BasicS::getWeather()->getWetDur(),getWetnessFunction())) >0 ? (cloudDensity * healthyAreaProportion * getInfectionEfficiency() *
                util.temperatureFavorability(BasicS::getWeather()->getTMean(),
                                             getTemperatureFavorabilitySet()) *
                util.wetnessFavorability(BasicS::getWeather()->getWetDur(),getWetnessFunction())) : 0;
            /*std::cout << " 1: " << newLesionsS << " 2: " << cloudDensity << " 3: " << healthyAreaProportion << " 4: " << getInfectionEfficiency() <<
                " 5: " << util.temperatureFavorability(BasicS::getWeather()->getTMean(),getTemperatureFavorabilitySet()) << " 6: " <<
                util.wetnessFavorability(BasicS::getWeather()->getWetDur()) << " 7: " << BasicS::getWeather()->getWetDur()<< std::endl; */
    //newLesionsS = newLesionsS * Utilities::runExpressionFunction(BasicS::getWeather()->getRh(),getRhFactor());
    //newLesionsS= newLesionsS *  Utilities::runExpressionFunction(BasicS::getWeather()->getRh(),getRhFactor());
    //std::cout<<newLesionsS<< " exp : "<<newLesionsS *  Utilities::runExpressionFunction(BasicS::getWeather()->getRh(),getRhFactor()) <<std::endl; 
   // std::cout<<"rhfacetor "<<getRhFactor()<<" RH : "<<BasicS::getWeather()->getRh()<<std::endl;
    }
    return newLesionsS;
}
