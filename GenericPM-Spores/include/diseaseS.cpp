#include "../../GenericPM/include/disease.h"
#include<cmath>
#include<iostream>
#include "../../FlexibleIO/Data/FlexibleIO.hpp"
#include "../../GenericPM/include/utilities.h"

std::vector<Disease*> Disease::listDiseasesS;

double Disease::getSporulationCrowdingFactor(double proportionDiseaseArea) {
    double a = (1 / (sporulationCrowdingFactorsSet[0] + sporulationCrowdingFactorsSet[1] * pow(proportionDiseaseArea, sporulationCrowdingFactorsSet[2])));
    return (fmin(a,1));
}

int Disease::newLesions(double cloudDensity, double healthyAreaProportion) {
    Utilities util;
    double newLesions = 0;
    double fitWetnessThreshold = getWetnessThreshold();

    if (healthyAreaProportion > 0 && Basic::getWeather()->getWetDur() >= fitWetnessThreshold) {
        newLesions = (cloudDensity * healthyAreaProportion * getInfectionEfficiency() *
                util.temperatureFavorability(Basic::getWeather()->getTMean(),
                                             getTemperatureFavorabilitySet()) *
                util.wetnessFavorability(Basic::getWeather()->getWetDur(),getWetnessFunction())) >0 ? (cloudDensity * healthyAreaProportion * getInfectionEfficiency() *
                util.temperatureFavorability(Basic::getWeather()->getTMean(),
                                             getTemperatureFavorabilitySet()) *
                util.wetnessFavorability(Basic::getWeather()->getWetDur(),getWetnessFunction())) : 0;
            /*std::cout << " 1: " << newLesions << " 2: " << cloudDensity << " 3: " << healthyAreaProportion << " 4: " << getInfectionEfficiency() <<
                " 5: " << util.temperatureFavorability(Basic::getWeather()->getTMean(),getTemperatureFavorabilitySet()) << " 6: " <<
                util.wetnessFavorability(Basic::getWeather()->getWetDur()) << " 7: " << Basic::getWeather()->getWetDur()<< std::endl; */
    //newLesions = newLesions * Utilities::runExpressionFunction(Basic::getWeather()->getRh(),getRhFactor());
    //newLesions= newLesions *  Utilities::runExpressionFunction(Basic::getWeather()->getRh(),getRhFactor());
    //std::cout<<newLesions<< " exp : "<<newLesions *  Utilities::runExpressionFunction(Basic::getWeather()->getRh(),getRhFactor()) <<std::endl; 
   // std::cout<<"rhfacetor "<<getRhFactor()<<" RH : "<<Basic::getWeather()->getRh()<<std::endl;
    }
    return newLesions;
}
