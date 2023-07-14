#include "lesioncohortS.h"
#include "cropinterfaceS.h"
#include<string>
#include<sstream>
#include<iostream>

int LesionCohortS::qtdS = 0;

void LesionCohortS::integration() {

    Disease *disease = cloudo->getDisease();

    if (getOrganHealthAreaProportion() > 0.01) {      
        if(dailyVisibleAreaGrow>0) {
            visibleArea = dailyVisibleAreaGrow * lesionsInThisCohort;
        }
        if(dailyInvisibleAreaGrow>0) {
            invisibleArea = (dailyInvisibleAreaGrow-dailyVisibleAreaGrow) * lesionsInThisCohort;
        }
        totalArea = visibleArea + invisibleArea;

        if (isLatentPeriodS()) {
            latentArea = totalArea;
            infectionArea = necroticArea = 0;
        } else if (isInfectionPeriodS()) {
            infectionArea = totalArea;
            latentArea = necroticArea = 0;
        } else {
            necroticArea = totalArea;
            infectionArea = latentArea = 0;
        }

        cloudo->addSporesCreatedS(newSpores);

        physiologicalDaysAcumm += physiologicalDay;

        std::ostringstream convert;
        convert << Basic::getWeather()->getYearDoy() << "," << totalArea << "," << lesionsInThisCohort << "," << getPhysiologicalDaysAcumm() << ","
                << getOrganDiseasedAreaProportion() << "," << latentArea << "," << infectionArea << "," << necroticArea << ","
                << newSpores << "," << util.temperatureFavorability(Basic::getWeather()->getTMean(),
                disease->getTemperatureFavorabilitySet()) << "," << dailyVisibleAreaGrow << "," << dailyInvisibleAreaGrow;
        Basic::output.push_back(convert.str());

        newSpores=0;
    }

}

int LesionCohortS::getVisibleLesions() {
    if (isInfectionPeriodS() || isNecroticPeriodS())
        return this->lesionsInThisCohort;
    else
        return 0;
}

void LesionCohortS::output() {
    std::ostringstream convert;
    convert << "Cpp_LesionCohort_" << getID() << ".txt";
    Basic::getOutput(convert.str());
}

void LesionCohortS::rate() {
    //double dailyVisibleGrowRate = 0, dailyInvisibleGrowRate = 0;
    Disease *disease = cloudo->getDisease();
    physiologicalDay = util.temperatureFavorability(
                            Basic::getWeather()->getTMean(),
                            disease->getTemperatureFavorabilitySet());
    // Thinking on: cumsum(runif(25, min = 0.01, max = 0.1))
    dailyInvisibleAreaGrow  = util.growthFunction(disease->getInvisibleGrowthFunction(), 
                                                  getPhysiologicalDaysAcumm()) 
/*                              * disease->getHostFactor() 
                              * totalArea 
                              * getOrganHealthAreaProportion()*/;
    dailyVisibleAreaGrow    = util.growthFunction(disease->getVisibleGrowthFunction(), 
                                                  getPhysiologicalDaysAcumm()) 
/*                              * disease->getHostFactor()
                              * totalArea
                              * getOrganHealthAreaProportion()*/;
    if (isLatentPeriodS()) { 
        dailyVisibleAreaGrow = 0;
    } else if (isNecroticPeriodS()) { 
        dailyVisibleAreaGrow = dailyInvisibleAreaGrow = 0;
    }
    
    newSpores = 0;
    if (getOrganHealthAreaProportion() > 0.01 && isInfectionPeriodS() &&
            Basic::getWeather()->getWetDur() >= disease->getWetnessThreshold()) // && Basic::getWeather()->getTMean() > 20
    {
        newSpores = (lesionsInThisCohort * disease->getDailySporeProductionPerLesion() * 
                     util.trapezoidalFunction(getAge(), disease->getCohortAgeSet()) *
                     disease->getSporulationCrowdingFactor(getOrganDiseasedAreaProportion()));
    }
    if(organHealthAreaProportion < 0.7) {
        dailyVisibleAreaGrow = dailyInvisibleAreaGrow = 0; // newSpores = 0;
    }
}

bool LesionCohortS::isLatentPeriodS() {
    Disease *disease = cloudo->getDisease();
    if (getPhysiologicalDaysAcumm() <= disease->getLatentPeriod()) {
        return true;
    }
    return false;
}

bool LesionCohortS::isInfectionPeriodS() {
    Disease *disease = cloudo->getDisease();
    if (getPhysiologicalDaysAcumm() > disease->getLatentPeriod() && getPhysiologicalDaysAcumm() <= (disease->getLatentPeriod() + disease->getInfectionPeriod())) {
        return true;
    }
    return false;
}

bool LesionCohortS::isNecroticPeriodS() {
    Disease *disease = cloudo->getDisease();
    if (getPhysiologicalDaysAcumm() > (disease->getLatentPeriod() + disease->getInfectionPeriod())) {
        return true;
    }
    return false;
}