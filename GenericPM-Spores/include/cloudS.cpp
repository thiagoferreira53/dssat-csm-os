#include "cloudS.h"
#include "diseaseS.h"
#include <cmath>
#include <iostream>
#include "simulatorS.h"
//#include "../Coupling-Interface/CPP/coupling_IO_cpp.h"
//#include "../../FlexibleIO/Data/FlexibleIO.hpp"

void CloudS::integration()
{
    int qtdS = 0;
    if (sporesToBeRemoved > 0)
    { // Remove from cloud spores used to infect tissue
        removeSporesCloudS(sporesToBeRemoved);
        /*        for (unsigned int i = 0; i < values.size() && sporesToBeRemoved>0; i++) {
            qtdS = fmin(values[i],sporesToBeRemoved);
            sporesToBeRemoved -= qtdS;
            values[i] = values[i] - qtdS;
        }*/
    }
    //std::cout<<"Rain: "<< BasicS::getWeather()->getRain()<<std::endl;
    //std::cout<<"New function: "<< Utilities::dispersalRainFunction(BasicS::getWeather()->getRain(),disease->getDispersionFreequency())<<std::endl;
    values.push_back(Utilities::dispersalRainFunction(BasicS::getWeather()->getRain(), disease->getDispersionFreequency()) * sporesCreated);
    //std::cout << disease->getDispersionFreequency() * sporesCreated << std::endl;
    sporesCreated = sporesToBeRemoved = 0;
}

double CloudS::getValue()
{
    double sum = 0;
    for (unsigned int i = 0; i < values.size(); i++)
        sum += values[i];

    return sum;
}

void CloudS::removeSporesCloudS(double toBeRemoved)
{
    double total = getValue();
    for (unsigned int i = 0; i < values.size() && total > 0; i++)
    {
        values[i] -= (toBeRemoved * (values[i] / total));
    }
}

void CloudS::removeSporesCloudByRainS(double percent)
{
    double oldValue = 0;
    for (unsigned int i = 0; i < values.size(); i++)
    {
        oldValue = values[i];
        values[i] = (oldValue * percent);
    }
}

void CloudS::removeSporesCloudFByAgeS(void)
{
    int yearDoy = SimulatorSpore::getInstance()->getCurrentYearDoy();

    for (int i = values.size(); i > 0; i--)
    {
        values[values.size() - i] = values[values.size() - i] * (1.61 * exp(-0.369 * (values.size() - i + 1))) > 0 ? values[values.size() - i] * (1.61 * exp(-0.369 * (values.size() - i + 1))) : 0;
    }
}
void CloudS::removeSporesCloudPByAgeS(void)
{

    for (int i = values.size(); i > 0; i--)
    {
        values[values.size() - i] = values[values.size() - i] * (1.40 * exp(-0.2030 * (values.size() - i + 1))) > 0 ? values[values.size() - i] * (1.40 * exp(-0.2030 * (values.size() - i + 1))) : 0;
    }
}
void CloudS::removeSporesCloudOByAgeS(void)
{

    for (int i = values.size(); i > 0; i--)
    {
        values[values.size() - i] = values[values.size() - i] * (1.4268 * exp(-0.2184 * (values.size() - i + 1))) > 0 ? values[values.size() - i] * (1.4268 * exp(-0.2184 * (values.size() - i + 1))) : 0;
    }
}

/*void CloudS::removeSporesCloudFByAgeUvS(void)
{
    int yearDoy = SimulatorSpore::getInstance()->getCurrentYearDoy();
    float tMax = 0;
    float tMin = 0;
    tMax = FlexibleIO::getInstance()->getRealYrdoy("WTH", std::to_string(yearDoy), "TMAX");
    tMin = FlexibleIO::getInstance()->getRealYrdoy("WTH", std::to_string(yearDoy), "TMIN");
   
    //std::cout <<yearDoy<< " uv index: "<<tMax<<" "<<tMin <<" "<< (10 / (tMax - tMin))<< std::endl;
    for (int i = values.size(); i > 0; i--)
    {
        //std::cout << "removed F by age i:" << i << " - v-i: " << values.size() - i << " - v-i+1: " << values.size() - i + 1 << " - " << values[values.size() - i];
        values[values.size() - i] = (values[values.size() - i] * (1.61 * exp(-0.369 * (values.size() - i + 1))))*(10 / (tMax - tMin)) > 0 ? (values[values.size() - i] * (1.61 * exp(-0.369 * (values.size() - i + 1))))*(10 / (tMax - tMin)) : 0;
        //std::cout << " - " << values[values.size() - i] << std::endl;
    }
}*/