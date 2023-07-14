#include <iostream>
#include "../GenericPM/include/simulator.h"
#include "../GenericPM/include/utilities.h"
#include <cmath>
#include <vector>
#include <fstream>
#include "../FlexibleIO/Data/FlexibleIO.hpp"
#include "cinterfaceS.h"

using namespace std;

double AREALF = Utilities::runExpressionFunction(1,"50000");

extern "C"
{
    int couplingInitSpore(int *YRDOY, int *YRPLT);
    int couplingRateSpore(int *YRDOY, float *SL1);
    double couplingIntegrationSpore(int *YRDOY, int *YRPLT);
    int couplingOutputSpore(int *doy);
}

int couplingInitSpore(int *YRDOY, int *YRPLT)
{
    // Set the start day for Disease Model
    // Get an instance of Simulator
    SimulatorSpore *sS = SimulatorSpore::newInstance();
    // Set the start day for Disease Model
    sS->setCurrentYearDoy(*YRDOY);
    // Set the sowing/planting date
    sS->getCropInterface()->setPlantingDate(*YRPLT);
    sS->getCropInterface()->setOrganArea(1, AREALF);

    return (1);
}

int couplingRateSpore(int *YRDOY, float *SL1)
{
    // Temporary variable used for computations 
    float temp = 0, newOrgan = 0;
    double CloudField = 0;
    // Get an instance of Simulator
    SimulatorSpore *sS = SimulatorSpore::getInstance();
    newOrgan = sS->getCropInterface()->getOrgansQtd()+1;

    sS->getCropInterface()->setOrganArea(1, AREALF);
    //sS->getCropInterface()->setSoilMoisture(*SL1);
    // Set the current YearDOY for next Disease step computation
    sS->updateCurrentYearDoy(*YRDOY);

    sS->getCropInterface()->setOrganArea(newOrgan, AREALF);

    
    // Feed the Disease Model with weather information
    Weather::getInstance()->update();
    // Disease Simulator Rate
    sS->rate();
    
    return (1);
}

double couplingIntegrationSpore(int *YRDOY,  int *YRPLT)
{
    // Temporary variable used for computations
    float dArea = 0, tArea = 0, pDArea = 0, sArea = 0, pclaCalc = 0;
    double CloudField = 0;

    // Get an instance of SimulatorSpore
    SimulatorSpore *sS = SimulatorSpore::getInstance();

    // Call the Disease Model Integration function
    sS->integration();

    if (sS->getPlants().size() > 0 && sS->getPlants()[0].getOrgans().size() > 0)
    {
        dArea = sS->getPlants()[0].getDiseaseArea();
        tArea = sS->getPlants()[0].getTotalArea();
        sArea = sS->getPlants()[0].getSenescenceArea();
        pDArea = (dArea / (tArea - sArea) * 100);
        if (pDArea > 99.9)
        {
            pDArea = 99.9;
        }
        for (int i = 0; i < sS->getPlants()[0].getOrgans().size(); i++)
        {
            if (sS->getPlants()[0].getOrgans().at(i).getVisibleLesions() /
                    (sS->getPlants()[0].getOrgans().at(i).getTotalArea() - sS->getPlants()[0].getOrgans().at(i).getSenescenceArea()) >=
                10)
            {
                pclaCalc += fmax(0, sS->getPlants()[0].getOrgans().at(i).getTotalArea() - sS->getPlants()[0].getOrgans().at(i).getSenescenceArea());
            }
        }
        CloudField = sS->getPlants()[0].getCloudsP()[0].getCloudF()->getValue();

    }
    std::ofstream out;    
    out.open("Daily_CloudF_"+ std::to_string(*YRPLT) +".txt", std::ofstream::out | std::ofstream::app);
    std::cout << *YRDOY<< " " << CloudField << std::endl;

    out << *YRDOY<< " " << CloudField << std::endl;
    out.close();
    
    return (CloudField);
}

int couplingOutputSpore(int *doy)
{
    // Get an instance of SimulatorSpore
    SimulatorSpore *sS = SimulatorSpore::getInstance();
    // Request disease outputs to be written in files
    sS->output();

    return (1);
}

