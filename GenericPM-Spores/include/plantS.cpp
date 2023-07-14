#include "plantS.h"
#include "simulatorS.h"

#include<sstream>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>

int PlantS::qtdS = 0;
int PlantS::firstOutputCallS = 0;

PlantS::PlantS() {
    std::vector<InitialConditionS> &vectIC = SimulatorSpore::getInstance()->getInitialConditions();
    InitialConditionS *ic;

    for (unsigned int i = 0; i < vectIC.size(); i++) {
        ic = &vectIC[i];
        cloudsP.emplace_back(ic->getCloud()->getDisease(), ic->getCloud());
    }
    //BasicS::output.push_back("Day, PlantArea, DiseaseArea, Density, Severity, LatentDArea, InfectionDArea, NecroticDArea,SenescenceArea");
}

void PlantS::integration() {
    totalArea = diseaseArea = latentDiseaseArea = infectionDiseaseArea = necroticDiseaseArea = visibleDiseaseArea = invisibleDiseaseArea = totalLesions = visibleLesions = senescenceArea = 0;
    int newOrgan = 0;
    double cloudOValue = 0, cloudPValue = 0, cloudFvalue = 0;
    OrganS *o;
    for (unsigned int i = 0; i < organs.size(); i++) {
        o = &organs[i];
        if(o->getSenescenceArea() < o->getTotalArea()) {
            o->integration();
            diseaseArea += o->getDiseaseArea();
            latentDiseaseArea += o->getLatentDiseaseArea();
            infectionDiseaseArea += o->getInfectionDiseaseArea();
            necroticDiseaseArea += o->getNecroticDiseaseArea();
            visibleDiseaseArea += o->getVisibleDiseaseArea();
            invisibleDiseaseArea += o->getInvisibleDiseaseArea();
            visibleLesions += o->getVisibleLesions();
            totalLesions += o->getTotalLesions();
        } else {
            //o->cloudIntegrationS();
            //printf("## OrganS: %d died - integration!!!\n",o->getOrganNumber());
        }
        totalArea += o->getTotalArea();
        senescenceArea += o->getSenescenceArea();
        cloudOValue += o->cloudAmountS();
    }

    CloudPS *cloud;
    for (unsigned int i = 0; i < cloudsP.size(); i++) {
        cloud = &cloudsP[i];
        cloud->integration();
    }

    newOrgan = SimulatorSpore::getInstance()->getCropInterface()->hasNewOrgan();
    if (newOrgan > 0) {
        //printf("Creating new organ: %i\n",newOrgan);
        organs.emplace_back(cloudsP, newOrgan, SimulatorSpore::getInstance()->getCropInterface()->getOrganArea(newOrgan));
    }

    cloudPValue = cloud->getValue();
    cloudFvalue = cloud->getCloudF()->getValue();

    std::ostringstream convert;
    //Plant, YearDoy, TotalArea, Senesced, Diseased, VisibleArea, InvisibleArea, TotalLesions, CloudOS, CloudPS, CloudFS
    convert << ID << "," << BasicS::getWeather()->getYearDoy() << "," << totalArea << "," << senescenceArea << "," << diseaseArea << "," 
            << visibleDiseaseArea << "," << invisibleDiseaseArea << "," << totalLesions << "," 
            << latentDiseaseArea << "," << infectionDiseaseArea << "," << necroticDiseaseArea << ","
            << Utilities::formatDouble(cloudOValue) << "," << Utilities::formatDouble(cloudPValue) << "," 
            << Utilities::formatDouble(cloudFvalue);
    BasicS::output.push_back(convert.str());
    //std::cout << ID << "," << BasicS::getWeather()->getYearDoy() << "," << totalArea << "," << senescenceArea << "," << diseaseArea << "," 
    //        << visibleDiseaseArea << "," << invisibleDiseaseArea << "," << totalLesions << "," 
    //        << latentDiseaseArea << "," << infectionDiseaseArea << "," << necroticDiseaseArea << ","
    //        << Utilities::formatDouble(cloudOValue) << "," << Utilities::formatDouble(cloudPValue) << "," 
    //        << Utilities::formatDouble(cloudFvalue)<<std::endl;
    BasicS::output.push_back(convert.str());

}

void PlantS::output() {
    std::ostringstream convert;
    //convert << "Cpp_Plant_" << getID() << ".txt";
    BasicS::getOutput("Cpp_Plant.txt", this->firstOutputCallS);
    this->firstOutputCallS++;

    // Speedup the model removing outputs
    //std::cout << "\nPlant " << getID() << ":\n";
    //for(unsigned int i=0; i<BasicS::output.size(); i++)
    //{
    //    std::cout << BasicS::output[i] << std::endl;
    //}

    OrganS *o;
    for (unsigned int i = 0; i < organs.size(); i++) {
        o = &organs[i];
        o->output();
    }

    CloudPS *cp;
    for (unsigned int i = 0; i < cloudsP.size(); i++) {
        cp = &cloudsP[i];
        cp->output();
    }
}

void PlantS::rate() {
    OrganS *o;
    for (unsigned int i = 0; i < organs.size(); i++) {
        o = &organs[i];
        if(o->getSenescenceArea() < o->getTotalArea()) {
            o->setProportionFromTotalArea(o->getTotalArea()/totalArea);
            o->rate();
        } else {
            //printf("## OrganS: %d died - rate!!!\n",o->getOrganNumber());
        }
    }
}
