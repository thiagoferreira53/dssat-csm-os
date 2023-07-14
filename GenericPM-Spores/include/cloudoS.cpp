#include "cloudoS.h"
#include<iostream>
#include<sstream>
#include<string>
//#include "../../FlexibleIO/Data/FlexibleIO.hpp"

//FlexibleIO *flexibleioS = FlexibleIO::getInstance();

int CloudOS::qtdS = 0;
int CloudOS::firstOutputCallS = 0;

void CloudOS::integration() {
    CloudS::integration();

    if (values.size() > (unsigned) disease->getVectorSizeCloudO()) {
        values.erase(values.begin());
    }
    CloudS::removeSporesCloudOByAgeS();


    if (getValue() > disease->getMaxSporeCloudsDensity()) {
        CloudS::removeSporesCloudS(getValue() - disease->getMaxSporeCloudsDensity()); 
    }
    if (BasicS::getWeather()->getRain() >= 20) {
        CloudS::removeSporesCloudByRainS(0.5);
    }

    std::ostringstream convert;
    convert << BasicS::getWeather()->getYearDoy() << "," << getValue();
    for (unsigned int i = 0; i < values.size(); i++) {
        convert << "," << values[i];
    }
    BasicS::output.push_back(convert.str());
}

void CloudOS::output() {
    CloudS::output();

    std::ostringstream convert;
    convert << "Cpp_CloudO_" << getID() << ".txt";
    BasicS::getOutput(convert.str(),this->firstOutputCallS);
    this->firstOutputCallS++;

    // Speedup the model removing outputs
    //std::cout << "\nCloudO" << ID << ":";
    //for(unsigned int i=0; i<BasicS::output.size(); i++)
    //    std::cout << BasicS::output[i] << std::endl;
}

void CloudOS::addSporesCreatedS(double sporesCreated) {
    CloudS::sporesCreated += (sporesCreated * (1 - disease->getProportionFromOrganToPlantCloud()));
    cloudP->addSporesCreatedS(sporesCreated * disease->getProportionFromOrganToPlantCloud());
}
