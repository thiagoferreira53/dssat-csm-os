#include "cloudpS.h"

#include<iostream>
#include<sstream>
//#include "../../FlexibleIO/Data/FlexibleIO.hpp"

//FlexibleIO *flexibleioS = FlexibleIO::getInstance();
int CloudPS::qtdS = 0;
int CloudPS::firstOutputCallS = 0;

void CloudPS::integration() {
    CloudS::integration();

    //std::cout << values.size()<<std::endl;
    if (values.size() > (unsigned) disease->getVectorSizeCloudP()) {
        values.erase(values.begin());
    }
    CloudS::removeSporesCloudPByAgeS();

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

void CloudPS::output() {
    CloudS::output();

    std::ostringstream convert;
    convert << "Cpp_CloudP_" << getID() << ".txt";
    BasicS::getOutput(convert.str(),this->firstOutputCallS);
    this->firstOutputCallS++;

    // Speedup the model removing outputs
    //std::cout << "\nCloudP" << getID() << ":";
    //for(unsigned int i=0; i<BasicS::output.size(); i++)
    //    std::cout << BasicS::output[i] << std::endl;
}

void CloudPS::addSporesCreatedS(double sporesCreated) {
    CloudS::sporesCreated += (sporesCreated * (1 - disease->getProportionFromPlantToFieldCloud()));
    cloudF->addSporesCreatedS(sporesCreated * disease->getProportionFromPlantToFieldCloud());
}
