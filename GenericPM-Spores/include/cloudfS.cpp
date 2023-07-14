#include "cloudfS.h"
#include<iostream>
#include<sstream>
//#include "../../FlexibleIO/Data/FlexibleIO.hpp"

//FlexibleIO *flexibleioS = FlexibleIO::getInstance();

int CloudFS::qtdS = 0;
int CloudFS::firstOutputCallS = 0;

void CloudFS::integration() {
    CloudS::integration();
    float porcent=0;
    if (values.size() > (unsigned) disease->getVectorSizeCloudF()) {
        values.erase(values.begin());
    }

    if (getValue() > disease->getMaxSporeCloudsDensity()) {
        CloudS::removeSporesCloudS(getValue() - disease->getMaxSporeCloudsDensity());
    }
    if (Basic::getWeather()->getRain() >= 20) {
        porcent = Basic::getWeather()->getRain() / 80;
        porcent = porcent>1?1:porcent;
        //std::cout<< "porcent : "<<porcent << std::endl;
        CloudS::removeSporesCloudByRainS(porcent);

        
    }

    std::ostringstream convert;
    convert << Basic::getWeather()->getYearDoy() << "," << getValue();
    for (unsigned int i = 0; i < values.size(); i++) {
        convert << "," << values[i];
    }
    Basic::output.push_back(convert.str());
}

void CloudFS::output() {
    CloudS::output();

    std::ostringstream convert;
    convert << "Cpp_CloudF_" << getID() << ".txt";
    Basic::getOutput(convert.str(),this->firstOutputCallS);
    this->firstOutputCallS++;

    // Speedup the model
    //std::cout << "\nCloudF" << ID << ":";
    //for(unsigned int i=0; i<Basic::output.size(); i++)
    //    std::cout << Basic::output[i] << std::endl;
}

