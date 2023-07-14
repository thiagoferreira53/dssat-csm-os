#ifndef CINTERFACES_H
#define CINTERFACES_H

#include <iostream>
#include "include/simulatorS.h"
#include "../GenericPM/include/utilities.h"
#include <cmath>
#include <vector>
#include <fstream>
#include "../FlexibleIO/Data/FlexibleIO.hpp"
#include "json.hpp"


class CinterfaceSpore {

private:
int YRDOY;
double CloudValue = 0;

public:
    double couplingIntegrationSpore(int YRDOY);

    void setcouplingCloudSpore(double CloudValue) {
        //std::cout<<"CloudValue H: "<<CloudValue<<std::endl;
        this->CloudValue = CloudValue;
    }

    double getcouplingCloudSpore() {
        return CloudValue;
    }

};

#endif // CINTERFACES_H


