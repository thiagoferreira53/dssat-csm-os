/**
 * @file utilities.cpp
 * 
 * @author Willingthon Pavan (wpavan.us@gmail.com)
 * @author Jose Mauricio Cunha Fernandes (jmauricio.fernandes@icloud.com)
 * 
 * @copyright Copyright (c) 2017–2025, DSSAT Foundation
 * @license BSD-3-Clause. See the LICENSE file in the root folder for details.
 */
#include "utilities.h"
#include<cmath>
#include<string>
#include<sstream>
#include<iomanip>
#include "disease.h"
#include<iostream>
#include "../../FlexibleIO/Data/FlexibleIO.hpp"


double Utilities::trapezoidalFunction(double value, double v[]) {
    if (value >= v[1] && value <= v[2])
        return 1;
    else if (value >= v[0] && value < v[1])
        return ((value - v[0]) / (v[1] - v[0]));
    else if (value > v[2] && value <= v[3])
        return 1 - ((value - v[2]) / (v[3] - v[2]));
    else
        return 0;
}

double Utilities::temperatureFavorability(double temp, double cardinalTemperatures[]) {
    double tempMax = cardinalTemperatures[0], 
           tempMin = cardinalTemperatures[1], 
           tempOpt = cardinalTemperatures[2];
    double tf, a, b;

    b = ((tempMax - tempOpt) / (tempOpt - tempMin));
    a = (1 / ((tempOpt - tempMin) * pow(tempMax - tempOpt, b)));
    tf = (a * (temp - tempMin) * pow(tempMax - temp, b));
    return (fmax(0,tf));
}

double Utilities::wetnessFavorability(double wetDuration, std::string wetnessFunction) {
//    double wf;
    double wf, k = 0.99374, b = 350.81064, r = 0.88591;
    //wf = (k / (1 + b * exp(-r * wetDuration)));
    wf = Utilities::growthFunction(wetnessFunction, wetDuration);
    //std::cout<<"WF: "<<wf<<std::endl;
    //std::cout<<"WF: "<<wf<<" wetnessFunction: "<<wetnessFunction<<" wetDuration: "<<wetDuration<<std::endl;

    //printf("WF: %f Wetness Function: %s WetDuration: %f",wetnessFunction, wf,wetDuration);  
    return wf;
}

double Utilities::dispersalRainFunction(double rain, std::string dispersionFreequency) {
    double dre;
    dre = Utilities::growthFunction(dispersionFreequency, rain);
    //std::cout<<"DRE: "<<dre<<" Rain: "<<rain<<" DisperifonF: "<<dispersionFreequency<<std::endl;
    return dre;
}
double Utilities::runExpressionFunction(double value, std::string expressionString) {
    return Utilities::runExpression(expressionString, value);
}

double Utilities::getHealthAreaProportion(double diseaseArea, double totalArea, double senescedArea) {
        double hap = 1 - ((diseaseArea+senescedArea) / totalArea);
        return (hap > 0 ? hap : 0);
    //double hap = (1 - (diseaseArea / totalArea));
    //return (hap > 0 ? hap : 0);
}

std::string Utilities::formatDouble(double value) {
    return formatDouble(value, 2);
}

std::string Utilities::formatDouble(double value, int decimals) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(decimals) << value;
    return ss.str();
}

int Utilities::addOneDay(int yearDoy) {
    yearDoy++;
    int year = yearDoy / 1000;
    int doy = yearDoy - (year * 1000);
    if (doy <= 365 || (doy == 366 && isLeapYear(year))) {
        return yearDoy;
    } else {
        return (1 + ((year + 1)*1000));
    }
}

bool Utilities::isLeapYear(int year) {
    return year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
}

double Utilities::temperatureFactor(double TAVG) {
    // Read cardinal temperature parameters from WHGEN048.PST via FlexibleIO
    // TFS field format: TFS = max, min, opt (e.g., "32,10,23")
    FlexibleIO *fio = FlexibleIO::getInstance();
    
    double tmax = fio->getReal("PST", "TFS1");   // Maximum temperature
    double tmin = fio->getReal("PST", "TFS2");   // Minimum temperature
    double topt = fio->getReal("PST", "TFS3");   // Optimal temperature
    
    // Check for invalid or non-finite inputs
    if (!std::isfinite(TAVG) || !std::isfinite(tmin) || !std::isfinite(topt) || !std::isfinite(tmax))
        return 0.0;

    // Degenerate parameter case
    if (tmin >= topt || topt >= tmax)
        return 1.0;

    if (TAVG <= tmin || TAVG >= tmax)
        return 0.0;

    // Avoid division by zero
    double denom1 = std::max(topt - tmin, 1e-9);
    double denom2 = std::max(tmax - topt, 1e-9);

    double a = denom1 / denom2;
    double b = denom2 / denom1;

    double x1 = std::max((TAVG - tmin) / denom1, 1e-12);
    double x2 = std::max((tmax - TAVG) / denom2, 1e-12);

    double f = std::pow(x1, a) * std::pow(x2, b);

    return std::min(std::max(f, 0.0), 1.0);
}