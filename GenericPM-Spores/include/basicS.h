#ifndef BASICS_H
#define BASICS_H

#include "../../GenericPM/include/weather.h"

#include<vector>
#include<string>
#include<fstream>

class BasicS {
protected:
    std::vector<std::string> output;
    Weather *weather;
    int firstCall = 1;

public:

    void getOutput(std::string fileName) {
        std::ofstream out;
//        if(firstCall) {
//            out.open(fileName);    
//            firstCall = 0;
//        } else {
            out.open(fileName, std::ofstream::out | std::ofstream::app);
//        }
        //out << fileName << std::endl;
        for (unsigned int i = 0; i < output.size(); i++)
            out << output[i] << std::endl;
        out.close();
        output.clear();
    }

    void getOutput(std::string fileName, int updateFile) {
        std::ofstream out;
        if(updateFile) {
            out.open(fileName, std::ofstream::out | std::ofstream::app);            
        } else {
            out.open(fileName);    
        }
        //out << fileName << std::endl;
        for (unsigned int i = 0; i < output.size(); i++)
            //out << output[i] << " " << BasicS::getWeather() << std::endl;
            out << output[i] << std::endl;
        out.close();
        output.clear();
    }

    Weather* getWeather() {
        return weather = Weather::getInstance();
    }

};

#endif // BASIC_H
