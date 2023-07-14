#ifndef INITIALCONDITIONS_H
#define INITIALCONDITIONS_H

#include "basicS.h"
#include "basicinterfaceS.h"
#include "cloudfS.h"

class InitialConditionS : public BasicS, virtual public BasicInterfaceS {
private:
    CloudFS cloudf;

protected:
    double acumulateFavorability = 0, dailyFavorability = 0;
    bool stop = false;
    int doc = BasicS::getWeather()->getDoy();
    static int qtdS;
    int ID = ++qtdS;

public:
    InitialConditionS(DiseaseS *disease) : cloudf{disease}
    {
        BasicS::output.push_back("Day of Simulation, Acumulated Favorability");
    }

    int getID() {
        return ID;
    }

    CloudFS* getCloud() {
        return &cloudf;
    }

    int getDoc() {
        return doc;
    }
    void rate();
    void integration();
    void output();

};

#endif // INITIALCONDITION_H
