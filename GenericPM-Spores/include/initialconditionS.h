#ifndef INITIALCONDITIONS_H
#define INITIALCONDITIONS_H

#include "../../GenericPM/include/basic.h"
#include "../../GenericPM/include/basicinterface.h"
#include "cloudfS.h"

class InitialConditionS : public Basic, virtual public BasicInterface {
private:
    CloudFS cloudf;

protected:
    double acumulateFavorability = 0, dailyFavorability = 0;
    bool stop = false;
    int doc = Basic::getWeather()->getDoy();
    static int qtdS;
    int ID = ++qtdS;

public:
    InitialConditionS(Disease *disease) : cloudf{disease}
    {
        Basic::output.push_back("Day of Simulation, Acumulated Favorability");
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
