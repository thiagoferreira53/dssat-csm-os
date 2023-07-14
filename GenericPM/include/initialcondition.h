#ifndef INITIALCONDITION_H
#define INITIALCONDITION_H

#include "basic.h"
#include "../../GenericPM/include/basicinterface.h"
#include "cloudf.h"

class InitialCondition : public Basic, virtual public BasicInterface {
private:
    CloudF cloudf;

protected:
    double acumulateFavorability = 0, dailyFavorability = 0;
    bool stop = false;
    int doc = Basic::getWeather()->getDoy();
    static int qtd;
    int ID = ++qtd;

public:
    InitialCondition(Disease *disease) : cloudf{disease}
    {
        Basic::output.push_back("Day of Simulation, Acumulated Favorability");
    }

    int getID() {
        return ID;
    }

    CloudF* getCloud() {
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
