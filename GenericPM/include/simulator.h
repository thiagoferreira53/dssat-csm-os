#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "../../GenericPM/include/basicinterface.h"
#include "cropinterface.h"
#include "initialcondition.h"
#include "plant.h"


#include<vector>

class Simulator : virtual public BasicInterface {
private:
    Utilities util;

protected:
    Simulator();
    static Simulator *instance;
    int currentYearDoy = 0;
    CropInterface *cropinterface;
    std::vector<InitialCondition> initialConditions;
    std::vector<Plant> plants;

public:
    static Simulator* getInstance();
    static Simulator* newInstance();
    void inputPST();
    void inicialization();
    void integration();
    void output();
    void rate();
    void updateCurrentYearDoy(int yearDoy);
    bool allPlantsSenesced();

    std::vector<Plant>& getPlants() {
        return plants;
    }

    std::vector<InitialCondition>& getInitialConditions() {
        return initialConditions;
    }

    CropInterface* getCropInterface() {
        return cropinterface;
    }

    void setCurrentYearDoy(int currentYearDoy) {
        this->currentYearDoy = currentYearDoy;
    }

    int getCurrentYearDoy() const {
        return currentYearDoy;
    }
};

#endif // SIMULATOR_H
