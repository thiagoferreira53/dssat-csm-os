#ifndef SIMULATORS_H
#define SIMULATORS_H

#include "basicinterfaceS.h"
#include "cropinterfaceS.h"
#include "initialconditionS.h"
#include "plantS.h"


#include<vector>

class SimulatorSpore : virtual public BasicInterfaceS {
private:
    Utilities util;

protected:
    SimulatorSpore();
    static SimulatorSpore *instance;
    int currentYearDoy = 0;

    CropInterfaceS *cropinterface;
    std::vector<InitialConditionS> initialConditions;
    std::vector<PlantS> plants;

public:
    static SimulatorSpore* getInstance();
    static SimulatorSpore* newInstance();
    void inputPSTS();
    void inicializationS();
    void integration();
    void output();
    void rate();
    void updateCurrentYearDoy(int yearDoy);
    bool allPlantsSenescedS();

    std::vector<PlantS>& getPlants() {
        return plants;
    }

    std::vector<InitialConditionS>& getInitialConditions() {
        return initialConditions;
    }

    CropInterfaceS* getCropInterface() {
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
