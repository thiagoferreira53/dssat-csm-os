#ifndef PLANTS_H
#define PLANTS_H

#include "../../GenericPM/include/basic.h"
#include "../../GenericPM/include/basicinterface.h"
#include "cloudpS.h"
#include "organS.h"

#include <vector>

class PlantS : public Basic, virtual public BasicInterface {
protected:
    int doc = Basic::getWeather()->getDoy();
    std::vector<OrganS> organs;
    std::vector<CloudPS> cloudsP;
    static int qtdS;
    int ID = ++qtdS;
    double totalArea = 0, diseaseArea = 0, latentDiseaseArea = 0, infectionDiseaseArea = 0,
            necroticDiseaseArea = 0, visibleDiseaseArea = 0, invisibleDiseaseArea = 0, visibleLesions = 0,
            totalLesions = 0, senescenceArea = 0;
    static int firstOutputCallS;
public:
    PlantS();

    int getID() {
        return ID;
    }
    void integration();
    void output();
    void rate();

    double getTotalArea() {
        return totalArea;
    }

    int getDoc() {
        return doc;
    }

    std::vector<OrganS>& getOrgans() {
        return organs;
    }

    std::vector<CloudPS>& getCloudsP() {
        return cloudsP;
    }

    bool isAlive() {
        if (organs.size() == 0)
            return true;
        for (int i = organs.size() - 1; i >= 0; i--) {
            OrganS *o = &organs[i];
            if (o->isAlive())
                return true;
        }
        return false;
    }

    double getDiseaseArea() {
        return diseaseArea;
    }

    double getLatentDiseaseArea() {
        return latentDiseaseArea;
    }

    double getInfectionDiseaseArea() {
        return infectionDiseaseArea;
    }

    double getNecroticDiseaseArea() {
        return necroticDiseaseArea;
    }

    double getVisibleDiseaseArea() {
        return visibleDiseaseArea;
    }

    double getInvisibleDiseaseArea() {
        return invisibleDiseaseArea;
    }

    double getVisibleLesions() {
        return visibleLesions;
    }

    double getTotalLesions() {
        return totalLesions;
    }

    double getSenescenceArea() {
        return senescenceArea;
    }
    
};

#endif // PLANT_H
