#ifndef PLANT_H
#define PLANT_H

#include "basic.h"
#include "../../GenericPM/include/basicinterface.h"
#include "cloudp.h"
#include "organ.h"

#include <vector>

class Plant : public Basic, virtual public BasicInterface {
protected:
    int doc = Basic::getWeather()->getDoy();
    std::vector<Organ> organs;
    std::vector<CloudP> cloudsP;
    static int qtd;
    int ID = ++qtd;
    double totalArea = 0, diseaseArea = 0, latentDiseaseArea = 0, infectionDiseaseArea = 0,
            necroticDiseaseArea = 0, visibleDiseaseArea = 0, invisibleDiseaseArea = 0, visibleLesions = 0,
            totalLesions = 0, senescenceArea = 0;
    static int firstOutputCall;
public:
    Plant();

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

    std::vector<Organ>& getOrgans() {
        return organs;
    }

    std::vector<CloudP>& getCloudsP() {
        return cloudsP;
    }

    bool isAlive() {
        if (organs.size() == 0)
            return true;
        for (int i = organs.size() - 1; i >= 0; i--) {
            Organ *o = &organs[i];
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
