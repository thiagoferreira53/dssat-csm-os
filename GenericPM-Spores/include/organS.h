#ifndef ORGANS_H
#define ORGANS_H

#include "../../GenericPM/include/basic.h"
#include "../../GenericPM/include/basicinterface.h"
#include "lesioncohortS.h"
#include "cloudoS.h"
#include "cloudpS.h"

#include <vector>

class OrganS : public Basic, virtual public BasicInterface {
protected:
    double totalArea = 0;
    double diseaseArea = 0;
    double visibleDiseaseArea = 0;
    double invisibleDiseaseArea = 0;
    double senescenceArea = 0;
    double latentDiseaseArea = 0;
    double infectionDiseaseArea = 0;
    double necroticDiseaseArea = 0;
    int organNumber = 0;
    int newLesionsS = 0;
    int newLesionsFromOrgan, newLesionsFromPlant, newLesionsFromField;
    int totalLesions = 0;
    int visibleLesions = 0;
    double dailyDiseaseArea = 0;
    double dailyVisibleDiseaseArea = 0;
    double dailySenescenceArea = 0;
    double dailyLatentDiseaseArea = 0;
    double dailyInfectionDiseaseArea = 0;
    double dailyNecroticDiseaseArea = 0;
    int dailyTotalLesions = 0;
    int dailyVisibleLesions = 0;
    bool suceptible = false;
    double organLastSize = 0;
    double healthAreaProportion = 0;
    double proportionFromTotalArea = 0;
    int doc = Basic::getWeather()->getDoy();
    double physiologicalLife = 0;
    std::vector<LesionCohortS> lesionCohorts;
    std::vector<CloudOS> cloudsO;
    static int firstOutputCallS;

public:

    OrganS(std::vector<CloudPS>& cloudsP, int organNumber, double totalArea) {
        //Basic::output.push_back("OrganS, YearDoy, TotalArea, Senesced, Diseased, VisibleArea, InvisibleArea, LesionDensity, Age, newLesionsS, TotalLesions, CloudOS, CloudPS, CloudFS, HealthAreaProportion");
        this->organNumber = organNumber;
        this->totalArea = totalArea;
        CloudPS *cloud;
        for (unsigned int i = 0; i < cloudsP.size(); i++) {
            cloud = &cloudsP[i];
            this->cloudsO.emplace_back(cloud->getDisease(), cloud);
        }
    }

    void integration();
    void cloudIntegrationS();
    float cloudAmountS();
    void output();
    void rate();

    float getAreaIfSeverity(int porcent) {
        if(diseaseArea>(totalArea*(porcent/100.0))) {
            return totalArea;
        }
        return 0.0;
    }
    bool isAlive() {
        if (this->totalArea > 0) {//this->senescenceArea) 
            return true;
        }
        else {
            return false;
        }
    }

    double getDensity() {
        return totalArea > 0 ? (visibleLesions / totalArea) : 0;
    }

    std::vector<LesionCohortS>& getlesionCohorts() {
        return lesionCohorts;
    }

    std::vector<CloudOS>& getCloudsO() {
        return cloudsO;
    }

    void setLesionCohorts(std::vector<LesionCohortS> lesionCohorts) {
        this->lesionCohorts = lesionCohorts;
    }

    std::vector<LesionCohortS> getLesionCohorts() const {
        return lesionCohorts;
    }

    void setPhysiologicalLife(double physiologicalLife) {
        this->physiologicalLife = physiologicalLife;
    }

    double getPhysiologicalLife() const {
        return physiologicalLife;
    }

    void setDoc(int doc) {
        this->doc = doc;
    }
    
    double getProportionFromTotalArea() {
        return proportionFromTotalArea;
    }

    void setProportionFromTotalArea(double proportionFromTotalArea) {
        this->proportionFromTotalArea = proportionFromTotalArea;
    }

    int getDoc() const {
        return doc;
    }

    void setHealthAreaProportion(double healthAreaProportion) {
        this->healthAreaProportion = healthAreaProportion;
    }

    double getHealthAreaProportion() const {
        return healthAreaProportion;
    }

    void setOrganLastSize(double organLastSize) {
        this->organLastSize = organLastSize;
    }

    double getOrganLastSize() const {
        return organLastSize;
    }

    void setSuceptible(bool suceptible) {
        this->suceptible = suceptible;
    }

    bool isSuceptible() const {
        return suceptible;
    }

    void setDailyVisibleLesions(int dailyVisibleLesions) {
        this->dailyVisibleLesions = dailyVisibleLesions;
    }

    int getDailyVisibleLesions() const {
        return dailyVisibleLesions;
    }

    void setDailyTotalLesions(int dailyTotalLesions) {
        this->dailyTotalLesions = dailyTotalLesions;
    }

    int getDailyTotalLesions() const {
        return dailyTotalLesions;
    }

    void setDailyNecroticDiseaseArea(double dailyNecroticDiseaseArea) {
        this->dailyNecroticDiseaseArea = dailyNecroticDiseaseArea;
    }

    double getDailyNecroticDiseaseArea() const {
        return dailyNecroticDiseaseArea;
    }

    void setDailyInfectionDiseaseArea(double dailyInfectionDiseaseArea) {
        this->dailyInfectionDiseaseArea = dailyInfectionDiseaseArea;
    }

    double getDailyInfectionDiseaseArea() const {
        return dailyInfectionDiseaseArea;
    }

    void setDailyLatentDiseaseArea(double dailyLatentDiseaseArea) {
        this->dailyLatentDiseaseArea = dailyLatentDiseaseArea;
    }

    double getDailyLatentDiseaseArea() const {
        return dailyLatentDiseaseArea;
    }

    void setDailySenescenceArea(double dailySenescenceArea) {
        this->dailySenescenceArea = dailySenescenceArea;
    }

    double getDailySenescenceArea() const {
        return dailySenescenceArea;
    }

    void setDailyVisibleDiseaseArea(double dailyVisibleDiseaseArea) {
        this->dailyVisibleDiseaseArea = dailyVisibleDiseaseArea;
    }

    double getDailyVisibleDiseaseArea() const {
        return dailyVisibleDiseaseArea;
    }

    void setDailyDiseaseArea(double dailyDiseaseArea) {
        this->dailyDiseaseArea = dailyDiseaseArea;
    }

    double getDailyDiseaseArea() const {
        return dailyDiseaseArea;
    }

    void setVisibleLesions(int visibleLesions) {
        this->visibleLesions = visibleLesions;
    }

    int getVisibleLesions() const {
        return visibleLesions;
    }

    void setTotalLesions(int totalLesions) {
        this->totalLesions = totalLesions;
    }

    int getTotalLesions() const {
        return totalLesions;
    }

    void setNewLesions(int newLesionsS) {
        this->newLesionsS = newLesionsS;
    }

    int getNewLesions() const {
        return newLesionsS;
    }

    void setOrganNumber(int organNumber) {
        this->organNumber = organNumber;
    }

    int getOrganNumber() const {
        return organNumber;
    }

    void setNecroticDiseaseArea(double necroticDiseaseArea) {
        this->necroticDiseaseArea = necroticDiseaseArea;
    }

    double getNecroticDiseaseArea() const {
        return necroticDiseaseArea;
    }

    void setInfectionDiseaseArea(double infectionDiseaseArea) {
        this->infectionDiseaseArea = infectionDiseaseArea;
    }

    double getInfectionDiseaseArea() const {
        return infectionDiseaseArea;
    }

    void setLatentDiseaseArea(double latentDiseaseArea) {
        this->latentDiseaseArea = latentDiseaseArea;
    }

    double getLatentDiseaseArea() const {
        return latentDiseaseArea;
    }

    void setSenescenceArea(double senescenceArea) {
        this->senescenceArea = senescenceArea;
    }

    double getSenescenceArea() const {
        return senescenceArea;
    }

    void setVisibleDiseaseArea(double visibleDiseaseArea) {
        this->visibleDiseaseArea = visibleDiseaseArea;
    }

    double getVisibleDiseaseArea() const {
        return visibleDiseaseArea;
    }

    double getInvisibleDiseaseArea() const {
        return invisibleDiseaseArea;
    }

    void setDiseaseArea(double diseaseArea) {
        this->diseaseArea = diseaseArea;
    }

    double getDiseaseArea() const {
        return diseaseArea;
    }

    void setTotalArea(double totalArea) {
        this->totalArea = totalArea;
    }

    double getTotalArea() const {
        return totalArea;
    }

    int getNewLesionsFromOrgan() {
        return newLesionsFromOrgan;
    }

    void setNewLesionsFromOrgan(int newLesionsFromOrgan) {
        this->newLesionsFromOrgan = newLesionsFromOrgan;
    }

    int getNewLesionsFromPlant() {
        return newLesionsFromPlant;
    }

    void setNewLesionsFromPlant(int newLesionsFromPlant) {
        this->newLesionsFromPlant = newLesionsFromPlant;
    }

    int getNewLesionsFromField() {
        return newLesionsFromField;
    }

    void setNewLesionsFromField(int newLesionsFromField) {
        this->newLesionsFromField = newLesionsFromField;
    }

};

#endif // ORGAN_H
