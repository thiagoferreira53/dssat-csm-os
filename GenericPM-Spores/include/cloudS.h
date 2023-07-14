#ifndef CLOUDS_H
#define CLOUDS_H

#include "../../GenericPM/include/basic.h"
#include "../../GenericPM/include/basicinterface.h"
#include "diseaseS.h"
#include<vector>

class CloudS : public Basic, virtual public BasicInterface {
protected:
    std::vector<double> values;
    DiseaseS *disease;
    int sporesCreated = 0;
    int sporesToBeRemoved = 0;
    
public:

    DiseaseS* getDisease() {
        return disease;
    }
    virtual void addSporesCreatedS(double sporesCreated) = 0;

    void rate() {
    }
    void integration();

    void output() {
    }
    double getValue();
    void removeSporesCloudS(double toBeRemove);
    void removeSporesCloudByRainS(double percent);
    void removeSporesCloudFByAgeS(void);
    void removeSporesCloudPByAgeS(void);
    void removeSporesCloudOByAgeS(void);


    //void removeSporesCloudFByAgeUvS(void);

    int getSporesToBeRemoved() {
        return sporesToBeRemoved;
    }

    void setSporesToBeRemoved(int sporesToBeRemoved) {
        this->sporesToBeRemoved = sporesToBeRemoved;
    }

    void addSporesToBeRemoved(int sporesToBeRemoved) {
        this->sporesToBeRemoved += sporesToBeRemoved;
    }
};

#endif // CLOUD_H
