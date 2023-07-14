#ifndef CLOUD_H
#define CLOUD_H

#include "basic.h"
#include "../../GenericPM/include/basicinterface.h"
#include "disease.h"
#include<vector>

class Cloud : public Basic, virtual public BasicInterface {
protected:
    std::vector<double> values;
    Disease *disease;
    int sporesCreated = 0;
    int sporesToBeRemoved = 0;
    
public:

    Disease* getDisease() {
        return disease;
    }
    virtual void addSporesCreated(double sporesCreated) = 0;

    void rate() {
    }
    void integration();

    void output() {
    }
    double getValue();
    void removeSporesCloud(double toBeRemove);
    void removeSporesCloudByRain(double percent);
    void removeSporesCloudFByAge(void);
    void removeSporesCloudPByAge(void);
    void removeSporesCloudOByAge(void);



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
