#ifndef CLOUDF_H
#define CLOUDF_H

#include "cloud.h"
#include "../../GenericPM/include/basicinterface.h"
#include "disease.h"
#include "cloud.h"

class CloudF : public Cloud, virtual public BasicInterface {
protected:
    static int qtd;
    int ID = ++qtd;
    double firstSporeCloud = 0;
    static int firstOutputCall;

public:

    CloudF(Disease *disease) {
        this->disease = disease;
    }

    int getID() {
        return ID;
    }
    void integration();

    double getValue() {
        return Cloud::getValue() + firstSporeCloud;
    }
    void output();

    void rate() {
        Cloud::rate();
    }

    void addSporesCreated(double sporesCreated) {
        this->sporesCreated = sporesCreated + Cloud::getValue();
    }

    void setFirstSporeCloud(double firstSporeCloud) {
        this->firstSporeCloud = firstSporeCloud;
    }

};

#endif // CLOUDF_H
