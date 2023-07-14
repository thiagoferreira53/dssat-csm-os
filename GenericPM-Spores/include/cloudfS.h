#ifndef CLOUDFS_H
#define CLOUDFS_H

#include "cloudS.h"
#include "basicinterfaceS.h"
#include "diseaseS.h"
#include "cloudS.h"

class CloudFS : public CloudS, virtual public BasicInterfaceS {
protected:
    static int qtdS;
    int ID = ++qtdS;
    double firstSporeCloud = 0;
    static int firstOutputCallS;

public:

    CloudFS(DiseaseS *disease) {
        this->disease = disease;
    }

    int getID() {
        return ID;
    }
    void integration();

    double getValue() {
        return CloudS::getValue() + firstSporeCloud;
    }
    void setValue(double value) {
        CloudS::values.clear();
        CloudS::values.push_back(value);
    }
    void output();

    void rate() {
        CloudS::rate();
    }

    void addSporesCreatedS(double sporesCreated) {
        this->sporesCreated += sporesCreated;
    }

    void setFirstSporeCloud(double firstSporeCloud) {
        this->firstSporeCloud = firstSporeCloud;
    }

};

#endif // CLOUDF_H
