#ifndef CLOUDOS_H
#define CLOUDOS_H

#include "cloudS.h"
#include "basicinterfaceS.h"
#include "cloudpS.h"

class CloudOS : public CloudS, virtual public BasicInterfaceS {
private:
    CloudPS *cloudP;

protected:
    static int qtdS;
    int ID = ++qtdS;
    static int firstOutputCallS;

public:

    CloudOS(DiseaseS *disease, CloudPS *cloudP) {
        this->disease = disease;
        this->cloudP = cloudP;
    }

    int getID() {
        return ID;
    }

    void rate() {
        CloudS::rate();
        cloudP->rate();
    }

    CloudPS* getCloudP() {
        return cloudP;
    }
    void integration();
    void output();
    void addSporesCreatedS(double sporesCreated);


};

#endif // CLOUDO_H
