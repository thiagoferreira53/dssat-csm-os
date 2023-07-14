#ifndef CLOUDPS_H
#define CLOUDPS_H

#include "cloudS.h"
#include "../../GenericPM/include/basicinterface.h"
#include "cloudfS.h"

class CloudPS : public CloudS, virtual public BasicInterface {
private:
    CloudFS *cloudF;

protected:
    static int qtdS;
    int ID = ++qtdS;
    static int firstOutputCallS;

public:

    CloudPS(Disease *disease, CloudFS *cloudF) {
        this->disease = disease;
        this->cloudF = cloudF;
    }

    int getID() {
        return ID;
    }

    CloudFS* getCloudF() {
        return cloudF;
    }

    void rate() {
        CloudS::rate();
        cloudF->rate();
    }
    void integration();
    void output();
    void addSporesCreatedS(double sporesCreated);
};

#endif // CLOUDP_H
