#ifndef CLOUDO_H
#define CLOUDO_H

#include "cloud.h"
#include "../../GenericPM/include/basicinterface.h"
#include "cloudp.h"

class CloudO : public Cloud, virtual public BasicInterface {
private:
    CloudP *cloudP;

protected:
    static int qtd;
    int ID = ++qtd;
    static int firstOutputCall;

public:

    CloudO(Disease *disease, CloudP *cloudP) {
        this->disease = disease;
        this->cloudP = cloudP;
    }

    int getID() {
        return ID;
    }

    void rate() {
        Cloud::rate();
        cloudP->rate();
    }

    CloudP* getCloudP() {
        return cloudP;
    }
    void integration();
    void output();
    void addSporesCreated(double sporesCreated);


};

#endif // CLOUDO_H
