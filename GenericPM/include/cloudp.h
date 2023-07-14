#ifndef CLOUDP_H
#define CLOUDP_H

#include "cloud.h"
#include "../../GenericPM/include/basicinterface.h"
#include "cloudf.h"

class CloudP : public Cloud, virtual public BasicInterface {
private:
    CloudF *cloudF;

protected:
    static int qtd;
    int ID = ++qtd;
    static int firstOutputCall;

public:

    CloudP(Disease *disease, CloudF *cloudF) {
        this->disease = disease;
        this->cloudF = cloudF;
    }

    int getID() {
        return ID;
    }

    CloudF* getCloudF() {
        return cloudF;
    }

    void rate() {
        Cloud::rate();
        cloudF->rate();
    }
    void integration();
    void output();
    void addSporesCreated(double sporesCreated);
};

#endif // CLOUDP_H
