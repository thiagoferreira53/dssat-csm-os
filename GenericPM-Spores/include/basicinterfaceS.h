#ifndef BASICINTERFACES_H
#define BASICINTERFACES_H

class BasicInterfaceS {
public:

    virtual ~BasicInterfaceS() {
    }
    virtual void rate() = 0;
    virtual void integration() = 0;
    virtual void output() = 0;

};

#endif // BASICINTERFACE_H
