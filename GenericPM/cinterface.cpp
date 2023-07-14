#include "include/simulator.h"
#include<cmath>
#include "../GenericPM-Spores/cinterfaceS.h"
using namespace std;

extern "C" {
    // Coupling Functions
    int couplingInit(int *YRDOY, int *YRPLT);
    int couplingRate(int *YRDOY,
            float *AREALF, float *CLW, float *CSW, float *PCLMT, float *PCSTMD,
            float *PDLA, float *PLFAD, float *PLFMD, float *PSTMD, float *PVSTGD,
            float *SLA, float *SLDOT, float *SSDOT, float *STMWT, float *TDLA,
            float *VSTGD, float *WLFDOT, float *WSTMD, float *WTLF,
            float *TLFAD, float *TLFMD, float *VSTAGE, float *WLIDOT,
            float *CLAI, float *CLFM, float *CSTEM, float *DISLA, float *DISLAP,
            float *LAIDOT, float *WSIDOT, float *SDWT, float *WSDD, 
            float *PSDD, int *DAS);
    int couplingIntegration(int *YRDOY,
            float *AREALF, float *CLW, float *CSW, float *PCLMT, float *PCSTMD,
            float *PDLA, float *PLFAD, float *PLFMD, float *PSTMD, float *PVSTGD,
            float *SLA, float *SLDOT, float *SSDOT, float *STMWT, float *TDLA,
            float *VSTGD, float *WLFDOT, float *WSTMD, float *WTLF,
            float *TLFAD, float *TLFMD, float *VSTAGE, float *WLIDOT,
            float *CLAI, float *CLFM, float *CSTEM, float *DISLA, float *DISLAP,
            float *LAIDOT, float *WSIDOT, float *SDWT, float *WSDD, 
            float *PSDD, int *DAS);
    int couplingOutput(int *doy);

}

float CLWp, SLAp, SDWTp;

// Coupling Functions Implementation 

int couplingInit(int *YRDOY, int *YRPLT) {
    // Set the start day for Disease Model
    // Get an instance of Simulator
    Simulator *s = Simulator::newInstance();
    // Set the start day for Disease Model
    s->setCurrentYearDoy(*YRDOY);
    // Set the sowing/planting date
    s->getCropInterface()->setPlantingDate(*YRPLT);

//    printf("YRDOY: %i YRPLT: %i\n",
//            *YRDOY, *YRPLT);
    CLWp=0; SLAp=0; SDWTp=0;
    return (1);
}

int couplingRate(int *YRDOY,
        float *AREALF, float *CLW, float *CSW, float *PCLMT, float *PCSTMD,
        float *PDLA, float *PLFAD, float *PLFMD, float *PSTMD, float *PVSTGD,
        float *SLA, float *SLDOT, float *SSDOT, float *STMWT, float *TDLA,
        float *VSTGD, float *WLFDOT, float *WSTMD, float *WTLF,
        float *TLFAD, float *TLFMD, float *VSTAGE, float *WLIDOT,
        float *CLAI, float *CLFM, float *CSTEM, float *DISLA, float *DISLAP,
        float *LAIDOT, float *WSIDOT, float *SDWT, 
        float *WSDD, float *PSDD, int *DAS) {
    // Temporary variable used for computations 
    float temp = 0, newOrgan = 0;
    double CloudField = 0;
    // Get an instance of Simulator
    Simulator *s = Simulator::getInstance();
    SimulatorSpore *sS = SimulatorSpore::getInstance();
    newOrgan = s->getCropInterface()->getOrgansQtd()+1;
    // Set the current YearDOY for next Disease step computation
    s->updateCurrentYearDoy(*YRDOY);
    
    CloudField = sS->getPlants()[0].getCloudsP()[0].getCloudF()->getValue();
    
    //std::cout<<"sS->getPlants()[0].getCloudsP()[0].getCloudF()->getValue(); "<<sS->getPlants()[0].getCloudsP()[0].getCloudF()->getValue()<<std::endl;

    //CinterfaceSpore spores;
    ////CloudField = spores.couplingIntegrationSpore(*YRDOY);
    //spores.couplingIntegrationSpore(*YRDOY);
    //CloudField = spores.getcouplingCloudSpore();
    // Set the current Leaf area for a specific organ (one big leaf for awhile)
    if(*SDWT-SDWTp > 0){ //&& CloudField > 0){
        s->getCropInterface()->setOrganArea(newOrgan, (*SDWT-SDWTp));
        //printf("SDWT %f SDWTp %f *SDWT-SDWTp %f\n", *SDWT, SDWTp, *SDWT-SDWTp);
        if(s->getPlants().size()>0) {
        //std::cout <<"Spores antes: "<<s->getPlants()[0].getCloudsP()[0].getCloudF()->getValue()<<std::endl;
        s->getPlants()[0].getCloudsP()[0].getCloudF()->addSporesCreated(CloudField);
        std::cout << "Data "<< *YRDOY<< " Alterando SPORES de: "<< s->getPlants()[0].getCloudsP()[0].getCloudF()->getValue() << " para: " << CloudField <<std::endl;
        }
    }


    SDWTp = *SDWT;

    // Feed the Disease Model with weather information
    Weather::getInstance()->update();
    // Disease Simulator Rate
    s->rate();
    //printf("Returning rate\n");
    return (1);
}

int couplingIntegration(int *YRDOY,
        float *AREALF, float *CLW, float *CSW, float *PCLMT, float *PCSTMD,
        float *PDLA, float *PLFAD, float *PLFMD, float *PSTMD, float *PVSTGD,
        float *SLA, float *SLDOT, float *SSDOT, float *STMWT, float *TDLA,
        float *VSTGD, float *WLFDOT, float *WSTMD, float *WTLF,
        float *TLFAD, float *TLFMD, float *VSTAGE, float *WLIDOT,
        float *CLAI, float *CLFM, float *CSTEM, float *DISLA, float *DISLAP,
        float *LAIDOT, float *WSIDOT, float *SDWT, 
        float *WSDD, float *PSDD, int *DAS) {
    // Temporary variable used for computations 
    float dArea = 0,tArea=0,sArea=0; //,pclaCalc=0,pDArea=0;
     int seedAge = 0;
    
    // Get an instance of Simulator
    Simulator *s = Simulator::getInstance();


    // Call the Disease Model Integration function
    s->integration();
    //printf("YRDOY: %i Plant size: %f Organ size: %f\n", *YRDOY, s->getPlants().size(),s->getPlants()[0].getOrgans().size());
    //std::cout<<"s->getPlants().size() "<<s->getPlants().size()<<std::endl;
    //std::cout<<s->getPlants()[0].getOrgans().size()<<std::endl;
    if(s->getPlants().size() > 0 && s->getPlants()[0].getOrgans().size() > 0) {
        dArea = s->getPlants()[0].getDiseaseArea();
        tArea = s->getPlants()[0].getTotalArea();
        sArea = s->getPlants()[0].getSenescenceArea();
        seedAge = s->getPlants()[0].getOrgans().size();
        //pDArea = (dArea/(tArea-sArea)*100);

        //*PSDD = (dArea/tArea*5);
        *PSDD= ((dArea/tArea)*15);
        //printf("ORIGINAL YRDOY: %i CloudF: %f PSDD %f\n",*YRDOY, s->getPlants()[0].getCloudsP()[0].getCloudF()->getValue(), PSDD);
        //printf("YRDOY: %i Plant Total Area: %f Disease Area: %f Senescence Area: %f AREALF: %f PDLA: %f PLFAD: %f\n", *YRDOY, tArea,dArea,sArea,*AREALF,*PDLA,*PLFAD);
        //printf("YRDOY: %i SDWT: %f PSDD: %f\n", *YRDOY, *SDWT, *PSDD);


    }
    
    //printf("YRDOY: %i PCLMT: %f\n", *YRDOY,*PCLMT);
    //printf("Returning integration\n");
    return (1);
}
int couplingOutput(int *doy) {
    // Get an instance of Simulator
    Simulator *s = Simulator::getInstance();
    // Request disease outputs to be written in files
    //s->output();
    
    return (1);
}