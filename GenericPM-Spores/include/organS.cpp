#include "organS.h"
#include "simulatorS.h"
#include "cloudoS.h"
#include "cloudpS.h"
#include "cloudfS.h"
#include "lesioncohortS.h"
#include "utilitiesS.h"
#include "basicS.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

int OrganS::firstOutputCallS = 0;

void OrganS::integrationS() {
    if (!isAlive()) return;
    
    CloudOS *cloudo;
    LesionCohortS *lc;

    float cloudDensity = 0;
    double cloudOValue = 0, cloudPValue = 0, cloudFvalue = 0;
    dailyDiseaseArea = diseaseArea;
    dailyVisibleDiseaseArea = visibleDiseaseArea;
    dailySenescenceArea = senescenceArea;
    dailyLatentDiseaseArea = latentDiseaseArea;
    dailyInfectionDiseaseArea = infectionDiseaseArea;
    dailyNecroticDiseaseArea = necroticDiseaseArea;
    dailyTotalLesions = totalLesions;
    dailyVisibleLesions = visibleLesions;
    diseaseArea = 0;
    visibleLesions = 0;
    visibleDiseaseArea = invisibleDiseaseArea = 0;
    latentDiseaseArea = infectionDiseaseArea = necroticDiseaseArea = 0;    newLesionsS = 0;

    if (suceptible) {
        for (unsigned int i = 0; i < cloudsO.size(); i++) {
            cloudo = &cloudsO[i];
            physiologicalLife += UtilitiesS::trapezoidalFunctionS(BasicS::getWeather()->getTMean(), cloudo->getDisease()->getCardinalTempPhysiologicalLife());

            cloudOValue = cloudAmountS();
            cloudPValue = cloudo->getCloudP()->getValueS();
            cloudFvalue = cloudo->getCloudP()->getCloudF()->getValueS();

            newLesionsFromOrgan = cloudo->getDisease()->newLesionsS(cloudOValue,healthAreaProportion);
            newLesionsFromPlant = cloudo->getDisease()->newLesionsS(cloudPValue,healthAreaProportion)*getProportionFromTotalArea();
            newLesionsFromField = cloudo->getDisease()->newLesionsS(cloudFvalue,healthAreaProportion)*getProportionFromTotalArea();
            //std::cout << "cloudF: " << cloudFvalue << std::endl;
            //std::cout << "cloudP: " << cloudPValue << std::endl;
            //std::cout << "cloudO: " << cloudOValue << std::endl;
            ////std::cout << "healthAreaProportion " << healthAreaProportion << " getProportionFromTotalArea " << getProportionFromTotalArea() <<
            //" newLesionsFromField " << newLesionsFromField << std::endl;
            if ((newLesionsFromOrgan+newLesionsFromPlant+newLesionsFromField) > 0 && physiologicalLife >= 5) {
                newLesionsS = newLesionsFromOrgan+newLesionsFromPlant+newLesionsFromField;
                lesionCohorts.emplace_back(newLesionsS, cloudo);
                totalLesions += newLesionsS;

                // Add Spores that will be removed because were used to infect the tissue
                cloudo->addSporesToBeRemoved(newLesionsFromOrgan);
                cloudo->getCloudP()->addSporesToBeRemoved(newLesionsFromPlant);
                cloudo->getCloudP()->getCloudF()->addSporesToBeRemoved(newLesionsFromField);
            } else {
                newLesionsFromOrgan = newLesionsFromPlant = newLesionsFromField = 0;
            }
        }
    }

    for (unsigned int i = 0; i < lesionCohorts.size(); i++) {
        lc = &lesionCohorts[i];
        lc->integrationS();
        diseaseArea += lc->getTotalArea();
        latentDiseaseArea += lc->getLatentArea();
        infectionDiseaseArea += lc->getInfectionArea();
        necroticDiseaseArea += lc->getNecroticArea();
        visibleLesions += lc->getVisibleLesions();
        visibleDiseaseArea += lc->getVisibleArea();
        invisibleDiseaseArea += lc->getInvisibleArea();
    }

    cloudIntegrationS();

    cloudOValue = cloudAmountS();
    cloudPValue = cloudo->getCloudP()->getValueS();
    cloudFvalue = cloudo->getCloudP()->getCloudF()->getValueS();

    dailyDiseaseArea = fmax(0,diseaseArea - dailyDiseaseArea);
    dailyVisibleDiseaseArea = fmax(0,visibleDiseaseArea - dailyVisibleDiseaseArea);
    dailySenescenceArea = fmax(0,senescenceArea - dailySenescenceArea);
    dailyLatentDiseaseArea = fmax(0,latentDiseaseArea - dailyLatentDiseaseArea);
    dailyInfectionDiseaseArea = fmax(0,infectionDiseaseArea - dailyInfectionDiseaseArea);
    dailyNecroticDiseaseArea = fmax(0,necroticDiseaseArea - dailyNecroticDiseaseArea);
    dailyTotalLesions = fmax(0,totalLesions - dailyTotalLesions);
    dailyVisibleLesions = fmax(0,visibleLesions - dailyVisibleLesions);
    //    printf(" Daily: DiseaseArea (%f),"
    //            "VisibleDiseaseArea (%f),"
    //            "SenescenceArea (%f),"
    //            "LatentDiseaseArea (%f),"
    //            "InfectionDiseaseArea (%f),"
    //            "NecroticDiseaseArea (%f)\n "
    //            "TotalLesions (%i),"
    //            "VisibleLesions (%i)\n", dailyDiseaseArea, dailyVisibleDiseaseArea,
    //            dailySenescenceArea,dailyLatentDiseaseArea,dailyInfectionDiseaseArea,
    //            dailyNecroticDiseaseArea,dailyTotalLesions,dailyVisibleLesions);
    std::ostringstream convert;
    convert << organNumber << "," 
            << BasicS::getWeather()->getYearDoy() << "," << totalArea << "," << senescenceArea << "," 
            << UtilitiesS::formatDoubleS(diseaseArea, 4) << ","
            << UtilitiesS::formatDoubleS(visibleDiseaseArea, 4) << "," 
            << UtilitiesS::formatDoubleS(invisibleDiseaseArea, 4) << ","
            << UtilitiesS::formatDoubleS(totalLesions / totalArea, 4) << ","
            << UtilitiesS::formatDoubleS(physiologicalLife, 4) << "," << dailyTotalLesions << "," 
            << totalLesions << "," 
            << UtilitiesS::formatDoubleS(cloudOValue) << "," << UtilitiesS::formatDoubleS(cloudPValue) << ","
            << UtilitiesS::formatDoubleS(cloudFvalue) << "," << healthAreaProportion << "," 
            << getProportionFromTotalArea() << "," << BasicS::getWeather()->getWetDur() << ","
            << newLesionsFromOrgan << "," << newLesionsFromPlant << "," << newLesionsFromField;
    BasicS::output.push_back(convert.str());
}
void OrganS::cloudIntegrationS() {
    for (unsigned int i = 0; i < cloudsO.size(); i++) {
        (&cloudsO[i])->integrationS();
    }
}
float OrganS::cloudAmountS() {
    float cloudOValue=0;
    for (unsigned int i = 0; i < cloudsO.size(); i++) {
        cloudOValue += (&cloudsO[i])->getValueS();
    }
    return cloudOValue;
}
void OrganS::outputS() {
    std::ostringstream convert;
    //convert << "Cpp_Organ_" << getOrganNumber() << ".txt";
    BasicS::getOutput("Cpp_Organ.txt",this->firstOutputCallS);
    this->firstOutputCallS++;
    // Speedup the model removing outputs
    //std::cout << "\nOrgan " << getID() << ":";
    //for(unsigned int i=0; i<BasicS::output.size(); i++)
    //{
    //    std::cout << BasicS::output[i] << std::endl;
    //}

    CloudOS *cloudo;
    for (unsigned int i = 0; i < cloudsO.size(); i++) {
        cloudo = &cloudsO[i];
        cloudo->outputS();
    }
    
//    LesionCohortS *lc;
//    for (unsigned int i = 0; i < lesionCohorts.size(); i++) {
//        lc = &lesionCohorts[i];
//        lc->outputS();
//    }

}

void OrganS::rateS() {
    // Calculate the ratio due senescence based on previews day
    float actualDisease=0, ratioSenescence = this->senescenceArea / this->totalArea;
    // Update the senescence area for the current day
    this->senescenceArea = SimulatorSpore::getInstance()->getCropInterface()->getSenescenceOrganArea(organNumber);
    // Recalculate the ratio due senescence and take the difference from previews ratio
    ratioSenescence = (this->senescenceArea / this->totalArea) - ratioSenescence;
    // Update the total organ area (current day)
    this->totalArea = SimulatorSpore::getInstance()->getCropInterface()->getOrganArea(organNumber);
    
    if (!suceptible && this->totalArea > 0) {
        suceptible = true;
    }
    if (!isAlive())
        return;
    CloudOS *cloudo;
    LesionCohortS *lc;
    //printf("OrganNumber: %i OrganTotalArea: %f healthAreaProportion: %f \n",this->organNumber, this->totalArea,healthAreaProportion);
    for (unsigned int i = 0; i < lesionCohorts.size(); i++) {
        lc = &lesionCohorts[i];
        // Must affect the disease area related with senescent area
        if (ratioSenescence > 0) {
            //printf("ID: %i Age: %i VisibleArea: %f InvisibleArea: %f RatioDueDefoliation: %f\n", lc->getID(), lc->getAge(), lc->getVisibleArea(),lc->getInvisibleArea(),CropInterfaceS::getInstance()->getRatioDueDefoliation(1));
            lc->setVisibleArea(lc->getVisibleArea() * (1-ratioSenescence));
            lc->setInvisibleArea(lc->getInvisibleArea() * (1-ratioSenescence));
            lc->setTotalArea(lc->getVisibleArea()+lc->getInvisibleArea());
            //printf("ID: %i Age: %i VisibleArea: %f InvisibleArea: %f \n", lc->getID(), lc->getAge(), lc->getVisibleArea(),lc->getInvisibleArea());
        }
        actualDisease = actualDisease + lc->getTotalArea();
    }
    // Updating the disease amount on organ
    setDiseaseArea(actualDisease);
    // Calculating the health area proportion
    healthAreaProportion = UtilitiesS::getHealthAreaProportionS(getDiseaseArea(), getTotalArea(), getSenescenceArea());

    for (unsigned int i = 0; i < lesionCohorts.size(); i++) {
        lc = &lesionCohorts[i];
        lc->setOrganHealthAreaProportion(healthAreaProportion);
        lc->rateS();
    }
}
