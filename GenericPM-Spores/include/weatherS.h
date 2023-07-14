#ifndef WEATHERS_H
#define WEATHERS_H

class Weather {
protected:
    Weather();
    static Weather* instance;
    int yearDoy = 0;
    int year = 0;
    int doy = 0;
    float tMin = 0;
    float tMax = 0;
    float tMean = 0;
    float sRad = 0;
    float rain = 0;
    float par = 0;
    float wetDur = 0;
    float hWetDur = 0;
    float rh = 0;
    float rhMax = 0;
    float rhMin = 0;
    float rhMean = 0;
    int hRH90 = 0;

public:
    static Weather* getInstance();

    float getWetDur() {
        if (wetDur == 0)
            return hRH90;
        else
            return wetDur;
    }
    bool update();

    void setRain(float rain) {
        this->rain = rain;
    }

    float getRain() const {
        return rain;
    }

    void setTMean(float tMean) {
        this->tMean = tMean;
    }

    float getTMean() const {
        return tMean;
    }

    void setDoy(int doy) {
        this->doy = doy;
    }

    int getDoy() const {
        return doy;
    }

    void setYearDoy(int yearDoy) {
        this->yearDoy = yearDoy;
    }

    int getYearDoy() const {
        return yearDoy;
    }

    float getRh() const {
        return rh;
    }


};
#endif // WEATHER_H
