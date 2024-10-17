!=======================================================================
!  Aloha_GROSUB, Subroutine
!
!  Maize growth routine
!-----------------------------------------------------------------------
!  Revision history
!  02/07/1993 PWW Header revision and minor changes   
!  02/07/1993 PWW Switch block added, etc
!  10/17/2017 CHP Adpated for CSM v4.7
!  09/05/2020 JVJ Stages changes for inclusion in Overview   
!-----------------------------------------------------------------------
!  INPUT  : NOUTDO,ISWNIT
!
!  LOCAL  : NSINK,NPOOL1,NPOOL2,NPOOL,NSDR,I,ICOLD,PCARB,PRFT,PC,
!           GRF,GROEAR,RGFILL,TTMP,GROGRN,SFAC,TFAC,RMNC,XNF,TNLAB,
!           RNLAB,RNOUT,SLFW,SLFN,SLFC,SLFT,PLAS,TI
!
!  OUTPUT :
!-----------------------------------------------------------------------
!  Called : PINE
!
!  Calls  : NFACTO NUPTAK
!-----------------------------------------------------------------------
!                         DEFINITIONS
!
!  GRF    :
!  GROGRN : Daily growth of the grain - g
!  I      : Loop counter
!  ICOLD  :
!  NOUTDO : File handle
!  NPOOL  : Total plant N available for translocation to grain (g/plant)
!  NPOOL1 : Tops N available for translocation to grain (g/plant)
!  NPOOL2 : Root N available for translocation to grain (g/plant)
!  NSDR   : Plant N supply/demand ratio used to modify grain N content
!  NSINK  : Demand for N associated with grain filling (g/plant/day)
!  PAR    : Daily photosynthetically active radiation, calculated as half
!           the solar radiation - MJ/square metre
!  PC     :
!  PCARB  : Daily amount of carbon fixed - g
!  PLAS   : The rate of senescence of leaf area on one plant - sq. cm/day
!  PRFT   : Photosynthetic reduction factor for low and high temperatures
!  RGFILL : Rate of grain fill - mg/day
!  RMNC   : Root minimum nitrogen concentration (g N/g root dry weight)
!  TI     : Fraction of a phyllochron interval which occurred as a fraction
!           of today's daily thermal time
!  TNLAB  :
!  TTMP   :
!=======================================================================

      SUBROUTINE Aloha_GROSUB (CONTROL, ISWITCH, 
     &    DTT, ISTAGE, NH4, NO3, SOILPROP, SW, SWFAC,!Input
     &    SUMDTT, TBASE, TURFAC, WEATHER, XSTAGE, DAP1, DAP3, DAP5, DAP7, DAP9, DAP13, DAP15, DAP17, DAP19, DAP21,             !Input
     &    AGEFAC, BASLFWT, BIOMAS, CRWNWT, EYEWT, FBIOM, MAXLAI, WEATHERFact, TMAXGROF, SRADGROF,      !Output
     &    FLRWT, FRTWT, FRUITS, GPP, GPSM, GRAINN, GRORT, SUMSRADGRO, SUMSRAD, SRADGRO, PARGRO, SUMPARGRO, SUMPAR,    !Output
     &    LAI, LFWT, LN, NSTRES, RLV, ROOTN, RTWT, SUMDTTGRO, SUMTMAXGRO, !Output
     &    SENESCE, SKWT, STMWT, STOVN, STOVWT,  TEMPM, SUMTMAX, TMAXGRO, GDDFR, BIOMAS4, LAI4, LN2, LN3, LN4, HIFact, !Output
     &    UNH4, UNO3, WTNUP, WTINITIAL, XGNP, YIELD, BIOMAS1, LAI1, BIOMAS13, LAI13, BIOMAS2, LAI2, BIOMAS3, LAI3)!Output 

      USE Aloha_mod
      USE Interface_SenLig_Ceres
      IMPLICIT  NONE
      SAVE

      INTEGER   ICOLD
      REAL      PCARB,PRFT,PC,TI,GRF,RGFILL,
     &          SLFW,SLFN,SLFC,SLFT,PLAS
      REAL      TABEX,PCO2,Y1, XLAT

      CHARACTER ISWNIT*1
      REAL    GROGRN,SFAC,TFAC,RMNC,XNF,TNLAB,RNLAB
      REAL    NSINK,NPOOL1,NPOOL2,NPOOL,NSDR,RNOUT
      INTEGER ICSDUR
      REAL    SEEDNI, ROOTN, STOVN, GRAINN, SEEDN, XANC
      REAL    APTNUP, RANC, GNP, NFAC, RCNP 
      REAL    TANC, VANC, VMNC, TMNC

      CHARACTER*6, PARAMETER :: ERRKEY ='GROSUB'
      CHARACTER*78 MSG(2)
      INTEGER I, ISTAGE, ISTAGE_old, IDURP, TIMDIF, YRDOY
      INTEGER STGDOY(20), YRPLT
      INTEGER DYNAMIC
      REAL    PLA, LAI, BIOMAS, LFWT, BASLFWT, STMWT, STOVWT, WTINITIAL, BIOMAS1, BIOMAS13, BIOMAS2, BIOMAS3, BIOMAS4, SEEDQLY 
      REAL    RLAE13, RLDW13, RBWTDW13, RSTMWT13, LAI13, PLA12, PLA13, LFWT12, LFWT13, BASLFWT12, BASLFWT13, STMWT12, STMWT13, ADJGDDF13, ADJTMAX13, ADJSRAD13  !CASE 13
      REAL    RLAE1, RLDW1, RBWTDW1, RSTMWT1, LAI1, PLA1, LFWT1, BASLFWT1, STMWT1, GRORTI, RGRORT1, GRORT1, TI1, ADJGDDF1, ADJTMAX1, ADJSRAD1  !CASE 1
      REAL    RLAE2, RLDW2, RBWTDW2, RSTMWT2, LAI2, LN2, PLA2, LFWT2, BASLFWT2, STMWT2, RGRORT2, GRORT2, RTI2, TI2 !CASE 2
      REAL    RLAE3, RLDW3, RBWTDW3, RSTMWT3, LAI3, LN3, PLA3, LFWT3, BASLFWT3, STMWT3, RGRORT3, GRORT3, RTI3, TI3 !CASE 3
      REAL    RLAE4, RLDW4, RBWTDW4, RSTMWT4, LAI4, LN4, PLA4, LFWT4, BASLFWT4, STMWT4, RGRORT4, GRORT4, RTI4, HIFact, YIELDFact !CASE 4
      REAL    RLAE6, RLDW6, RBWTDW6, RSTMWT6, LAI6, LN6, PLA6, LFWT6, BASLFWT6, STMWT6, RGRORT6, GRORT6 !CASE 6
      REAL    RLAE7, RLDW7, RBWTDW7, RSTMWT7, LAI7, LN7, PLA7, LFWT7, BASLFWT7, STMWT7, RGRORT7, GRORT7 !CASE 7
      REAL    RLAE8, RLDW8, RBWTDW8, RSTMWT8, LAI8, LN8, PLA8, LFWT8, BASLFWT8, STMWT8, RGRORT8, GRORT8 !CASE 8
      REAL    RLAE9, RLDW9, RBWTDW9, RSTMWT9, LAI9, LN9, PLA9, LFWT9, BASLFWT9, STMWT9, RGRORT9, GRORT9 !CASE 9-10
      REAL    DAP1, DAP2, DAP3, DAP4, DAP5, DAP6, DAP7, DAP8, DAP9, DAP10, DAP11, DAP12, DAP13, DAP14, DAP15, DAP16, DAP17, DAP18, DAP19, DAP20, DAP21
      REAL    PLAG, RTWT, FLRWT, GROSTM, SENLA, SLAN, GRORT, GDDFR, TMAXGRO, PLACASE12, PLTadj, WEATHERFact, TMAXGROF, SRADGROF
      REAL    GDDFRS6, GDDFRS7, GDDFRS8, GDDFRS9, TMAXGROS6, TMAXGROS7, TMAXGROS8, TMAXGROS9, SRADGROS6, SRADGROS7, SRADGROS8, SRADGROS9
      REAL    PARGROS6, PARGROS7, PARGROS8, PARGROS9
      REAL    GROBSL, GROLF, CUMPH, LN, CUMDEP, SUMP, PLAMX, GROFLR
      REAL    GROCRWN, GROFRT, FRTWT, CRWNWT, SKWT, GROSK, PTF, EYEWT
      REAL    SWMAX, SWMIN, NDEF3, NSTRES, AGEFAC, LIFAC, SRADGRO, PARGRO
      REAL    PAR, CC, TRF2, CARBO, SWFAC, TEMPM  !,TRNU, 
      REAL    DTT, TURFAC, XN, CMF, TOTPLTWT, SUMDTT, GPP, SUMTMAX, SUMDTTGRO, SUMTMAXGRO
      REAL    PDWI, PGRORT, DM, FBIOM, MAXLAI, PHOTOSYNEYE, FRUITS, SUMSRADGRO, SUMSRAD, SUMPARGRO, SUMPAR
      REAL    YIELD, GPSM, XSTAGE, YIELD1, TEMPAd  !, FDMC

      REAL    CO2, SRAD, TMIN, TMAX
      REAL    PLTPOP, SDWTPL, PLANTSIZE
      REAL    G2, G3, TC, P1, P2, P3, P4, P8, PHINT, TBASE      !G2, G3, P4, PHINT, TBASE G2, G3, P7, PHINT, TBASE
      INTEGER PMTYPE, NDOF, NFORCING
      REAL    GRNWT, SDWTAH, SDWTAM, WTNUP, BWAH
      REAL    WTNLF, WTNST, WTNSH, WTNRT, WTNLO
      REAL    NDEF4, ANFAC, ATANC, TCNP, XGNP, GNUP, TOTNUP
      REAL    CUMDTT, CANNAA, CANWAA
      REAL    PLIGLF, PLIGRT

      REAL, DIMENSION(10) :: CO2X, CO2Y
      REAL, DIMENSION(NL) :: RLV, NO3, NH4, SW, UNH4, UNO3

      TYPE (ControlType) CONTROL
      TYPE (SwitchType)  ISWITCH
      TYPE (WeatherType) WEATHER
      TYPE (SoilType) SOILPROP
      TYPE (ResidueType) SENESCE

      DYNAMIC = CONTROL % DYNAMIC
      CO2  = WEATHER % CO2
      SRAD = WEATHER % SRAD
      TMIN = WEATHER % TMIN
      TMAX = WEATHER % TMAX
      XLAT = WEATHER % XLAT

!=======================================================================
      SELECT CASE (DYNAMIC)
!=======================================================================
      CASE (RUNINIT)
!=======================================================================
      ISWNIT     = ISWITCH % ISWNIT
      PLA        = 0.0
      PLACASE12  = 0.0
      LAI        = 0.0
      BIOMAS     = 0.0
      LFWT       = 0.0
      BASLFWT    = 0.0
      STMWT      = 0.0
      STOVWT     = 0.0
      SWFAC      = 1.0
      TURFAC     = 1.0
      LN     = 0.0
      SUMDTTGRO = 0.0
      SUMTMAXGRO = 0.0
      SUMTMAX = 0.0
      SUMSRADGRO = 0.0
      SUMSRAD = 0.0
      SUMPARGRO = 0.0
      SUMPAR = 0.0
      FLRWT  = 0.0
      FRTWT  = 0.0
      CRWNWT = 0.0
      SKWT   = 0.0
      GROSK  = 0.0
      YIELD  = 0.0
      SENLA  = 0.0
      SLAN   = 0.0
      CARBO  = 0.0
      GRNWT  = 0.0  !Not ever given a value, but used to compute SDWT
      RTWT   = 0.0
      SDWTAH = 0.0  !Only used for output in OPHarv
      SDWTAM = 0.0
      BWAH   = 0.0
      WTNLF  = 0.0
      WTNST  = 0.0
      WTNSH  = 0.0
      WTNRT  = 0.0
      WTNLO  = 0.0
      GPSM   = 0.0
      GPP    = 0.0
      PTF    = 0.0
      FRUITS = 0.0
      PLTadj = 0.0
      

      DO I = 1, NL
         RLV(I) = 0.0
      END DO

      BIOMAS = 0.0
      LAI    = 0.0
      XN     = 0.0
      SWFAC  = 1.0
      TURFAC = 1.0
      NDEF4  = 1.0
      ANFAC  = 0.0
      ATANC  = 0.0
      VANC   = 0.0
      VMNC   = 0.0
      SEEDNI = 0.0
      GRAINN = 0.0
      GNP    = 0.0
      XGNP   = 0.0
      APTNUP = 0.0
      GNUP   = 0.0
      TOTNUP = 0.0
      CUMDTT = 0.0
      SUMDTT = 0.0
      DTT    = 0.0
      CANNAA = 0.05
      CANWAA = 0.0
      SUMDTTGRO = 0.0
      SUMTMAXGRO = 0.0
      SUMTMAX = 0.0
      SUMSRADGRO = 0.0
      SUMSRAD = 0.0
      SUMPARGRO = 0.0
      SUMPAR = 0.0

      PLAG    = 0.0   ! PLAG (cm^2) is daily green leaf area growth
      GROSTM  = 0.0   ! GROSTM (g/plant/day) is daily stem growth
      GRORT   = 0.0   ! GRORT (g/plant/day) is daily root growth
      GROBSL  = 0.0   ! GROBSL (g/plant/day) is daily basal leaf growth
      GROLF   = 0.0   ! GROLF (g/plant/day) is daily green leaf growth
      GROFLR  = 0.0
      GROCRWN = 0.0
      GROFRT  = 0.0

      CALL Aloha_NFACTO (DYNAMIC, 
     &    ISTAGE, TANC, XSTAGE,                           !Input
     &    AGEFAC, NDEF3, NFAC, NSTRES, RCNP, TCNP, TMNC)  !Output

      CALL Aloha_NUPTAK(CONTROL, ISWITCH, 
     &    ISTAGE, NO3, NH4, PDWI, PGRORT, PLIGRT,         !Input
     &    PLTPOP, PTF, RANC, RCNP, RLV, RTWT, SOILPROP,   !Input
     &    STOVWT, SW, TCNP, XSTAGE,                       !Input
     &    ROOTN, SENESCE, STOVN, TANC, UNH4, UNO3, WTNUP) !Output

!=======================================================================
      CASE (SEASINIT)
!=======================================================================
      SDWTPL    = PLANTING % SDWTPL
      PLTPOP    = PLANTING % PLTPOP
      PMTYPE    = PLANTING % PMTYPE
      NFORCING  = PLANTING % NFORCING
      PLANTSIZE = PLANTING % PLANTSIZE
      NDOF     = Planting % NDOF

      CO2X = SPECIES % CO2X
      CO2Y = SPECIES % CO2Y
      CC   = SPECIES % CONV      
      CMF  = Species % CMFC
      LIFAC= Species % LIFAC

      G2  = CULTIVAR % G2
      G3  = CULTIVAR % G3
      P8  = CULTIVAR % P8                   
      PHINT = CULTIVAR % PHINT
      TC  = CULTIVAR % TC
      P1  = CULTIVAR % P1
      P2  = CULTIVAR % P2
      P3  = CULTIVAR % P3
      P4  = CULTIVAR % P4

      PLA        = 0.0
      LAI        = 0.0
      BIOMAS     = 0.0
      LFWT       = 0.0
      BASLFWT    = 0.0
      STMWT      = 0.0
      STOVWT     = 0.0
      LN     = 0.0

!     Calculate initial SEED N
      SEEDNI = (ROOTN+STOVN+GRAINN+SEEDN)*PLTPOP

      ISTAGE_OLD = 0

!     Initialize senescence variables
      CALL SenLig_Ceres(PLIGLF=PLIGLF, PLIGRT=PLIGRT)

      CALL Aloha_NFACTO (DYNAMIC, 
     &    ISTAGE, TANC, XSTAGE,                           !Input
     &    AGEFAC, NDEF3, NFAC, NSTRES, RCNP, TCNP, TMNC)  !Output

      CALL Aloha_NUPTAK(CONTROL, ISWITCH, 
     &    ISTAGE, NO3, NH4, PDWI, PGRORT, PLIGRT,         !Input
     &    PLTPOP, PTF, RANC, RCNP, RLV, RTWT, SOILPROP,   !Input
     &    STOVWT, SW, TCNP, XSTAGE,                       !Input
     &    ROOTN, SENESCE, STOVN, TANC, UNH4, UNO3, WTNUP) !Output

!=======================================================================
      CASE (RATE)
!=======================================================================
      !TEMPM = (WEATHER % TMAX + WEATHER % TMIN) / 2.

      IF (TMIN .GT. TBASE .AND. TMAX .LT. 38.0) THEN
             IF (XLAT .LT. 21.0 .and. XLAT .GT. -21.0) THEN
                TEMPM = 0.6*TMIN+0.4*TMAX
              ELSE
                TEMPM = (TMAX+TMIN)/2
             ENDIF
             ENDIF

      IF (ISWNIT .NE. 'N') THEN
!       Top actual N concentration (g N/g Dry weight)
        XANC   = TANC*100.0               
        APTNUP = STOVN*10.0*PLTPOP

        IF (ISTAGE .LT. 10) THEN                      
          CALL Aloha_NFACTO (DYNAMIC, 
     &      ISTAGE, TANC, XSTAGE,                           !Input
     &      AGEFAC, NDEF3, NFAC, NSTRES, RCNP, TCNP, TMNC)  !Output
        ENDIF
      ENDIF

      IF (ISTAGE .GT. 10) RETURN                     

!-----------------------------------------------------------------

       IF (TMAX .GE. 33) THEN
       TEMPAd = TMAX
       ELSE
       TEMPAd = (TMAX - 25)
       ENDIF

       Yieldfact =  EXP(-TEMPAd*(TEMPAd/20000))

      !HIFact = Yieldfact !0.935 + ((0.082*EXP(0.1*TEMPM))/(20))  !factor fruit
      YIELD1 = Yieldfact * 0.0520 !0.065 + ((0.082*EXP(0.1*TEMPM))/(20))  !Factor crown
      PAR   = 0.5*SRAD
      Y1    = EXP(-LIFAC*LAI)                       ! Beer's law
      
      PLTadj = -0.0172 * PLTPOP**2 + 0.2162 * PLTPOP + 0.9308    ! File DSSAT ecuations y = -0.0172x2 + 0.2162x + 0.9308 R² = 0.9005 
       
      PCARB = CC*PAR/(PLTPOP/PLTadj)*(1.0-Y1)                ! on per plant basis
      !
      ! Calculate Photosynthetic Response to CO2
      !
      PCO2  = TABEX (CO2Y,CO2X,CO2,10)
      PCARB = PCARB*PCO2

      !TEMPM = 0.6*TMIN + 0.4*TMAX
      IF (TMIN .GT. TBASE .AND. TMAX .LT. 38.0) THEN
             IF (XLAT .LT. 21.0 .and. XLAT .GT. -21.0) THEN
                TEMPM = 0.6*TMIN+0.4*TMAX
              ELSE
                TEMPM = (TMAX+TMIN)/2
             ENDIF
             ENDIF
      SELECT CASE (ISTAGE)
        CASE (1,2,3,4,11,12,13)                      
          IF (TEMPM .LE. 20.0) THEN
             PRFT = 1.0-0.015*(TEMPM-25.0)**2    !  PRFT   : Photosynthetic reduction factor for low and high temperatures
           ELSEIF (TEMPM .LE. 25.0) THEN
             PRFT = 1.0-0.001*(TEMPM-25.0)**2    !  PRFT   : Photosynthetic reduction factor for low and high temperatures
           ELSEIF (TEMPM .LT. 29.0) THEN
             PRFT = 1.0-0.056*(TEMPM-25.0)**2
           ELSE
             PRFT = 0.1                                   
          ENDIF
        CASE (5,6,7,8,9,10)                                   
          IF (TEMPM .LE. 20.0) THEN
             PRFT = 1.0-0.190*(TEMPM-20.0)**2    !  PRFT   : Photosynthetic reduction factor for low and high temperatures
           ELSEIF (TEMPM .LT. 25.0) THEN
             PRFT = 1.0-0.005*(TEMPM-25.0)**2
           ELSEIF (TEMPM .LT. 29.0) THEN
             PRFT = 1.0-0.0504*(TEMPM-26.)**2
           ELSE
             PRFT = 0.1                                   
          ENDIF
          
          
          
          
          PRFT = AMAX1 (PRFT,0.0)
      END SELECT
!-----------------------------------------------------------------
      !
      ! Temperature factor
      !
      IF (TEMPM .LT. 20.0) THEN                               ! 
         TRF2 = 0.35                                          ! 
       ELSEIF (TEMPM .GE. 20.0 .AND. TEMPM .LT. 26.0) THEN    ! 
         TRF2 = 1.13314845306683*EXP(-0.005*TEMPM)   !1-0.005*(TEMPM-25) !
       ELSE
         TRF2 = 0.98                                          ! 
      ENDIF

      IF (ISTAGE .GE. 4 .AND. ISTAGE .LT. 13) THEN                
         CARBO = TRF2*PCARB*AMIN1(PRFT,TRF2*SWFAC,NSTRES)      
       ELSE
         CARBO = TRF2*PCARB*AMIN1(PRFT,SWFAC,NSTRES)
      ENDIF
      DTT = AMAX1 (DTT,0.0)
      
!----------------------------------------------------------------- Leaf Number modification
      IF (ISTAGE .LE. 4 .OR. ISTAGE .EQ. 13) THEN                                   
!                                                                 
!        Calculate leaf emergence                                 
!         
        ! PC Used to compute fraction of phyllochron interval occurring today
         IF (XN .LE. 13.0) THEN                                               
            PC = 1.0 - (PHINT*0.0001)*CUMPH                   
            
            ELSEIF (XN .GT. 13 .AND. XN .LE. 26.0) THEN   
                PC = 0.9 - (PHINT*0.0001)*CUMPH      !0.4              
          
                ELSEIF (XN .GT. 26.0 .AND. XN .LE. 39.0) THEN                            
                     PC = 0.85 - (PHINT*0.0001)*CUMPH     !0.5         
                     
                     ELSEIF (TMIN .LE. TBASE .OR. TMAX .GE. 38.0) THEN 
                        IF (TMAX .LT. TBASE) THEN     
                          PC = 1.85 - (0.25/20)*CUMPH
               ENDIF
               ENDIF
!         
!        TI is the fraction of leaf emerged for a day.  It is calculated from
!        the following equations
!        
!        Correcting water stress effect and effect due to shading.
!                                                            
         IF (ISTAGE .LE. 4) THEN                                
            IF (TEMPM .GE. TBASE) THEN
               IF ((LN*PLTPOP) .LE. (13*PLTPOP)) THEN        
                    TI1       = -0.07260309 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 + 0.29041612 * (LOG(GDDFR/(TMAXGRO/SRADGRO))) - 0.15937426  !y = -0.07260309x2 + 0.29041612x - 0.15937426 R² = 0.97140069
                    TI     = TI1/PC
                 
                   ENDIF

               IF ((LN*PLTPOP) .LE. (26*PLTPOP)) THEN
                    TI2       = -0.07260309 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 + 0.29041612 * (LOG(GDDFR/(TMAXGRO/SRADGRO))) - 0.15937426  !y = -0.07260309x2 + 0.29041612x - 0.15937426 R² = 0.97140069
                    TI     = TI2/PC
                   
                   ENDIF

                 IF ((LN*PLTPOP) .LE. (39*PLTPOP)) THEN                   
                    TI3       = -0.07260309 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 + 0.29041612 * (LOG(GDDFR/(TMAXGRO/SRADGRO))) - 0.15937426  !y = -0.07260309x2 + 0.29041612x - 0.15937426 R² = 0.97140069
                    TI     = TI3/PC
               
                
                 ENDIF
                  ENDIF     
                  
          
          IF (ISTAGE .EQ. 5 .OR. TMIN .LE. TBASE) THEN     !This stage is forcing therefore no more leaves are produced.
           TI = 0.0    
           ENDIF
         ENDIF

!        CUMPH is number of expanded leaves. It is updated daily

         CUMPH = CUMPH + TI
                                           
!        XN is leaf number of the oldest expanding leaf
         XN    = CUMPH +  1 
         LN    = XN         ! LN is leaf number
      ENDIF  ! CIERRA IF (ISTAGE .LE. 4)     

!-----------------------------------------------------------------
!  ISTAGE Definition
!    11 Start simulation to planting
!    12 Planting to Root Initiation
!    13 Root Initiation to First New Leaf
!     1 First new leaf emergence to foliar cycle 1
!     2,3,4 Foliar cycle 1 to foliar cycle 2,3 and forcing 
!     5 Forcing to Open Heart
!     6 Open Heart to EarlyAnthesis
!     7 Early Anthesis to Last Anthesis
!     8 Last Anthesis to Physiological maturity
!     9 Physiology to Harvest
!    10 Harvest
!-----------------------------------------------------------------
      SELECT CASE (ISTAGE)
!-----------------------------------------------------------------

       CASE (12)  
        !
        !      ! Planting to Root Initiation
        !
          LFWT        = ((LFWT12/TC)*DTT)     + LFWT                     
          BASLFWT     = ((BASLFWT12/TC)*DTT)  + BASLFWT               
          STMWT       = ((STMWT12/TC)*DTT)    + STMWT  
        

        !
        ! Check the balance of supply and demand
        !
        GRORT = CARBO - GROLF - GROBSL - GROSTM
        IF (GRORT .LT. 0.15*CARBO) THEN
           IF (GROLF .GT. 0.0 .OR. GROBSL .GT. 0.0 .OR.
     &         GROSTM .GT. 0.0) THEN
              GRF   = CARBO*0.9/(GROLF+GROBSL+GROSTM)
              GRORT = CARBO*0.1
            ELSE
              GRF = 1.0
           ENDIF
        ENDIF
   !   ENDIF

                                        
!-----------------------------------------------------------------         

       CASE (13)  
        !
        !      !Root Initiation to First New leaf 
        !
          LFWT        = ((LFWT13/P1)*DTT)     + LFWT                
          BASLFWT     = ((BASLFWT13/P1)*DTT)  + BASLFWT           
          STMWT       = ((STMWT13/P1)*DTT)    + STMWT  
        
        !
        ! Check the balance of supply and demand
        !
        GRORT = CARBO - GROLF - GROBSL - GROSTM
        IF (GRORT .LT. 0.15*CARBO) THEN
           IF (GROLF .GT. 0.0 .OR. GROBSL .GT. 0.0 .OR.
     &         GROSTM .GT. 0.0) THEN
              GRF   = CARBO*0.9/(GROLF+GROBSL+GROSTM)
              GRORT = CARBO*0.1
            ELSE
              GRF = 1.0
           ENDIF
        ENDIF
         
!-----------------------------------------------------------------

      CASE (1) 
        !
        !      !First New leaf  to Foliar Cycle 1
        !
          
        
!      IF ((LN) .LE. (13)) THEN

          LFWT        = ((LFWT1/P2)*DTT)     + LFWT            
          BASLFWT     = ((BASLFWT1/P2)*DTT)  + BASLFWT                  
          STMWT       = ((STMWT1/P2)*DTT)    + STMWT
          GRORT       = ((GRORT1/P2)*DTT)    + GRORT              
          MAXLAI      = ((LAI1/P2)*DTT)    + MAXLAI                  
        
!-----------------------------------------------------------------
      CASE (2)

        !
        !      Leaf cycle 1 to Leaf cycle 2
        !




          LFWT        = ((LFWT2/P3)*DTT)     + LFWT                      
          BASLFWT     = ((BASLFWT2/P3)*DTT)  + BASLFWT                   
          STMWT       = ((STMWT2/P3)*DTT)    + STMWT
          GRORT       = ((GRORT2/P3)*DTT)    + GRORT  
          MAXLAI      = ((LAI2/P3)*DTT)      + MAXLAI  
            
         
        !
        ! Check the balance of supply and demand
        !
        GRORT = CARBO - GROLF - GROBSL - GROSTM
        IF (GRORT .LT. 0.15*CARBO) THEN
           IF (GROLF .GT. 0.0 .OR. GROBSL .GT. 0.0 .OR.
     &         GROSTM .GT. 0.0) THEN
              GRF   = CARBO*0.9/(GROLF+GROBSL+GROSTM)
              GRORT = CARBO*0.1
            ELSE
              GRF = 1.0
           ENDIF
         ENDIF

!-----------------------------------------------------------------
      CASE (3)
        !
        !      ! Leaf cycle 2  to Leaf cycle 3
        !

          LFWT        = ((LFWT3/P4)*DTT)     + LFWT                      
          BASLFWT     = ((BASLFWT3/P4)*DTT)  + BASLFWT                   
          STMWT       = ((STMWT3/P4)*DTT)    + STMWT
          GRORT       = ((GRORT3/P4)*DTT)    + GRORT                     
          MAXLAI      = ((LAI2/P3)*DTT)      + MAXLAI  
                                                                                   
        !
        ! Check the balance of supply and demand
        !
        GRORT = CARBO - GROLF - GROBSL - GROSTM
        IF (GRORT .LT. 0.15*CARBO) THEN
           IF (GROLF .GT. 0.0 .OR. GROBSL .GT. 0.0 .OR.
     &         GROSTM .GT. 0.0) THEN
              GRF   = CARBO*0.9/(GROLF+GROBSL+GROSTM)
              GRORT = CARBO*0.1
            ELSE
              GRF = 1.0
           ENDIF  
        ENDIF

!-----------------------------------------------------------------
C      CASE (4)
        !
        !       Leaf cycle 2  to Leaf cycle 3
        !

!        IF ((LN) .GT. (26) . AND. (LN) .LE. (39)) THEN


C          LFWT        = ((LFWT3/P4)*DTT)     + LFWT                      
C          BASLFWT     = ((BASLFWT3/P4)*DTT)  + BASLFWT                   
C          STMWT       = ((STMWT3/P4)*DTT)    + STMWT
C          GRORT       = ((GRORT3/P4)*DTT)    + GRORT                    
C          MAXLAI      = ((LAI2/P3)*DTT)      + MAXLAI  
                                                                                   
        !
        ! Check the balance of supply and demand
        !
C        GRORT = CARBO - GROLF - GROBSL - GROSTM
C        IF (GRORT .LT. 0.15*CARBO) THEN
C          IF (GROLF .GT. 0.0 .OR. GROBSL .GT. 0.0 .OR.
C     &         GROSTM .GT. 0.0) THEN
C              GRF   = CARBO*0.9/(GROLF+GROBSL+GROSTM)
C              GRORT = CARBO*0.1
C            ELSE
C              GRF = 1.0
C           ENDIF  
C        ENDIF


!-----------------------------------------------------------------

       
        
      CASE (5,6,7)     
        !
        ! Forcing to EarlyAnthesis
        !

C-----------------------------------------------------------------
C Flower and fruit growth factor
c ----------------------------------------------------------------


        GRORT  = 0.05 * GROLF
        GRORT  = AMAX1 (GRORT,0.0)
        GROFLR = (1.26-0.17*PLTPOP+0.0075*PLTPOP**2)*DTT/20.5   
     &           *AMIN1(AGEFAC,TURFAC)
        GROFLR = AMAX1 (GROFLR,0.0)
        GROSTM = CARBO - GROLF - GROBSL - GRORT - GROFLR
        IF (GROSTM .LT. 0.16*CARBO) THEN
           IF (GROLF .GT. 0.0 .OR. GROBSL .GT. 0.0 .OR. GRORT .GT. 0.0
     &         .OR. GROFLR .GT. 0.0) THEN
              GRF    = CARBO*0.84/(GROLF+GROBSL+GRORT+GROFLR)
              GROSTM = CARBO*0.16
            ELSE
              GRF    = 1.0
           ENDIF

        ENDIF

        !STMWT   = STMWT   + GROSTM
        FLRWT    = FLRWT   + GROFLR
        

        SUMP  = SUMP  + CARBO !Total biomass cumulated during the stage
        IDURP = IDURP + 1     !Duration of the stage        
!-----------------------------------------------------------------
 
      CASE (8)         
        !
        ! Fruit growth -----> Last Anthesis to Physiological maturity
        !

        IF (SWMAX .LE. 0.0) THEN
           IF (XSTAGE .GE. 8.0) THEN
              SWMAX = STMWT
              SWMIN = 0.65*SWMAX
           ENDIF
        ENDIF

        IF (ISTAGE .GT. P8) THEN            
           GROSK = CARBO*0.1
           STMWT = STMWT
           SKWT  = SKWT + GROSK
           GO TO 2400
        ENDIF


        RGFILL       = (-0.23311286 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 
     & + 4.47249682 * (LOG(GDDFR/(TMAXGRO/SRADGRO))) - 7.27511652)*0.40
        !RGFILL  = 1-0.0025*(TEMPM-26.)**2                  
        

           GROFRT = RGFILL*GPP*G3*0.001*(0.7+0.2*SWFAC+1.30)

        GROCRWN = GROFRT * Yieldfact
        GRORT   = CARBO  * YIELD1

        IF (TOTPLTWT .GT. 600.0) then
           GROSK  = CARBO*0.09
           GROSK  = AMAX1 (GROSK,0.0)
           GROSTM = CARBO - GROFRT - GROCRWN - GRORT - GROSK
           IF (GROSTM.LT.0.0) GO TO 1700
           IF (GROSTM .GT. 0.15*CARBO) THEN
              GROSTM = 0.15*CARBO
           ENDIF
            SKWT   = SKWT   + GROSK
           !STMWT  = STMWT  + GROSTM
           CRWNWT = CRWNWT + GROCRWN
           GO TO 1900
        ENDIF

        GROSTM = CARBO - GROFRT - GROCRWN - GRORT
        IF (GROSTM.LT.0.0) GO TO 1700
        IF (GROSTM .GT. 0.15*CARBO) THEN
           GROSTM = 0.15*CARBO
        ENDIF

        GO TO 1900

1700    SELECT CASE (PMTYPE)
          CASE (1:13)                         
            IF (SUMDTT .LT. P8) THEN       
               IF (SRAD .LT. 6.0) THEN
                  GROSTM  = CARBO
                ELSEIF (SRAD .LT. 13.0) THEN
                  GROSTM  = CARBO*((13.-SRAD)/7.)
                  GROFRT  = (CARBO - GROSTM)* Yieldfact
                  GROCRWN = (CARBO - GROSTM)* YIELD1
                ELSE
                   GROSTM  = 0.0
                   GROFRT  = Yieldfact  * CARBO
                   GROCRWN = YIELD1 * CARBO
               ENDIF
               STMWT  = STMWT  + GROSTM
               CRWNWT = CRWNWT + GROCRWN
             ELSE
               IF (SRAD .LT. 6.0) THEN
                  GROSTM  = CARBO
                ELSEIF (SRAD .LT. 13.0) THEN
                  GROSTM  = CARBO*((13.0-SRAD)/7.0)
                  GROFRT  = (CARBO - GROSTM)*Yieldfact
                  GROCRWN = (CARBO - GROSTM)*YIELD1
                  STMWT   = STMWT  + GROSTM
                  CRWNWT  = CRWNWT + GROCRWN
                ELSE
                  STMWT   = STMWT + CARBO - GROFRT - GROCRWN - GRORT
                  IF (STMWT .LT. SWMIN) THEN
                     STMWT   = SWMIN
                     GROFRT  = Yieldfact  * CARBO /0.9
                     GROCRWN = YIELD1  * CARBO /0.9
                     CRWNWT  = CRWNWT + GROCRWN
                  ENDIF
               ENDIF
            ENDIF
          CASE (0)
            GROSTM  = 0.0
            GROFRT  = Yieldfact  * CARBO / 0.9
            GROCRWN = YIELD1  * CARBO / 0.9
            !STMWT   = STMWT  + GROSTM
            CRWNWT  = CRWNWT + GROCRWN
        END SELECT

 1900   IF (ISWNIT .EQ. 'Y') THEN
           !
           ! Grain N allowed to vary between .01 and .018.
           ! High temp., low soil water, and high N increase grain N
           !
           SFAC  = 1.125 - 0.1250*TURFAC
           TFAC  = 0.690 + 0.0125*TEMPM
           GNP   = (0.004+0.013*NFAC)*AMAX1(SFAC,TFAC)
           NSINK = GROGRN*GNP   
 
           IF (NSINK .GT. 0.0) THEN
              RMNC   = 0.75*RCNP
              RANC   = AMAX1  (RANC,RMNC)
              VANC   = STOVN / STOVWT
              VANC   = AMAX1  (VANC,VMNC)
              NPOOL1 = STOVWT*(VANC-VMNC)
              NPOOL2 = RTWT  *(RANC-RMNC)
              XNF    = 0.15  + 0.25*NFAC
              TNLAB  = XNF   * NPOOL1
              RNLAB  = XNF   * NPOOL2
              NPOOL  = TNLAB + RNLAB
              IF (ICSDUR .EQ. 1) THEN
                 GPP = AMIN1(GPP*NDEF3,(NPOOL/(0.062*.0095)))
              ENDIF
              NSDR = NPOOL/NSINK
              IF (NSDR .LT. 1.0) THEN
                 NSINK = NSINK*NSDR
              ENDIF
              IF (NSINK .GT. TNLAB) THEN
                 STOVN = STOVN - TNLAB
                 RNOUT = NSINK - TNLAB
                 ROOTN = ROOTN - RNOUT
                 RANC  = ROOTN / RTWT
               ELSE
                 STOVN = STOVN - NSINK
                 VANC  = STOVN / STOVWT
              ENDIF
           ENDIF
 
           GRAINN = GRAINN + NSINK
        ENDIF
        !
        ! Update fruit weight
        !
        FRTWT = FRTWT + GROFRT
        CRWNWT = CRWNWT + GROCRWN
        
        IF (SUMDTT .GT. P8) THEN       !IF (SUMDTT .GT. 0.8*P4) THEN   
           !STMWT = AMIN1 (STMWT,SWMAX)
        ENDIF

!-----------------------------------------------------------------
      CASE (9)      
!        
!       Physiological maturity
!

          IF (SWMAX .LE. 0.0) THEN
           IF (XSTAGE .GE. 8.0) THEN
              SWMAX = STMWT
              SWMIN = 0.65*SWMAX
           ENDIF
        ENDIF

        IF (ISTAGE .GE. P8) THEN            
           GROSK = CARBO*0.1
           STMWT = STMWT
           SKWT  = SKWT + GROSK
           GO TO 2400
        ENDIF

        GROCRWN = Yieldfact*GROFRT
        GRORT   = CARBO*YIELD1

        IF (TOTPLTWT .GT. 600.0) then
           GROSK  = CARBO*0.09
           GROSK  = AMAX1 (GROSK,0.0)
           GROSTM = CARBO - GROFRT - GROCRWN - GRORT - GROSK
           IF (GROSTM.LT.0.0) GO TO 1700
           IF (GROSTM .GT. 0.15*CARBO) THEN
              GROSTM = 0.15*CARBO
           ENDIF
            SKWT   = SKWT   + GROSK
           !STMWT  = STMWT  + GROSTM
           CRWNWT = CRWNWT + GROCRWN
           GO TO 1900
        ENDIF

        GROSTM = CARBO - GROFRT - GROCRWN - GRORT
        IF (GROSTM.LT.0.0) GO TO 1700
        IF (GROSTM .GT. 0.15*CARBO) THEN
           GROSTM = 0.15*CARBO
        ENDIF
        CRWNWT = CRWNWT + GROCRWN
        STMWT  = STMWT  + GROSTM


        GO TO 1900




!-----------------------------------------------------------------

      END SELECT
!-----------------------------------------------------------------

      IF (CARBO .EQ. 0.0) THEN
         CARBO = 0.001                 ! Make sure that carbo is not 0.
      ENDIF
!     PDWI (g/plant/day) is potential shoot growth
      PDWI   = PCARB*(1.0-GRORT/CARBO) 
!     Pgrort is potential root growth
      PGRORT = PCARB*GRORT/CARBO       
!    
!     Calculation of zero-to-unity factors for leaf senescence due to drought
!     stress (SLFW), competition for light (SLFC), and low temperature (SLFT).
!    
 2400 SLFW = 1.0
      SLFN = 0.95+0.05*AGEFAC
      SLFC = 1.0
      IF (ISTAGE .GT. 2 .AND. ISTAGE .LT. 10) THEN   
         IF (LAI .GT. 6.0) THEN
            SLFC = 1.0-0.0005*(LAI-6.0)
         ENDIF
      ENDIF
     
      SLFT = 1.0
      IF (TEMPM .LE. 4.0) THEN
         SLFT = 1.0-(4.0-TEMPM)/4.0
      ENDIF

      IF (TMIN .GT. 0.0) THEN
         ICOLD = 0
       ELSE
         SLFT  = 0.0
         ICOLD = ICOLD + 1
      ENDIF
!      
!     Leaf area senescence on a day (PLAS) and LAI is calculated for stage 1 to 5
!      
      SLFT  = AMAX1 (SLFT,0.0)
      PLAS  = (PLA-SENLA)*(1.0-AMIN1(SLFW,SLFC,SLFT))
      SENLA = SENLA + PLAS
      SENLA = AMAX1 (SENLA,SLAN)
      SENLA = AMIN1 (SENLA,PLA)
      LAI   = (PLA-SENLA)*PLTPOP*0.0001     
        
      

      IF (LN .GT. 3 .AND .LAI .LE. 0.0 .AND. ISTAGE .LE. 5) THEN    
         WRITE (MSG(1),   2800)
         CALL WARNING(1, ERRKEY, MSG)

         ISTAGE = 7                                      
       ELSE
         IF (ICOLD .GE. 7) THEN
            WRITE (MSG(1),   2800)
            CALL WARNING(1, ERRKEY, MSG)

           ISTAGE = 8                                     
         ENDIF
      ENDIF
!      
!     Half GRORT is used for respiration and 0.5% of root is lost due to senescence
!      
      RTWT = RTWT + 0.45*GRORT - 0.0025*RTWT
!      
!     Finally, total biomass per unit area (BIOMAS g/m2), total plant weight,
!     Total plant dry weight per hectare (DM kg/ha) and Plant top fraction
!     (PTF) are calculated
!      
!     When fruit development starts, the fruit population (FRUITS) is
!       less than the plant population (PLTPOP). Need to differentiate
!       for consistency with daily and seasonal outputs.
      SELECT CASE(ISTAGE)
      CASE(8,9,10)                                            
!       In this case we need sum the crown
        BIOMAS   = (LFWT + STMWT + BASLFWT + SKWT)*PLTPOP 
     &                + (FRTWT * FRUITS) + (CRWNWT * FRUITS)
      CASE(5,6,7)                                            
!       In this case FLRWT is fruit + crown
        BIOMAS   = (LFWT + STMWT + BASLFWT + SKWT)*PLTPOP 
     &                + (FLRWT * FRUITS)
      CASE DEFAULT
        BIOMAS   = (LFWT + STMWT + FLRWT + BASLFWT + SKWT)*PLTPOP ! It is BIOMAS at forcing.
      END SELECT 

      TOTPLTWT =  LFWT + STMWT + FLRWT + BASLFWT + SKWT
      DM       = BIOMAS*10.0 
      STOVWT   = LFWT + STMWT
      PTF      = (LFWT+BASLFWT+STMWT+FLRWT+SKWT) /
     &           (LFWT+BASLFWT+STMWT+FLRWT+SKWT+RTWT)
      
      IF (ISWNIT .NE. 'N') THEN
        CALL Aloha_NUPTAK(CONTROL, ISWITCH, 
     &    ISTAGE, NO3, NH4, PDWI, PGRORT, PLIGRT,         !Input
     &    PLTPOP, PTF, RANC, RCNP, RLV, RTWT, SOILPROP,   !Input
     &    STOVWT, SW, TCNP, XSTAGE,                       !Input
     &    ROOTN, SENESCE, STOVN, TANC, UNH4, UNO3, WTNUP) !Output
      ENDIF

!-----------------------------------------------------------------
      RETURN

C-----------------------------------------------------------------------
C     Format Strings
C-----------------------------------------------------------------------

 2800 FORMAT (2X,'Crop failure growth program terminated ')

!=======================================================================
!     Integration
!-----------------------------------------------------------------------
      CASE (INTEGR)
!=======================================================================
!     This code used to be in PhaseI subroutine. Put here to make timing match 
!     with old code.
!     Some of the code was removed to other subroutines.
      IF (ISTAGE /= ISTAGE_old) THEN
        ISTAGE_OLD = ISTAGE

!       New stage initialization
!-----------------------------------------------------------------------
        SELECT CASE (ISTAGE)
!-----------------------------------------------------------------------
         CASE (12)    !Planting initial variables      
          YRDOY   = CONTROL % YRDOY
          
         
          WTINITIAL = SDWTPL/(PLTPOP*10.0)        ! kg/ha  --> g/plt
          
          IF (WTINITIAL .LE. 61.0) THEN                               ! 
                                                   ! 
          PLA        = EXP(-0.18039357 * (LOG(WTINITIAL))**2 + 2.60612446 
     & * (LOG(WTINITIAL)) + 0.11938255) ! y = -0.18039357x2 + 2.60612446x + 0.11938255 R² = 0.9802 
          PLA12      = PLA
          LAI        = PLTPOP*PLA*0.0001     !         
          LFWT       = EXP(-0.13682972 * (LOG(WTINITIAL))**2 + 2.27562675 
     & * (LOG(WTINITIAL)) - 3.34521207)  ! y = -0.13682972x2 + 2.27562675x - 3.34521207 R² = 0.9844 
          LFWT12     = LFWT 
          BASLFWT    = EXP(0.247306 * (LOG(WTINITIAL))**2 - 1.189199 
     & * (LOG(WTINITIAL)) + 3.244101)  ! y = 0.247306x2 - 1.189199x + 3.244101R² R² = 0.6950 
          BASLFWT12  = BASLFWT 
          STMWT      = EXP(0.74674847 * (LOG(WTINITIAL))**2 - 4.79095024 
     & * (LOG(WTINITIAL)) + 8.93546032)  ! y = 0.74674847x2 - 4.79095024x + 8.93546032 R² = 0.9032 
          STMWT12    = STMWT
          STOVWT     = WTINITIAL
          
          ELSE
                                                   ! 
          PLA        = EXP(-11.36267709 * (LOG(WTINITIAL))**2 + 96.88737658 
     & * (LOG(WTINITIAL)) - 198.54329897) ! y = -11.36267709x2 + 96.88737658x - 198.54329897 R² = 0.2313 
          PLA12      = PLA
          LAI        = PLTPOP*PLA*0.0001     !         
          LFWT       = EXP(-35.88863437 * (LOG(WTINITIAL))**2 + 305.89682709 
     & * (LOG(WTINITIAL)) - 647.87752188)  ! y = -35.88863437x2 + 305.89682709x - 647.87752188 R² = 0.9516 
          LFWT12     = LFWT 
          BASLFWT    = EXP(106.46751403 * (LOG(WTINITIAL))**2 - 904.07142367 
     & * (LOG(WTINITIAL)) + 1921.62542892)  ! y = 106.46751403x2 - 904.07142367x + 1921.62542892 R² = 0.9348 
          BASLFWT12  = BASLFWT 
          STMWT      = EXP(66.59416831 * (LOG(WTINITIAL))**2 - 564.70708020 
     & * (LOG(WTINITIAL)) + 1198.91190261)  ! y = 66.59416831x2 - 564.70708020x + 1198.91190261 R² = 0.9816 
          STMWT12    = STMWT
          STOVWT     = WTINITIAL
          
         ENDIF 
          BIOMAS= (LFWT + STMWT + BASLFWT)*PLTPOP
          YRPLT   = YRDOY
!-----------------------------------------------------------------------    
         CASE (13)   !Planting to Root initiation

          GRORTI  = 0.01
          RTWT    = 0.10                !Root weight         
          FLRWT   = 0.0
          FLRWT   = 0.0
          GROSTM  = 0.0                 
          SENLA   = 0.0                 
          SLAN    = 0.0                 
          GRORT   = 0.0                 
          GROBSL  = 0.0                 
          GROLF   = 0.0                 
          CUMPH   = 0.514               
          LN      = 1                   
          CUMDEP  = 0.0

          YRDOY   = CONTROL % YRDOY    
          NDOF = TIMDIF(YRPLT, YRDOY)  
          DAP1 = NDOF                  
          GDDFR   = (SUMDTTGRO - SUMDTT)/(DAP1)   
          
          TMAXGRO = (SUMTMAXGRO - SUMTMAX)/(DAP1)   
          SRADGRO = (SUMSRADGRO - SUMSRAD)/(DAP1)  
          PARGRO =  (SUMPARGRO - SUMPAR)/(DAP1)    
      
          ADJGDDF13 = (SUMDTTGRO - SUMDTT)
          ADJTMAX13 = (SUMTMAXGRO - SUMTMAX)
          ADJSRAD13 = (SUMSRADGRO - SUMSRAD)



         
          
C          IF (TEMPM .LT. 20.0) THEN                               ! 
C         TRF2 = 0.98                                          ! 
C       ELSEIF (TEMPM .GE. 20.0 .AND. TEMPM .LT. 26.0) THEN    ! 
C         TRF2 = 1.13314845306683*EXP(-0.005*TEMPM)   !1-0.005*(TEMPM-25) !
C       ELSE
C         TRF2 = 0.98                                          ! 
C      ENDIF
          
C          IF ((LOG(GDDFR/(TMAXGRO/SRADGRO))).LT. 1.0) THEN                
C         CARBO = TRF2*PCARB*AMIN1(PRFT,TRF2*SWFAC,NSTRES)

C       ELSEIF ((LOG(GDDFR/(TMAXGRO/SRADGRO))).LT. 1.93) THEN   
C         CARBO = TRF2*PCARB*AMIN1(PRFT,SWFAC,NSTRES)

C         ELSEIF ((LOG(GDDFR/(TMAXGRO/SRADGRO))).GT. 1.93) THEN   
C         CARBO = TRF2*PCARB*AMIN1(PRFT,SWFAC,NSTRES)
C       ENDIF 

C          PLA13 = LFWT * 0.5
C          LFWT13 = CARBO * 10
C          BASLFWT13 = CARBO * 10
C          STMWT13 = CARBO * 10




         RLAE13    = 0.01289488 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 0.02086888 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 0.00524122     !y = 0.01289488x2 - 0.02086888x + 0.00524122 R² = 0.8379  
         PLA13     = PLA12*EXP(RLAE13*(DAP1))                  
         
         RLDW13    = 16.64364072 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 32.63404677 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 13.62813905     !y = 16.64364072x2 - 32.63404677x + 13.62813905 R² = 0.6795 
         LFWT13    = LFWT12*EXP((RLDW13/1000)*(DAP1)) 
          
         RBWTDW13  = 81.57054008 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 339.67610159 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 350.05047327   !y = 81.57054008x2 - 339.67610159x + 350.05047327 R² = 0.5642 
         BASLFWT13 = BASLFWT12*EXP((RBWTDW13/1000)*(GDDFR/(TMAXGRO/SRADGRO)))
         
         RSTMWT13  = -24.00888219 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 + 80.30332258 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) - 42.97510499   !y = -24.00888219x2 + 80.30332258x - 42.97510499 R² = 0.1005 
         STMWT13   = STMWT12*EXP((RSTMWT13/1000)*(DAP1))
          
         LAI13    = PLTPOP*PLA13*0.0001

          LAI         = LAI13                     
          LFWT        = LFWT13                   
          BASLFWT     = BASLFWT13                
          STMWT       = STMWT13                 
          STOVWT= STMWT
          BIOMAS= (LFWT + STMWT + BASLFWT)*PLTPOP
          BIOMAS13 = BIOMAS
          LAI13 = LAI
          MAXLAI = LAI

!-----------------------------------------------------------------------                 
        CASE (1)        !Root initiation to First New Leaf
          PLAG    = 0.0

         
          GRORTI  = 0.01
          RTWT    = 0.20                !Root weight         
          FLRWT   = 0.0
          FLRWT   = 0.0
          GROSTM  = 0.0                 
          SENLA   = 0.0                 
          SLAN    = 0.0                 
          GRORT   = 0.0                 
          GROBSL  = 0.0                 
          GROLF   = 0.0                 
          CUMPH   = 0.514               
          LN      = 1                   
          CUMDEP  = 0.0
          
         
         
         YRDOY   = CONTROL % YRDOY   ! Root initiation date
         NDOF = TIMDIF(YRPLT, YRDOY)  
         DAP2 = NDOF
         DAP3     = NDOF - DAP1      
          
         ADJGDDF1 = (SUMDTTGRO - SUMDTT)
         ADJTMAX1 = (SUMTMAXGRO - SUMTMAX)
         ADJSRAD1 = (SUMSRADGRO - SUMSRAD)

         GDDFR   = (ADJGDDF13 + ADJGDDF1)/(DAP3+DAP1)
         TMAXGRO = (ADJTMAX13 + ADJTMAX1)/(DAP3+DAP1)
         SRADGRO = (ADJSRAD13 + ADJSRAD1)/(DAP3+DAP1)
          
         PARGRO =  (SUMPARGRO - SUMPAR)/(DAP3)     
         SEEDQLY = ((1-(BASLFWT12/LFWT12))*1.3) !  Seed quality adjusts LFWT1
         
         
         RLAE1    = 0.01341697 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 0.02224724 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 0.00607536     !y = 0.01341697x2 - 0.02224724x + 0.00607536 R² = 0.8349  
         PLA1     = PLA12*EXP(RLAE1*(DAP3+DAP1))                                                                                  
         
         RLDW1    = 16.64364072 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 32.63404677 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 13.62813905  !y = 16.64364072x2 - 32.63404677x + 13.62813905 R² = 0.6795 
         LFWT1    = LFWT12*EXP((RLDW1/1000)*(DAP3+DAP1))/SEEDQLY
        
         RBWTDW1  = 81.57054008 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 339.67610159 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 350.05047327   !y = 81.57054008x2 - 339.67610159x + 350.05047327 R² = 0.5642 
         BASLFWT1 = BASLFWT12*EXP((RBWTDW1/1000)*(GDDFR/(TMAXGRO/SRADGRO)))/SEEDQLY
         
         RSTMWT1  = -39.90744022 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 + 129.64536249 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) - 69.90080305   !y = -39.90744022x2 + 129.64536249x - 69.90080305 R² = 0.1494 
         STMWT1   = STMWT12*EXP((RSTMWT1/1000)*(DAP3+DAP1))/SEEDQLY                                                                 
         
         RGRORT1 = 110.22488765 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 72.34063846 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 94.37884405   !y = 110.22488765x2 - 72.34063846x + 94.37884405 R² = 0.8617  
         GRORT1  = GRORTI*EXP((RGRORT1/1000)*(DAP3+DAP1))

         
         
          LAI1    = PLTPOP*PLA1*0.0001
          
          LAI         = LAI1                       
          LFWT        = LFWT1                      
          BASLFWT     = BASLFWT1                   
          STMWT       = STMWT1                 
          STOVWT= STMWT
          BIOMAS= (LFWT + STMWT + BASLFWT)*PLTPOP
          BIOMAS1 = BIOMAS
          LAI1 = LAI 
          MAXLAI = LAI
          
!-----------------------------------------------------------------------         
             
    
      CASE (2)       ! First New Leaf to Leaf Cycle 1
          !GROSTM = 0.0  ! Daily stem growth (g/plant/day)
          YRDOY   = CONTROL % YRDOY
          NDOF = TIMDIF(YRPLT, YRDOY)
          DAP4 = NDOF
          DAP5    = NDOF-DAP2
          GDDFR   = (SUMDTTGRO - SUMDTT)/(DAP5)  
          TMAXGRO = (SUMTMAXGRO - SUMTMAX)/(DAP5)
          SRADGRO = (SUMSRADGRO - SUMSRAD)/(DAP5)
          PARGRO  = (SUMPARGRO - SUMPAR)/(DAP5)
   
             RLAE2    = -0.00042225 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 + 0.00730812 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) - 0.00551540     !y = -0.00042225x2 + 0.00730812x - 0.00551540 R² = 0.8942  
             PLA2     = PLA1*EXP(RLAE2*(DAP5))                                                                                                                                                                                                                                      
             
             RLDW2    = 0.56659753 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 + 4.55511443 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) - 2.97791272     !y = 0.56659753x2 + 4.55511443x - 2.97791272 R² = 0.8372 
             LFWT2    = LFWT1*EXP((RLDW2/1000)*(DAP5))                                                                                   
         
             RBWTDW2  = 265.95155057 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 1019.47038606 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 1044.29203812   !y = 265.95155057x2 - 1019.47038606x + 1044.29203812 R² = 0.9186 (used without ZIP1) 
             BASLFWT2 = BASLFWT1*EXP((RBWTDW2/1000)*(GDDFR/(TMAXGRO/SRADGRO)))                                                             
         
             RSTMWT2  = 252.16180637 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 1189.05973879 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 1451.11634867   !y = 252.16180637x2 - 1189.05973879x + 1451.11634867 R² = 0.8903 
             STMWT2   = STMWT1*EXP((RSTMWT2/1000)*(GDDFR/(TMAXGRO/SRADGRO)))                                                                

             RGRORT2 = 3.41064008 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 + 20.16657317 
     & * (LOG(GDDFR/(TMAXGRO/SRADGRO))) - 8.25952711   !y = 3.41064008x2 + 20.16657317x - 8.25952711 R² = 0.8636  
             GRORT2  = GRORT1*EXP((RGRORT2/1000)*(DAP5))                                                                                             

             LAI2    = PLTPOP*PLA2*0.0001
              

             LAI         = LAI2                         
             LFWT        = LFWT2                        
             BASLFWT     = BASLFWT2                     
             STMWT       = STMWT2                       
             STOVWT= STMWT
             BIOMAS= (LFWT + STMWT + BASLFWT)*PLTPOP         
             BIOMAS2 = BIOMAS
             LAI2 = LAI
             LN2  = LN 
             MAXLAI = LAI
            
!-----------------------------------------------------------------------
      
        CASE (3)       !Leaf Cycle 1 to Leaf cycle 2

          YRDOY   = CONTROL % YRDOY
          NDOF = TIMDIF(YRPLT, YRDOY)
          DAP6 = NDOF
          DAP7 = NDOF-DAP4

          GDDFR   = (SUMDTTGRO - SUMDTT)/(DAP7)  
          TMAXGRO = (SUMTMAXGRO - SUMTMAX)/(DAP7)
          SRADGRO = (SUMSRADGRO - SUMSRAD)/(DAP7)
          PARGRO  = (SUMPARGRO - SUMPAR)/(DAP7)

          
          
             RLAE3    = 0.42242459 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 1.70953080 
     &  * (LOG(GDDFR/(TMAXGRO/SRADGRO)))  + 1.80173219     !y = 0.42242459x2 - 1.70953080x + 1.80173219 R² = 0.57828588  
             PLA3     = PLA2*EXP(RLAE3*(GDDFR/(TMAXGRO/SRADGRO)))                                                                                         
             
             RLDW3    = 321.20939543 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 1369.37265384 
     &  * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 1536.95773441      !y = 321.20939543x2 - 1369.37265384x + 1536.95773441 R² = 0.7055 
             LFWT3    = LFWT2*EXP((RLDW3/1000)*(GDDFR/(TMAXGRO/SRADGRO)))                                                                               
         
             RBWTDW3  = -81.06861974 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 + 184.67892539 
     &  * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 53.80929308  !y = -81.06861974x2 + 184.67892539x + 53.80929308 R² = 0.90920630 
             BASLFWT3 = BASLFWT2*EXP((RBWTDW3/1000)*(GDDFR/(TMAXGRO/SRADGRO)))                                                             
         
             RSTMWT3  = 1286.19102531 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 4995.56061201 
     &  * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 4932.93231043   !y = 1286.19102531x2 - 4995.56061201x + 4932.93231043 R² = 0.74936713 
             STMWT3   = STMWT2*EXP((RSTMWT3/1000)*(GDDFR/(TMAXGRO/SRADGRO)))                                                                

             RGRORT3 = 1946.55069451 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 7623.83657710 
     &  * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 7519.46534541   !y = 1946.55069451x2 - 7623.83657710x + 7519.46534541 R² = 0.5264 
             GRORT3  = GRORT2*EXP((RGRORT3/1000)*(GDDFR/(TMAXGRO/SRADGRO)))                                                                                            

             LAI3    = PLTPOP*PLA3*0.0001
              

             LAI         = LAI3                         
             LFWT        = LFWT3                        
             BASLFWT     = BASLFWT3                     
             STMWT       = STMWT3                       
             STOVWT= STMWT
             BIOMAS= (LFWT + STMWT + BASLFWT)*PLTPOP  
             BIOMAS3 = BIOMAS
             LAI3 = LAI
             LN3  = LN 
             MAXLAI = LAI

!-------------------------------------------------------- 
         CASE (4)           !Leaf cycle 2 to Leaf cycle 3

             
          PLANTSIZE = TOTPLTWT
          
          YRDOY   = CONTROL % YRDOY
          NDOF = TIMDIF(YRPLT, YRDOY)
          DAP8 = NDOF
          DAP9 = NDOF-DAP6

          GDDFR   = (SUMDTTGRO - SUMDTT)/(DAP9)  
          TMAXGRO = (SUMTMAXGRO - SUMTMAX)/(DAP9)
          SRADGRO = (SUMSRADGRO - SUMSRAD)/(DAP9)
          PARGRO  = (SUMPARGRO - SUMPAR)/(DAP9)
          
             RLAE4    = 0.02279067 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 0.09618739 
     &  * (LOG(GDDFR/(TMAXGRO/SRADGRO)))  + 0.15232236     !y = 0.42157206x2 - 1.70530517x + 1.79648701 R² = 0.4883  
             PLA4     = PLA3*EXP(RLAE4*(GDDFR/(TMAXGRO/SRADGRO)))
          
             RLDW4    = 200.57284063 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 774.01981366 
     &  * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 804.12791119      !y = 200.57284063x2 - 774.01981366x + 804.12791119 R² = 0.4563 
             LFWT4    = LFWT3*EXP((RLDW4/1000)*(GDDFR/(TMAXGRO/SRADGRO)))

             RBWTDW4  = 204.37991743 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 788.67484148 
     &  * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 816.44955507  !y = 204.37991743x2 - 788.67484148x + 816.44955507 R² = 0.30922162 
             BASLFWT4 = BASLFWT3*EXP((RBWTDW4/1000)*(GDDFR/(TMAXGRO/SRADGRO)))                                                             
         
             RSTMWT4  = 443.62453629 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 - 1745.81995161 
     &  * (LOG(GDDFR/(TMAXGRO/SRADGRO))) + 1821.34841003   !y = 443.62453629x2 - 1745.81995161x + 1821.34841003 R² = 0.28684196 
             STMWT4   = STMWT3*EXP((RSTMWT4/1000)*(GDDFR/(TMAXGRO/SRADGRO)))
               
             RGRORT4 = -380.88231019 * (LOG(GDDFR/(TMAXGRO/SRADGRO)))**2 + 1533.01000755 
     &  * (LOG(GDDFR/(TMAXGRO/SRADGRO))) - 1424.86934789   !y = -380.88231019x2 + 1533.01000755x - 1424.86934789 R² = 0.32845892 
             GRORT4  = GRORT3*EXP((RGRORT4/1000)*(GDDFR/(TMAXGRO/SRADGRO)))                                                                                            

             LAI4    = PLTPOP*PLA4*0.0001
  
             LAI         = LAI4                         
             LFWT        = LFWT4                        
             BASLFWT     = BASLFWT4                     
             STMWT       = STMWT4                       
             STOVWT= STMWT
             BIOMAS= (LFWT + STMWT + BASLFWT)*PLTPOP  
             BIOMAS4 = BIOMAS
             LAI4 = LAI
             LN4  = LN
  
!--------------------------------------------------------   

         CASE (5)                           

             PLANTSIZE = TOTPLTWT

          
          YRDOY   = CONTROL % YRDOY
          NDOF = TIMDIF(YRPLT, YRDOY)
          DAP10 = NDOF
          DAP11 = NDOF-DAP8

          GDDFR   = (SUMDTTGRO - SUMDTT)/(DAP11)  
          TMAXGRO = (SUMTMAXGRO - SUMTMAX)/(DAP11)
          SRADGRO = (SUMSRADGRO - SUMSRAD)/(DAP11)
          PARGRO  = (SUMPARGRO - SUMPAR)/(DAP11)
          
                  
          FBIOM  = BIOMAS               
          SUMP   = 0.0                  
          IDURP  = 0                    
          PLAMX  = PLA                  
          GROFLR = 0.0
          GROCRWN= 0.0
          GROFRT = 0.0
          FLRWT  = 0.0
          FRTWT  = 0.0
          CRWNWT = 0.0
          LAI     = AMAX1(LAI1, LAI) 
          LAI     = AMAX1(LAI2, LAI) 
          LAI     = AMAX1(LAI3, LAI) 
          LAI     = AMAX1(LAI4, LAI) 
          BIOMAS= (LFWT + STMWT + BASLFWT + FLRWT)*PLTPOP
         
!---------------------------------------------------------NEW END          


       CASE (6)                  
          YRDOY   = CONTROL % YRDOY
          NDOF = TIMDIF(YRPLT, YRDOY)
          DAP12 = NDOF
          DAP13 = NDOF-DAP10

          GDDFR   = SUMDTTGRO/DAP13  
          TMAXGRO = SUMTMAXGRO/DAP13
          SRADGRO = SUMSRADGRO/(DAP13)
          PARGRO  = SUMPARGRO/(DAP13)


C         Average photosysnthesis rate of fruit eye
          PHOTOSYNEYE = SUMP*1000./IDURP  

C         G2 is genetic coefficient for potential eye number
          GPP    = G2*(PHOTOSYNEYE/12000+0.43)*
     &             (0.7+0.3*PLANTSIZE/G2)          
          GPP    = AMIN1 (GPP,G2)                
          GPP    = AMAX1 (GPP,0.0)




!CHP 10/14/2017          FLRWT  =  0.1*STMWT           ! FLRWT stands for the weight ofwhole inflorescence STMWT is stem weight.  Both are in gram/plant.
          SKWT   =  0.0
          GROSK  =  0.0
          PTF    =  1.0                 
          EYEWT  =  0.0                 
          VANC   = TANC                 
          VMNC   = TMNC          

          LAI     = AMAX1(LAI1, LAI) 
          LAI     = AMAX1(LAI2, LAI) 
          LAI     = AMAX1(LAI3, LAI) 
          LAI     = AMAX1(LAI4, LAI)
          BIOMAS= (LFWT + STMWT + BASLFWT + FLRWT)*PLTPOP   
!-----------------------------------------------------------------------
   

          CASE (7)       

          YRDOY   = CONTROL % YRDOY
          NDOF = TIMDIF(YRPLT, YRDOY)
          DAP14 = NDOF
          DAP15 = NDOF-DAP12

          GDDFR   = (SUMDTTGRO - SUMDTT)/(DAP15)  
          TMAXGRO = (SUMTMAXGRO - SUMTMAX)/(DAP15)
          SRADGRO = (SUMSRADGRO - SUMSRAD)/(DAP15)
          PARGRO  = (SUMPARGRO - SUMPAR)/(DAP15)


C         Average photosysnthesis rate of fruit eye
          PHOTOSYNEYE = SUMP*1000./IDURP  

C         G2 is genetic coefficient for potential eye number
          GPP    = G2*(PHOTOSYNEYE/12000+0.43)*
     &             (0.7+0.3*PLANTSIZE/G2)          
          GPP    = AMIN1 (GPP,G2)                
          GPP    = AMAX1 (GPP,0.0)



          SKWT   =  0.0
          GROSK  =  0.0
          PTF    =  1.0                 
          EYEWT  =  0.0                 
          VANC   = TANC                 
          VMNC   = TMNC          


             LAI     = AMAX1(LAI1, LAI) 
             LAI     = AMAX1(LAI2, LAI) 
             LAI     = AMAX1(LAI3, LAI) 
             LAI     = AMAX1(LAI4, LAI)
             BIOMAS= (LFWT + STMWT + BASLFWT + FLRWT)*PLTPOP
!-----------------------------------------------------------------------              
 
        CASE (8)                
          YRDOY   = CONTROL % YRDOY
          NDOF = TIMDIF(YRPLT, YRDOY)
          DAP16 = NDOF
          DAP17 = NDOF-DAP14

         
          GDDFR   = (SUMDTTGRO - SUMDTT)/(DAP17)  
          TMAXGRO = (SUMTMAXGRO - SUMTMAX)/(DAP17)
          SRADGRO = (SUMSRADGRO - SUMSRAD)/(DAP17)
          PARGRO  = (SUMPARGRO - SUMPAR)/(DAP17)

          GDDFRS8   = GDDFR
          TMAXGROS8 = TMAXGRO
          SRADGROS8 = SRADGRO
          PARGROS8  = PARGRO


             LAI     = AMAX1(LAI1, LAI) 
             LAI     = AMAX1(LAI2, LAI) 
             LAI     = AMAX1(LAI3, LAI) 
             LAI     = AMAX1(LAI4, LAI)  

 
          FRUITS = PLTPOP*(1.-0.10*PLTPOP/65.0)  
C         There will be some loss of mass when going from flower mass  !  
C           to fruit + crown because FRUITS (#/m2) < PLTPOP (#/m2)     !  
C                                                                      ! 
C FRTWT (g/plant) is fruit weight.  It is assumed to be 50% of inflorescence at begining of the stage
C CRWNWT (g/plant) is crown weight which is assumed to be 20% of inflorescence at the begining of the stage
C         FRTWT  = FLRWT*0.5            
C         CRWNWT = FLRWT*0.2            
C 10/14/2017 CHP 50% to fruit and 20% to crown causes 30% of flower mass to be lost. 
C     change ratios to add up to 1, maintaining approximately the same ratio.
          FRTWT  = FLRWT*1       
          CRWNWT = FLRWT*0.14    

          SWMAX  = 0.0
          SWMIN  = 0.0
          BIOMAS= (LFWT + STMWT + BASLFWT + FLRWT)*PLTPOP
          
        CASE (9)       
          
          YRDOY   = CONTROL % YRDOY
          NDOF = TIMDIF(YRPLT, YRDOY)
          DAP18 = NDOF
          DAP19 = NDOF-DAP16

          GDDFR   = (SUMDTTGRO - SUMDTT)/(DAP19)  
          TMAXGRO = (SUMTMAXGRO - SUMTMAX)/(DAP19)
          SRADGRO = (SUMSRADGRO - SUMSRAD)/(DAP19)
          PARGRO  = (SUMPARGRO - SUMPAR)/(DAP19)
         
          GDDFRS9   = GDDFR
          TMAXGROS9 = TMAXGRO
          SRADGROS9 = SRADGRO
          PARGROS9  = PARGRO

 
             LAI     = AMAX1(LAI1, LAI) 
             LAI     = AMAX1(LAI2, LAI) 
             LAI     = AMAX1(LAI3, LAI) 
             LAI     = AMAX1(LAI4, LAI)
             !BIOMAS= (LFWT + STMWT + BASLFWT + FLRWT)*PLTPOP  

        
            STGDOY (ISTAGE) = YRDOY

!-----------------------------------------------------------------------

       CASE (10)       
          
          YRDOY   = CONTROL % YRDOY
          NDOF = TIMDIF(YRPLT, YRDOY)
          DAP20 = NDOF
          DAP21 = NDOF-DAP18
          !DAP21 = (DAP13 + DAP15 + DAP17 + DAP19) 

           STGDOY (ISTAGE) = YRDOY

             LAI     = AMAX1(LAI1, LAI) 
             LAI     = AMAX1(LAI2, LAI) 
             LAI     = AMAX1(LAI3, LAI) 
             LAI     = AMAX1(LAI4, LAI)  

c         YIELD = FRTWT*10.0*FRUITS                         ! Smooth Cayenne yield 'only fruit weight'
          YIELD = (FRTWT*10.0*FRUITS) + (CRWNWT*10.0*FRUITS) ! Only fruit but MD-2 yield is fruit + crown
             
  
!-----------------------------------------------------------------------
        CASE (11) 
        YRDOY   = CONTROL % YRDOY

          IF (ISWNIT .NE. 'N') THEN
             IF (FRTWT .GT. 0.0) THEN
                XGNP = (GRAINN/FRTWT)*100.0
                GNUP = GRAINN*FRUITS*10.0
             ENDIF
             TOTNUP = GNUP + APTNUP
          ENDIF

       
!-----------------------------------------------------------------------
        END SELECT
!-----------------------------------------------------------------------        
      ENDIF

      CALL Aloha_NUPTAK (CONTROL, ISWITCH, 
     &    ISTAGE, NO3, NH4, PDWI, PGRORT, PLIGRT,         !Input
     &    PLTPOP, PTF, RANC, RCNP, RLV, RTWT, SOILPROP,   !Input
     &    STOVWT, SW, TCNP, XSTAGE,                       !Input
     &    ROOTN, SENESCE, STOVN, TANC, UNH4, UNO3, WTNUP) !Output

!=======================================================================
      END SELECT
!=======================================================================
      RETURN
      END SUBROUTINE Aloha_GROSUB

! PLAG (cm^2) is daily green leaf area growth    ! O sea debo modificar PLAG para hacer que aumentar el area foliar 
! leaf area index (m2 leaf/m2 ground)
! LFWT (g/plant) is green leaf weight which is assumed to be 53% of initial crown weight
! RTWT (g/plant) is root weight
! STMWT is 115% of initial crown weight
! Basal white leaf weight is 66% of green leaf weight
! Inflorescence weight is set to 0.0
! STOVWT (g/plant) is stover weight
! GROSTM (g/plant/day) is daily stem growth
! SENLA (cm2/plant) is area of leaf senesces due to stress on a given day
! SLAN (cm2/plant) is total normal leaf senescence since emergence.
! GRORT (g/plant/day) is daily root growth
! GROBSL (g/plant/day) is daily basal leaf growth
! GROLF (g/plant/day) is daily green leaf growth
! CUMPH (leaves/plant) is number of leaves emerged
! LN (leaves/plant) is leaf number
!PEYEWT =  Eye weight (mg/eye)
!GPSM   =  Number of eyes per square meter
!STOVER =  Total plant weight except fruit
!YIELD  =  Dry fruit yield (kg/ha)
!YIELDB =  Fresh fruit yield (lb/acre)
!LAI     = leaf area index (m2 leaf/m2 ground)
 
!LFWT    = LFWT (g/plant) is green leaf weight which is assumed to be 53% of initial crown weight
!BASLFWT = Basal white leaf weight is 66% of green leaf weight
!STMWT   = STMWT is 115% of initial crown weight
!STOVWT  = STOVWT (g/plant) is stover weight
!FBIOM  =  Record biomass at forcing
!SUMP   =  SUMP is the total weight of biomass cumulated in Istage 4.
!IDURP  =  Duration of stage 3 (days)
!PLAMX  =  PLAMX (cm2/plant) is maximal green leaf area.  PLA is total green leaf area.
! PTF is plant top fraction in gram/plant.
! EYEWT (G/eye) is weight of the eye
!FRUITS =  number of fruits=PLTPOP/m2*FRUITING%
!YIELD = fruit dry weight yield (kg/ha)
! TURFAC      !Soil water stress effect on expansion (0-1), 1 is no stress, 0 is full stress
! PC          !Used to compute fraction of phyllochron interval occurring today      


