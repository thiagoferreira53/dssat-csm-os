C=======================================================================
! combination of previous PHENOL and INPHEN subroutines
!  SUBROUTINE Aloha2_PHENOL
C  Determines phenological stage
C-----------------------------------------------------------------------
C  Revision history
C
C  1. Written
C  2. Header revision and minor changes           P.W.W.      2-7-93
C  3. Added switch block, code cleanup            P.W.W.      2-7-93
C  4. Modified TT calculations to reduce line #'s P.W.W.      2-7-93
C  5. Modified for MILLET model                   W.T.B.      MAY 94
C  6. Stages changes for inclusion in Overview    J.V.J.      9-5-20      
C=======================================================================

      SUBROUTINE Aloha2_PHENOL (CONTROL, ISWITCH,
     &    SW, WEATHER, SOILPROP, YRPLT, SUMDTTGRO, SUMTMAXGRO,
     &    SUMTMAX, DTT, EDATE, ISDATE, ISTAGE, MDATE, PMDATE, 
     &    SUMSRADGRO, SUMSRAD, SUMPARGRO, SUMPAR, STGDOY, SUMDTT, 
     &    TBASE, TEMPM, XSTAGE, EDATE12, EDATE13, EDATE1, EDATE2,
     &    EDATE3, EDATE5, EDATE6, EDATE7)               

      USE Aloha2_mod
      IMPLICIT    NONE
      SAVE

      INTEGER     STGDOY(20),YRDOY,I,NDAS,L,L0, TIMDIF, YRPLT

      REAL        TTMP,SWSD,XLAT,ROOTINGTIME

!     REAL        YIELDB,PHOTOSYNEYE,PEYEWT,LAI, BIOMAS, MAXLAI, SUMP
!     INTEGER     IDURP, ICSDUR
!     REAL        STMWT, APTNUP, RTDEP, 
!     REAL        FRUITS, SWMAX, SWMIN, YIELD, EYEWT, GPSM, STOVER
!     REAL        FDMC, HBIOM, XGNP, GNUP, TOTNUP
!     REAL        CSD1, CSD2, CNSD1, CNSD2
!     REAL, DIMENSION(NL) :: FBIOM
!     REAL, DIMENSION(20) :: SI1, SI2, SI3, SI4

      INTEGER      DYNAMIC, EDATE, MDATE, HAREND, EDATE12, EDATE13
      INTEGER      EDATE1, EDATE2, EDATE3, EDATE5, EDATE6, EDATE7
      REAL         XSTAGE
!TEMP      REAL         GRAINN

      CHARACTER*1 ISWWAT, IDETO, ISWNIT
      INTEGER     ISTAGE, NLAYR, NOUTDO, ISDATE, FHDATE, PMDATE
      REAL        TBASE
      REAL        DTT, TEMPM
!      REAL        TBASV, TOPTV, TTOPV, TBASR, TOPTR, TTOPR
      REAL        TMFAC1(8)
      REAL        TMIN, TMAX, TEMPFMX, SUMDTT, CUMDEP, GPP, SRAD, PAR
      REAL        FRTWT, TEMPFM, TOTPLTWT
      REAL        TC, P1, P2, P3, P4, P5, P6, P7, P8, G1
      REAL        TBASE1, TBASE2
      REAL        CUMDTT, SUMDTTGRO, SUMTMAX, SUMTMAXGRO, SUMSRADGRO
      REAL        SUMSRAD, SUMPARGRO, SUMPAR
      REAL, DIMENSION(NL) :: SW, LL, DLAYR

      
      REAL PLTPOP, SDEPTH, PLANTSIZE
      INTEGER NFORCING, NDOF

      TYPE (CONTROLTYPE) CONTROL
      TYPE (SwitchType) ISWITCH
      TYPE (SOILTYPE) SOILPROP
      TYPE (WEATHERTYPE) WEATHER

      DYNAMIC = CONTROL % DYNAMIC
      YRDOY   = CONTROL % YRDOY

      TMIN = WEATHER % TMIN
      TMAX = WEATHER % TMAX
      
      SRAD = WEATHER % SRAD
      PAR  = 0.5*SRAD

      LL    = SOILPROP % LL
      NLAYR = SOILPROP % NLAYR
      DLAYR = SOILPROP % DLAYR

!  ISTAGE Definition
!      11 Start simulation to planting
!      12 Planting to Root Initiation
!      13 Root Initiation to First New Leaf
!       1 First new leaf emergence to foliar cycle 1
!   2,3,4 Foliar cycle 1 to foliar cycle 2,3 and forcing 
!       5 Forcing to Open Heart
!       6 Open Heart to Early Anthesis
!       7 Early Anthesis to Last Anthesis
!       8 Last Anthesis to Physiological maturity
!       9 Physiology to Harvest
!      10 Harvest


!=================================================================
      SELECT CASE(DYNAMIC)
!=================================================================
      CASE (RUNINIT, SEASINIT)
!-----------------------------------------------------------------

      XLAT = WEATHER % XLAT
      ISWWAT = ISWITCH % ISWWAT
      ISWNIT = ISWITCH % ISWNIT


      ISTAGE = 11                                   
      XSTAGE = 0.1

      STGDOY(14) = CONTROL%YRSIM                    
      MDATE      = -99
      HAREND     = -99
      EDATE      = 9999999

      TBASE      = 13.0
      !TBASV = SPECIES % TBASV
      !TOPTV = SPECIES % TOPTV
      !TTOPV = SPECIES % TTOPV
      !TBASR = SPECIES % TBASR
      !TOPTR = SPECIES % TOPTR
      !TTOPR = SPECIES % TTOPR

!TEMP
!      IF (ISWNIT .NE. 'Y') THEN
!         TANC = 0.0
!      ENDIF
!      
!      ! Calculate initial SEED N
!      !
!      SEEDNI = (ROOTN+STOVN+GRAINN+SEEDN)*PLTPOP

      DO I = 1, 8
         TMFAC1(I) = 0.931 + 0.114*I-0.0703*I**2+0.0053*I**3
      END DO

      SDEPTH   = PLANTING % SDEPTH
      NFORCING = PLANTING % NFORCING
      NDOF     = PLANTING % NDOF
      PLTPOP   = PLANTING % PLTPOP
      
      
      TC = Cultivar % TC
      P1 = Cultivar % P1
      P2 = Cultivar % P2
      P3 = Cultivar % P3
      P4 = Cultivar % P4
      P5 = Cultivar % P5
      P6 = Cultivar % P6
      P7 = Cultivar % P7
      P8 = Cultivar % P8
      G1 = Cultivar % G1
      
      TBASE1  = 13. 
      TBASE2  = 13.

!=================================================================
      CASE (RATE)
!-----------------------------------------------------------------

!moved to grosub      XANC   = TANC*100.0               ! Top actual N concentration (g N/g Dry weight)
!moved to grosub      APTNUP = STOVN*10.0*PLTPOP
!from FileX           SDEPTH = 5.0
      
      DTT    = TEMPM - TBASE
      

      SELECT CASE (ISTAGE)

c !      CASE (2,3,4)
c          IF (YRDOY .GT. ISDATE) THEN
c                DTT = 0
c                ELSE
c                TEMPM = 0.6*TMIN+0.4*TMAX  
c          ENDIF
      
      
        CASE (1,2,3, 11,12,13)        
 
        IF (TMIN .GT. TBASE .AND. TMAX .LT. 38.0) THEN
             IF (XLAT .LT. 21.0 .and. XLAT .GT. -21.0) THEN
                TEMPM = 0.6*TMIN+0.4*TMAX
              ELSE
                TEMPM = (TMAX+TMIN)/2
             ENDIF
             DTT = TEMPM - TBASE
             ELSEIF (TMIN .LE. TBASE .OR. TMAX .GE. 38.0) THEN 
             IF (TMAX .LT. TBASE) THEN
             DTT = 0.0
                ENDIF

                          
c            IF (XLAT .LT. 21.0 .and. XLAT .GT. -21.0) THEN

c               IF (TMIN .GT. TBASE .AND. TMAX .LT. 33.5) THEN
c                    TEMPM = 0.6*TMIN+0.4*TMAX
c                  ELSE
c                  TEMPM = 0.9*TMIN+0.1*TMAX
c                ENDIF
c             DTT = TEMPM - TBASE
c             ELSEIF (TMIN .LE. TBASE .OR. TMAX .GE. 33.5) THEN 
c             IF (TMAX .LT. TBASE) THEN
c             DTT = 0.0
c                ENDIF

             IF (DTT .NE. 0.0) THEN                          
                DTT = 0.0
                DO I = 1, 8                                  
                   TTMP = TMIN + TMFAC1(I)*(TMAX-TMIN)       
                   IF (TTMP .GT. TBASE .AND. TTMP .LE. 31.0) THEN
                      DTT = DTT + (TTMP-TBASE)/8.0
                      ENDIF
                   IF (TTMP .GT. 31.0 .AND. TTMP .LT. 45.0) THEN
                      DTT = DTT + 
     &                 (31.0-TBASE)*(1.0-(TTMP-31.0)/10.0)/8.
                      ENDIF
                END DO
             ENDIF
             ENDIF
           
!-----------------------------------------------------------------
!       Reproductive Phase

        CASE (4,5,6,7,8,9,10)
          IF (TMAX .LT. TBASE) THEN     
             DTT = 0.0
          ENDIF
          IF (DTT .GT. 0.0) THEN
             
!            Correcting fruit temperature and higher temperature effect
             IF (TMAX .GT. 20.0 .AND. TMAX .LT. 33.0) THEN
                TEMPFMX = 4.32*EXP(0.078*TMAX)
              ELSEIF (TMAX .GE. 33.0 .AND. TMAX .LT. 50.0) THEN
                TEMPFMX = TMAX*(1.715-(TMAX-33.0)/35.3)
              ELSEIF (TMAX .GE. 50.0) THEN
                TEMPFMX = 62.0
              ELSE
                TEMPFMX = TMAX
             ENDIF
             IF (TMIN .GT. TBASE .AND. TEMPFMX .LT. 42.0) THEN
                IF (XLAT .LT. 21.0 .AND. XLAT .GT. -21.0) THEN
                   TEMPFM = 0.6*TMIN+0.4*TEMPFMX
                 ELSE
                   TEMPFM = (TEMPFMX+TMIN)/2
                ENDIF
                DTT = TEMPFM-TBASE
                GO TO 20
          
             ENDIF

             IF (TEMPFMX .LT. TBASE) THEN
                DTT = 0.0
             ENDIF

             IF (DTT .GT. 0.0) THEN
                DTT = 0.0
                DO I = 1, 8
                   TTMP = TMIN + TMFAC1(I)*(TEMPFMX-TMIN)
                   IF (TTMP .GT. TBASE .AND. TTMP .LE. 42.0) THEN
                      DTT = DTT + (TTMP-TBASE)/8.0
                    ELSEIF (TTMP .GT. 42.0 .AND. TTMP .LT. 62.0) THEN
                      DTT = DTT + (42.0-TBASE)*(1.0-((TTMP-42.0)/
     &                      (62.0-42.0)))/8.0
                    ELSE
                      DTT = DTT
                   ENDIF
                END DO
             ENDIF
          ENDIF
           END SELECT

    


   20 SUMDTT  = SUMDTT  + DTT
      SUMTMAX = SUMTMAX + TMAX
      SUMSRAD = SUMSRAD + SRAD
      SUMPAR  = SUMPAR  + PAR
      

      
!-----------------------------------------------------------------
!  ISTAGE Definition
!    11  Start simulation to planting
!    12  Planting to Root Initiation
!    13  Root Initiation to First New Leaf
!     1  First New Leaf to Ciclo 1,
! 2,3,4  Foliar cycle 1 to foliar cycle 2,3 and forcing
!     5  Forcing to Open Heart 
!     6  Open heart to EarlyAnthe
!     7  EarlyAnthe to LastAnthe
!     8  LastAnthe to Physiological maturity
!     9   Physiological maturity to Harvest
!    10   Harvest
!-----------------------------------------------------------------             
      SELECT CASE (ISTAGE)
!-----------------------------------------------------------------
        CASE (11)         
          !
          ! Stage 11 >> Preplanting
          !
          STGDOY(ISTAGE) = YRDOY
          NDAS           = 0
         
 !        CALL PHASEI (ISWWAT,ISWNIT)
          
          SUMDTTGRO= SUMDTT
          SUMTMAXGRO= SUMTMAX
          SUMSRADGRO= SUMSRAD
          SUMPARGRO= SUMPAR

          ISTAGE = 12         
          SUMDTT =  DTT
          SUMTMAX= TMAX
          SUMSRAD= SRAD
          SUMPAR = PAR
          

c !          IF (ISWWAT .EQ. 'N') RETURN
c           CUMDEP = 0.0
c           DO L = 1, NLAYR
c              CUMDEP = CUMDEP + DLAYR(L)
c              IF (SDEPTH .LT. CUMDEP) EXIT
c           END DO
c           L0 = L
c           RETURN

          
          
!-----------------------------------------------------------------
        CASE (12)    
          !
          ! Stage 12 >> Planting to root initiation
          !

c !         Check for soil too dry for rooting
c           IF (ISWWAT .NE. 'N') THEN
c              IF (SW(L0) .LE. LL(L0)) THEN
c                  SWSD = (SW(L0)-LL(L0))*0.65+(SW(L0+1)-LL(L0+1))*0.35
c                  NDAS = NDAS + 1
c                  IF (SWSD .LT. 0.02) RETURN
c              ENDIF
c           ENDIF

!         After 140 days, give up
          IF (NDAS .GT. 140) THEN  
             ISTAGE = 13       
             PLTPOP = 0.0 
             GPP    = 1.0
             FRTWT  = 0.0
             WRITE (     *,1399)
             IF (IDETO .EQ. 'Y') THEN
                WRITE (NOUTDO,1399)
             ENDIF
            RETURN
          ENDIF
           
          IF (SUMDTT .LT. (TC)) THEN
             
              RETURN                      
          ENDIF          
          ROOTINGTIME = SUMDTT / TBASE 
          SUMDTTGRO= SUMDTT            
          SUMTMAXGRO= SUMTMAX
          SUMSRADGRO= SUMSRAD
          SUMPARGRO= SUMPAR
          STGDOY(ISTAGE) = YRDOY
          EDATE12 = YRDOY
          
          !        CALL PHASEI (ISWWAT,ISWNIT)

          ISTAGE =  13                  
          
          SUMDTT =  DTT                 ! Cumulative growing degree days set to 0.0 
          CUMDTT  = 0.0                 ! CUMDTT is also cumulative growing degree days but it is set to 0.0 only at root initiation 
          TBASE  = 13.0                 
          SUMTMAX= TMAX
          SUMSRAD= SRAD
          SUMPAR = PAR
          RETURN

!-----------------------------------------------------------------
        CASE (13) 
          !
          ! Stage 13 >> Root initiation to first new leaf emergence
          !
          NDAS   = NDAS + 1
 !
          IF (SUMDTT .LT. (P1)) THEN   
             RETURN                       
         
             ENDIF          
         
                 
          STGDOY(ISTAGE) = YRDOY            
          EDATE = YRDOY                   
          EDATE13 = YRDOY
          SUMDTTGRO= SUMDTT               
          SUMTMAXGRO= SUMTMAX
          SUMSRADGRO= SUMSRAD
          SUMPARGRO= SUMPAR
          !        CALL PHASEI (ISWWAT,ISWNIT)
          ISTAGE  = 1
          TBASE   = TBASE1              ! Tbase1 used for calibration
          SUMDTT  = DTT                           
          SUMTMAX= TMAX
          SUMSRAD= SRAD
          SUMPAR = PAR
          RETURN

!-----------------------------------------------------------------
      CASE (1)            
          !
          ! Stage 1 >> First new leaf emergence to foliar cycle 1
          !
          NDAS   = NDAS + 1
           
          IF (YRDOY .EQ. PLANTING % ForcingYRDOY
     7     .OR. (NDAS) .GE. 650) THEN
             GO TO 21               
            ELSE  
             IF (SUMDTT .LT. (P2)) THEN  
          RETURN 
             ENDIF
          
             ENDIF 
!         Ready for next stage                                               

          STGDOY(ISTAGE) = YRDOY
          EDATE = YRDOY 
          EDATE1 = YRDOY
          SUMDTTGRO= SUMDTT              
          SUMTMAXGRO= SUMTMAX
          SUMSRADGRO= SUMSRAD
          SUMPARGRO= SUMPAR

          ISTAGE = 2
          
          TBASE  = TBASE2                 
          SUMDTT =  DTT                 
          SUMTMAX= TMAX
          SUMSRAD= SRAD
          SUMPAR = PAR
           RETURN 
!-----------------------------------------------------------------
!-----------------------------------------------------------------
        CASE (2) 
          !
          !  Stage 2 >>   foliar cycle 1 to foliar cycle 2
          ! 
          NDAS   = NDAS + 1 
          
         IF (YRDOY .EQ. PLANTING % ForcingYRDOY 
     &    .OR. (NDAS) .GE. 650)  THEN
             GO TO 21
           ELSE
             IF (SUMDTT .LT. (P3) ) THEN
          
              RETURN    
           
          ENDIF
               
          ENDIF

!         Ready for next stage
 
          STGDOY(ISTAGE) = YRDOY
          EDATE = YRDOY                  
          EDATE2 = YRDOY 
          SUMDTTGRO= SUMDTT              
          SUMTMAXGRO= SUMTMAX
          SUMSRADGRO= SUMSRAD
          SUMPARGRO= SUMPAR

          ISTAGE = 3
      
          TBASE  = TBASE1               
          SUMDTT =  DTT                
          SUMTMAX= TMAX
          SUMSRAD= SRAD
          SUMPAR = PAR
!----------------------------------------------------------------- 
          
          CASE (3) 
          !
          ! Stage 3 >>   Foliar cycle 2 to foliar cycle 3
          !
          NDAS   = NDAS + 1
          
          IF (YRDOY .EQ. PLANTING % ForcingYRDOY .OR. 
     &     (NDAS) .GE. 650) THEN  
             GO TO 21                                                        
                                                                             
          ELSE   
             IF (SUMDTT .LT. (P4) ) THEN
          
              RETURN

             ENDIF
          
          ENDIF

!         Ready for next stage
          STGDOY(ISTAGE) = YRDOY
          EDATE3 = YRDOY

          ISTAGE = 4
          SUMDTTGRO= SUMDTT               
          SUMTMAXGRO= SUMTMAX
          SUMSRADGRO= SUMSRAD
          SUMPARGRO = SUMPAR
          TBASE  = TBASE1                 
          SUMDTT =  DTT                 
          SUMTMAX= TMAX
          SUMSRAD= SRAD
          SUMPAR = PAR       
 !----------------------------------------------------------------- 
          
 
        CASE (4)       
          !
          !  Stage 4 >> Foliar Cycle 3 growth to forcing 
          !
          NDAS   = NDAS + 1
         
                
          IF (PLANTING % NFORCING .GE. 2) THEN  !
           !NDOF = TIMDIF(YRPLT, PLANTING % ForcingYRDOY) -ROOTINGTIME              
           NDOF = TIMDIF(YRPLT, PLANTING % ForcingYRDOY) -
     &      FLOOR (ROOTINGTIME) + 1  
                                                                                   
                                                                                   
                                                                                   
            
        ENDIF
            
          IF (NFORCING .GE. 2) THEN
             !
             ! Forcing by number of days after planting
             !

             IF (YRDOY .LT. PLANTING % ForcingYRDOY) THEN
                RETURN
             ENDIF

           ELSE
              !
              ! Forcing by Plant Size (200 to 350 grams usually)
              !
              IF (TOTPLTWT .LT. PLANTSIZE) THEN
                 RETURN
              ENDIF
          ENDIF
21        STGDOY(ISTAGE) = YRDOY        

          ISTAGE = 4                    
          SUMDTTGRO= SUMDTT               
          SUMTMAXGRO= SUMTMAX
          SUMSRADGRO= SUMSRAD
          SUMPARGRO = SUMPAR
          TBASE  = TBASE1                 
          SUMDTT =  DTT                 
          SUMTMAX= TMAX
          SUMSRAD= SRAD
          SUMPAR = PAR
          ISDATE = YRDOY                ! Record forcing date.

!         Ready for next stage
          STGDOY(ISTAGE) = YRDOY

          ISTAGE = 5
          SUMDTTGRO= SUMDTT              
          SUMTMAXGRO= SUMTMAX
          SUMSRADGRO= SUMSRAD
          SUMPARGRO = SUMPAR
          TBASE  = 2.00                
          SUMDTT = DTT                 
          SUMTMAX= TMAX
          SUMSRAD= SRAD
          SUMPAR = PAR 

!-----------------------------------------------------------------
      CASE (5)         
          !
          !  Stage 5 >> Forcing to Open Heart 
          !
          IF (SUMDTT .LT. (P5)) THEN
             RETURN                       
          ENDIF

!         Ready for next stage
          STGDOY(ISTAGE) = YRDOY
          EDATE5 = YRDOY

          ISTAGE = 6
          SUMDTTGRO= SUMDTT             
          SUMTMAXGRO= SUMTMAX
          SUMSRADGRO= SUMSRAD
          SUMPARGRO= SUMPAR
          TBASE  = 2.0                 
          SUMDTT =  DTT                      
          SUMTMAX= TMAX
          SUMSRAD= SRAD
          SUMPAR = PAR 
!-----------------------------------------------------------------          
        CASE (6)            
          !
          !  Stage 6 >> Open Heart to EarlyAnthesis 
          !
          IF (SUMDTT .LT. P6) THEN        
             RETURN                      
          ENDIF

!         Ready for next stage
          STGDOY(ISTAGE) = YRDOY
          EDATE6 = YRDOY

          ISTAGE = 7                    
          SUMDTTGRO= SUMDTT               
          SUMTMAXGRO= SUMTMAX
          SUMSRADGRO= SUMSRAD
          SUMPARGRO= SUMPAR
          TBASE  = 2.0                 
          SUMDTT = DTT                 
          SUMTMAX= TMAX
          SUMSRAD= SRAD
          SUMPAR = PAR 

     
!-----------------------------------------------------------------
        CASE (7)             
          !
          !  Stage 7 >> Early Anthesis to Last Anthesis 
          !
          IF (SUMDTT .LT. P7) THEN        
             RETURN                       
          ENDIF

!         Ready for next stage
          STGDOY(ISTAGE) = YRDOY
          EDATE7 = YRDOY

          ISTAGE = 8                    
          SUMDTTGRO= SUMDTT               
          SUMTMAXGRO= SUMTMAX
          SUMSRADGRO= SUMSRAD
          SUMPARGRO= SUMPAR
          TBASE  = 2.0                 
          SUMDTT = DTT                  
          SUMTMAX= TMAX
          SUMSRAD= SRAD
          SUMPAR = PAR

!-----------------------------------------------------------------
        CASE (8)               
          !
          ! Stage 8 Last Anthesis to Physiological maturity
          !
          !  
 
          IF (SUMDTT .LT. P8) THEN 
             RETURN                        
          ENDIF
          
          STGDOY(ISTAGE) = YRDOY
          PMDATE = YRDOY                   ! Physiological maturity date PMDATE = YRDOY

!         Ready for next stage

          ISTAGE = 9                  
          SUMDTTGRO= SUMDTT               
          SUMTMAXGRO= SUMTMAX
          SUMSRADGRO= SUMSRAD
          SUMPARGRO= SUMPAR
          TBASE  = 2.0
          SUMDTT = DTT
          SUMTMAX= TMAX
          SUMSRAD= SRAD
          SUMPAR = PAR 

!-----------------------------------------------------------------
        CASE (9)
          !
          ! Stage 6 >> Physiological maturity to Harvest
          !
          IF (SUMDTT .LT. (G1)) THEN
             RETURN
          ENDIF


          
          STGDOY(ISTAGE) = YRDOY
          FHDATE = YRDOY                  ! harvest date FHDATE = YRDOY 
          MDATE  = YRDOY                  ! Set MDATE to stop model

!         Ready for next stage

           ISTAGE = 10
      
          SUMDTTGRO= SUMDTT               
          SUMTMAXGRO= SUMTMAX
          SUMSRADGRO= SUMSRAD
          SUMPARGRO= SUMPAR
          TBASE  = 2.0
          SUMDTT = DTT
          SUMTMAX= TMAX
          SUMSRAD= SRAD
          SUMPAR = PAR 

!-----------------------------------------------------------------


!-----------------------------------------------------------------
        CASE (10)             
          !
          ! Stage 10 Harvest
          !
          
          IF (SUMDTT .LT. (G1)) THEN    
             RETURN
          ENDIF

          
          STGDOY(ISTAGE) = YRDOY
          FHDATE = YRDOY                  ! harvest date FHDATE = YRDOY 
          MDATE  = YRDOY                  ! Set MDATE to stop model


      END SELECT
!-----------------------------------------------------------------

  !    IF (ISTAGE .NE. 6) THEN
  !       CALL PHASEI (ISWWAT,ISWNIT)
  !       RETURN
  !    ENDIF

!=================================================================
      END SELECT
!=================================================================


      RETURN

C-----------------------------------------------------------------------
C     Format Strings
C-----------------------------------------------------------------------

1399  FORMAT ('Crop failure because of lack of root initiation',
     1       ' within 140 days of planting')
2380  FORMAT ('Crop failure - Growth program terminated')
3600  FORMAT (1X,'Crop failure because of lack of germination ',
     1           'within 40 days of sowing')

      END SUBROUTINE Aloha2_PHENOL
!=================================================================