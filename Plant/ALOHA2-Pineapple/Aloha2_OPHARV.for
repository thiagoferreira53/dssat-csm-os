!=======================================================================
!  Aloha2_OPHARV, Subroutine C.H. Porter
!  Generates data for use by OPSUM and OVERVIEW for Pineapple.
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  06/24/2017 CHP Written, based on MZ_OPHARV
!  09/05/2020 JVJ Modifying the PlantStres% NSTAGES names     
!=======================================================================

      SUBROUTINE Aloha2_OPHARV(CONTROL, ISWITCH,
     &   AGEFAC, BIOMAS, CNAM, CRWNWT, EYEWT, FBIOM,      
     &   FRTWT, FRUITS, GPSM, GPP, HARVFRAC, ISDATE,      
     &   ISTAGE, LAI, LN, MDATE, NSTRES, PLTPOP, PMDATE,  
     &   PSTRES1, PSTRES2, STGDOY, STOVER, SWFAC,         
     &   TURFAC, VWATM, WTINITIAL, WTNCAN, WTNGRN,        
     &   WTNUP, YIELD, YRDOY, YRPLT, EDATE12, EDATE13,
     &   EDATE1, EDATE2, EDATE3, EDATE5, EDATE6, EDATE7,
     &   BIOMAS1, LAI1, BIOMAS13, LAI13, BIOMAS2, LAI2,
     &   BIOMAS3, LAI3, BIOMAS4, LAI4, LN2, LN3, LN4, HIFact, MAXLAI) 

!-----------------------------------------------------------------------
      USE Aloha2_mod
      IMPLICIT NONE
      SAVE

      CHARACTER*1  IDETO, IDETS, IPLTI, RNMODE
      CHARACTER*2  CROP
      CHARACTER*6  SECTION
      CHARACTER*6, PARAMETER :: ERRKEY = 'OPHARV'
      CHARACTER*12 FILEA
      CHARACTER*30 FILEIO
	    CHARACTER*80 PATHEX


      INTEGER DMAT
      INTEGER DNR7, DYNAMIC, ERRNUM, FOUND
      INTEGER IMAT, IFORC, IHARV, DFR1, DFORC, DHARV, HDAP
      INTEGER DROOT, IROOT, RIDAP, EDATE12
      INTEGER DFNL, IFNL, LEDAP, EDATE13
      INTEGER DLC1, ILC1, L1DAP, EDATE1
      INTEGER DLC2, ILC2, L2DAP, EDATE2
      INTEGER DLC3, ILC3, L3DAP, EDATE3
      INTEGER DROH, IROH, OHDAF, EDATE5
      INTEGER DREA, IREA, EADAF, EDATE6
      INTEGER DRLA, IRLA, LADAF, EDATE7
      INTEGER ISDATE, ISENS, LINC, LNUM, LUNIO, MDATE, ISTAGE, RUN
      INTEGER TIMDIF, TRTNUM, YRNR1, YRNR2, YRNR3
      INTEGER YRDOY, YREMRG, YRNR5, YRSIM, YRPLT
      INTEGER DNR3, PMDATE
      INTEGER TRT_ROT
      INTEGER STGDOY(20)
      
      REAL AGEFAC, BWAH
      REAL CRWNWT, EYEWT, FBIOM, FBTONS, FRTWT, FRUITS
      REAL GPP, GPSM, HI, RIBIO, BIOMAS1, RILAI, LAI1
      REAL LEBIO, BIOMAS13, LELAI, LAI13
      REAL L1BIO, BIOMAS2, L1LAI, LAI2, L1LN, LN2
      REAL L2BIO, BIOMAS3, L2LAI, LAI3, L2LN, LN3
      REAL L3BIO, BIOMAS4, L3LAI, LAI4, L3LN, LN4, HIFact
      REAL MAXLAI, NSTRES, PBIOMS, PEYEWT, PSDWT, PLTPOP 
      REAL Pstres1, Pstres2   
      REAL SDRATE
      REAL SDWT, SDWTAH, SDWTAM, STOVER   !, StovSenes  
      REAL SWFAC, BIOMAS, Biomass_kg_ha, TURFAC
      REAL WTINITIAL, WTNCAN, WTNGRN, WTNUP, LAI, LN
      REAL YIELD, YIELDB, YieldFresh
      REAL VWATM, CNAM, BWAM, HWAH, StovSenes

      REAL, DIMENSION(2) :: HARVFRAC

!     Arrays which contain data for printing in SUMMARY.OUT file
      INTEGER, PARAMETER :: SUMNUM = 18
      CHARACTER*4, DIMENSION(SUMNUM) :: LABEL
      REAL, DIMENSION(SUMNUM) :: VALUE

!     Arrays which contain Simulated and Measured data for printing
!       in OVERVIEW.OUT and EVALUATE.OUT files (OPVIEW subroutine)
      INTEGER ACOUNT
      CHARACTER*6, DIMENSION(EvaluateNum) :: OLAB, OLAP !OLAP in dap
      CHARACTER*12 X(EvaluateNum)
      CHARACTER*8 Simulated(EvaluateNum), Measured(EvaluateNum)
      CHARACTER*50 DESCRIP(EvaluateNum)

      TYPE (ControlType) CONTROL
      TYPE (SwitchType)  ISWITCH
      TYPE (ResidueType) SENESCE

!     Variables added for environmental and stress factors output
      TYPE (PlStresType) PlantStres

      DYNAMIC = CONTROL % DYNAMIC

!-----------------------------------------------------------------------
      ACOUNT = 36  !Number of FILEA headings.
!     Headings in FILEA for Measured data
      DATA OLAB /
     &    'RIDAP',   !1Pineapple root initiation (was not)
     &    'LEDAP',   !2Pineapple first new leaf (was not)
     &    'L1DAP',   !3Pineapple leaf cycle 1 (was not)
     &    'L2DAP',   !4Pineapple leaf cycle 2 (was not)
     &    'L3DAP',   !5Pineapple leaf cycle 3 (was not)
     &    'FDAT',    !61  forcing date
     &    'OHDAF',   !7Pineapple open heart (was not)
     &    'EADAF',   !8Pineapple early anthesys (was not)
     &    'LADAF',   !9Pineapple last anthesys (was not)
     &    'MDAT',    !10 2  maturity date (was 4) 
     &    'HDAT',    !11 3  harvest date (was 2)
     &    'RIBIO',   !12Pineapple biomass at root initiation (was not)
     &    'RILAI',   !13Pineapple leaf area index at root initiation (was not)
     &    'LEBIO',   !14Pineapple biomass at first new leaf (was not)
     &    'LELAI',   !15Pineapple leaf area index at first new leaf (was not)
     &    'L1BIO',   !16Pineapple biomass at Leaf Cycle 1 (was not)
     &    'L1LAI',   !17Pineapple leaf area index at Leaf Cycle 1 (was not)
     &    'L1LN',    !18Pineapple leaf number al Leaf Cycle 1 (eas not)
     &    'L2BIO',   !19Pineapple biomass at Leaf Cycle 2 (was not)
     &    'L2LAI',   !20Pineapple leaf area index at Leaf Cycle 2 (was not)
     &    'L2LN',    !21Pineapple leaf number al Leaf Cycle 2 (eas not)
     &    'L3BIO',   !22Pineapple biomass at Leaf Cycle 3 (was not)
     &    'L3LAI',   !23Pineapple leaf area index at Leaf Cycle 3 (was not)
     &    'L3LN',    !24Pineapple leaf number al Leaf Cycle 3 (was not)
     &    'FWAH',    !25 4  yield fresh weight, t/ha (was 5) 
     &    'YDWAH',   !26 5  yield dry weight, t/ha (new) 
     &    'BADMF',   !27 6  biomass at forcing, t/ha (was 10)
     &    'BADMH',   !28 7  biomass at harvest, t/ha (was 11)
     &    'VWATM',   !29 8  veg dry matter @ mat, t/ha (was 6)
     &    'LAIX',    !30 9  max LAI (was 12) 
     &    'L#SM',    !31 10 Leaf # (was 19)
     &    'HIAM',    !32 11 HI (was 13)
     &    'E#AM',    !33 12 eye number/m2 (was 7)
     &    'E#UM',    !34 13 eye number/fruit (was 9)
     &    'EWUM',    !35 14 eye unit weight (was 8)
     &    'CNAM',    !36 15 tops N, kg/ha (was 16)


   
!    &    'PDFT',    !not used (was 3)
!    &    'THAM',    !not used (was 14)
!    &    'FRNAM',   !not used (was 15)
!    &    'FRN%M',   !not used (was 18)
!    &    'VNAM',    !not used - same as CNAM (was 17)

     &    4*'    '/ 
!                                                                                       Old    New 
! ORIGINAL Aloha2 MODEL OUTPUT:                                           Simul   Meas   FileA  FileA
! 400 FORMAT (6X, 'FORCING DATE (DAP)             ',2X,I6,7X,I6,/,  !  1 DNR1,   DFLR   FDAT - OK
!    &        6X, 'HARVEST DATE (DAP)             ',2X,I6,7X,I6,/,  ! 2? DNR3,   PMAT   HDAT - OK
!    &        6X, 'PHYSIO. MATURITY DATE          ',2X,I6,7X,I6,/,  ! 4? DNR7,   DMAT   MDAT - OK
!    &        6X, 'FRESH FRUIT YIELD (T/HA)       ',F8.2, 7X,A6,/,  !  5 FYIELD, XGWT   FWAH - OK
!    &        6X, 'EYE WEIGHT (G)                 ',F8.3, 7X,A6,/,  !  8 EYEWT,  XGWU   EWUM - OK
!    &        6X, 'EYE PER SQ METER               ',F8.0, 7X,A6,/,  !  7 GPSM,   XNOGR  E#AM - OK
!    &        6X, 'EYE PER FRUIT                  ',F8.2, 7X,A6,/,  !  9 GPP,    XNOGU  E#UM - OK
!    &        6X, 'MAX. LAI                       ',F8.2, 7X,A6,/,  ! 12 MAXLAI, XLAM   LAIX - OK
!    &        6X, 'DRY BIOM @ FORC. (T/HA)        ',F8.2, 7X,A6,/,  ! 10 FBTONS, XCWT   CWAM - BADMF, aerial biomass at forcing, t/ha
!    &        6X, 'DRY BIOMASS (T/HA)             ',F8.2, 7X,A6,/,  ! 11 PBIOMS, XSWT   BWAM - BADMH, aerial biomass at harvest, t/ha
!    &        6X, 'DRY VEG. WT. (T/HA)            ',F8.2, 7X,A6,/,  !  6 PVEGWT, XPDW   PWAM - VWATM, veg wt @ maturity, t/ha 
!    &        6X, 'FRUIT N%                       ',F8.2, 7X,A6,/,  ! 18 XGNP,   XNPS   GN%M - FRN%M, fruit N concentration at maturity, %
!    &        6X, 'TOT N UPTAKE (KG N/HA)         ',F8.1, 7X,A6,/,  ! 16 TOTNUP, XNTP   CNAM - OK
!    &        6X, 'VEGETATIVE N UPTAKE            ',F8.1, 7X,A6,/,  ! 17 APTNUP, XNST   SNAM - VNAM, veg N at maturity, kg/ha
!    &        6X, 'FRUIT N UPTAKE                 ',F8.1, 7X,A6)    ! 15 GNUP,   XNGR   GNAM - FRNAM, fruit N at maturity, kg/ha

!***********************************************************************
!***********************************************************************
!     RUN INITIALIZATION
!***********************************************************************
      IF (DYNAMIC .EQ. RUNINIT) THEN
!-----------------------------------------------------------------------
      CROP    = CONTROL % CROP
      YRSIM   = CONTROL % YRSIM
      RNMODE  = CONTROL % RNMODE
      FILEIO  = CONTROL % FILEIO
      RUN     = CONTROL % RUN
      IPLTI   = ISWITCH % IPLTI

!       Read FILEIO
        CALL GETLUN('FILEIO', LUNIO)
        OPEN (LUNIO, FILE = FILEIO, STATUS = 'OLD', IOSTAT=ERRNUM)
        IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,0)

        READ (LUNIO,'(55X,I5)', IOSTAT=ERRNUM) ISENS; LNUM = 1
        IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)

        READ (LUNIO,'(3(/),15X,A12,1X,A80)', IOSTAT=ERRNUM) FILEA,
     &        PATHEX
        LNUM = LNUM + 4  
        IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)
  
        SECTION = '*TREAT'
        CALL FIND(LUNIO, SECTION, LINC, FOUND) ; LNUM = LNUM + LINC
        IF (FOUND .EQ. 0) THEN
          CALL ERROR(SECTION, 42, FILEIO, LNUM)
        ELSE
          READ(LUNIO, '(I3)', IOSTAT=ERRNUM) TRTNUM ; LNUM = LNUM + 1
          IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)
        ENDIF

        CLOSE (LUNIO)

!     Assign descriptions to Measured and Simulated data 
!         from DATA.CDE.
      CALL GETDESC(ACOUNT, OLAB, DESCRIP)
      OLAP = OLAB

      Pstres1 = 1.0
      Pstres2 = 1.0

      IDETO = ISWITCH % IDETO

      Biomass_kg_ha = BIOMAS * 10. !Convert from g/m2 to kg/ha
      
      CALL OPVIEW(CONTROL, 
     &    Biomass_kg_ha, ACOUNT, DESCRIP, IDETO, LN, 
     &    Measured, PlantStres, Simulated, STGDOY, 
     &    STNAME, WTNCAN, LAI, NINT(YIELD), YRPLT, ISTAGE)

!***********************************************************************
!***********************************************************************
!     SEASONAL INITIALIZATION
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
      CROP    = CONTROL % CROP
      YRSIM   = CONTROL % YRSIM
      RNMODE  = CONTROL % RNMODE
      FILEIO  = CONTROL % FILEIO
      RUN     = CONTROL % RUN

      Simulated = ' '
      Measured  = ' '

!     Establish # and names of stages for environmental & stress summary
      PlantStres % ACTIVE = .FALSE.
      PlantStres % NSTAGES = 5

      PlantStres % StageName(0) = 'Planting to Harvest  '   
      PlantStres % StageName(1) = 'PLT - LC1 '         
      PlantStres % StageName(2) = 'LC1 - LC2 '              
      PlantStres % StageName(3) = 'LC2 - Forcing'            
      PlantStres % StageName(4) = 'Forcing - LastAnthes'    
      PlantStres % StageName(5) = 'LastAnthes - PhMaturity' 
      

      Biomass_kg_ha = BIOMAS * 10. !Convert from g/m2 to kg/ha

      CALL OPVIEW(CONTROL, 
     &    Biomass_kg_ha, ACOUNT, DESCRIP, IDETO, LN, 
     &    Measured, PlantStres, Simulated, STGDOY, 
     &    STNAME, WTNCAN, LAI, NINT(YIELD), YRPLT, ISTAGE)

      MAXLAI = 0.0

!***********************************************************************
!***********************************************************************
!     DAILY OUTPUT
!***********************************************************************
      ELSE IF (DYNAMIC .EQ. OUTPUT) THEN
!-----------------------------------------------------------------------
      MAXLAI = AMAX1 (MAXLAI,LAI)      ! Maximum LAI season

      PlantStres % W_grow = TURFAC 
      PlantStres % W_phot = SWFAC  
      PlantStres % N_grow = AGEFAC 
      PlantStres % N_phot = NSTRES 
      PlantStres % P_grow = PSTRES2
      PlantStres % P_phot = PSTRES1

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

!     01/08/2024 TF - Added CASEs to include new stages
      PlantStres % ACTIVE = .FALSE.
      SELECT CASE(ISTAGE)
      CASE(12,13,1)         
        PlantStres % ACTIVE(1) = .TRUE.
      CASE(2)                 
        PlantStres % ACTIVE(2) = .TRUE.
      CASE(3,4)                 
        PlantStres % ACTIVE(3) = .TRUE.
      CASE(5,6,7)                 
        PlantStres % ACTIVE(4) = .TRUE.
      CASE(8)                 
        PlantStres % ACTIVE(5) = .TRUE.
      
      END SELECT


      YRDOY = CONTROL % YRDOY
      IF (YRDOY >= YRPLT) THEN
        IF (MDATE < 0 .OR.
     &     (MDATE > 0 .AND. YRDOY < MDATE)) THEN
          PlantStres % ACTIVE(0) = .TRUE.
        ELSE
          PlantStres % Active(0) = .FALSE.
        ENDIF
      ENDIF

      Biomass_kg_ha = BIOMAS * 10. !Convert from g/m2 to kg/ha

!     Send data to Overview.out data on days where stages occur
      CALL OPVIEW(CONTROL, 
     &    Biomass_kg_ha, ACOUNT, DESCRIP, IDETO, LN, 
     &    Measured, PlantStres, Simulated, STGDOY, 
     &    STNAME, WTNCAN, LAI, NINT(YIELD), YRPLT, ISTAGE)

!***********************************************************************
!***********************************************************************
!     Seasonal Output 
!***********************************************************************
      ELSE IF (DYNAMIC .EQ. SEASEND) THEN
!-----------------------------------------------------------------------
!     Transfer dates for model stages.
      YRNR1  = ISDATE
      YRNR2  = STGDOY(2)
      YRNR3  = STGDOY(3)
      YRNR5  = STGDOY(5)
!      YRNR7  = MDATE
!      WTNUP  = TOTNUP/10.0
!      WTNCAN = TOTNUP/10.0
!      WTNSD  = GNUP  /10.0

!     For now the senesced leaf mass is not tracked, so this will be zero
      StovSenes = SENESCE % ResWt(0)
C-----------------------------------------------------------------------
C     Adjust dates since Pineapple grows after harvest
C-----------------------------------------------------------------------

!      MDATE = FHDATE

!-----------------------------------------------------------------------
!     Actual yield harvested (default is 100 %)
!     YIELD variable is in kg/ha
!     FRUITS = Fruits by squere meter
!-----------------------------------------------------------------------
      IF (PLTPOP .GE. 0.0) THEN
         IF (GPP .GT. 0.0) THEN           ! GPP = eyes per fruit
            EYEWT = (FRTWT+CRWNWT)/GPP    ! Eye weight (g/eye) in MD-2 or fresh fruit pineapple crown is included
            PEYEWT = EYEWT * 1000.0       ! Eye weight (mg/eye)
         ENDIF
         GPSM   = GPP*FRUITS              ! # eyes/m2
         STOVER = BIOMAS*10.0-YIELD       ! Total plant weight except fruit (g/m2)
         !YIELDFresh  = YIELD / Species % FDMC  ! Fresh fruit yield (kg/ha)
         YIELDFresh  = ((FRTWT*10.0*FRUITS) + (CRWNWT*10.0*FRUITS)) / Species % FDMC  ! Fresh fruit yield (kg/ha)
         YIELDB = YIELDFresh / 0.8914          ! Fresh fruit yield (lb/acre)
      ENDIF

      SDWT   = YIELD  / 10.0      !g/m2
      SDWTAM = YIELD / 10.0       !g/m2
      SDWTAH = SDWT * HARVFRAC(1) !g/m2

      SDRATE = WTINITIAL
     

!     CHP 10/17/2017
!     Oddly, PSDWT = SDWT/GPSM is exactly equal to EYEWT = FRTWT/GPP
!       Actually, not so odd:
!       PSDWT = SDWT/GPSM = FRTWT*FRUITS/GPSM = (FRTWT*FRUITS)/(GPP*FRUITS) = FRTWT/GPP = EYEWT
      PSDWT = 0.0
      IF (GPSM .GT. 0.0 .AND. SDWT  .GE. 0.0) THEN
         PSDWT = SDWT/GPSM            !Unit weight of eye g/unit = (g/m2) / (#/m2)
      ENDIF
      
      IF (BIOMAS .GT. 0.0 .AND. YIELD .GE. 0.0 
     & .AND. FBIOM .GT. 0.0) THEN
          !HI = (YIELD + (CRWNWT*10*FRUITS))/(FBIOM*10)
          !HI = ((FRTWT*10.0*FRUITS) + (CRWNWT*10.0*FRUITS))/(FBIOM*10)
          !HI = ((YIELD/1000)/FBIOM*10)*10
          HI = (YIELD/FBIOM)/10
       ELSE
         HI = 0.0
      ENDIF

!-----------------------------------------------------------------------
!     Actual yield harvested (default is 100 %)
      HWAH = YIELD * HARVFRAC(1)

!     Actual byproduct harvested (default is 0 %)
!     Byproduct not harvested is incorporated
      BWAM = BIOMAS*10. - YIELD 
      BWAH = (BWAM + StovSenes) * HARVFRAC(2) 
!-----------------------------------------------------------------------

      RIBIO = BIOMAS13*10
      RILAI = LAI13
      
      
      LEBIO = BIOMAS1*10
      LELAI = LAI1

      L1BIO = BIOMAS2*10
      L1LAI = LAI2
      L1LN  = LN2
      
      L2BIO = BIOMAS3*10
      L2LAI = LAI3
      L2LN  = LN3

      L3BIO = BIOMAS4*10
      L3LAI = LAI4
      L3LN  = LN4
!-----------------------------------------------------------------------
  
      IF ((INDEX('YE',IDETO) > 0 .OR. INDEX('IAEBCGDT',RNMODE) .GT. 0) 
     &  .OR. (INDEX('AY',IDETS) .GT. 0 .AND. CROP .NE. 'FA')) THEN
         IF (INDEX('FQ',RNMODE) > 0) THEN
           TRT_ROT = CONTROL % ROTNUM
         ELSE
           TRT_ROT = TRTNUM
         ENDIF
         !CALL READA (FILEA, PATHEX,OLAB, TRT_ROT, YRSIM, X)
         CALL READA_Y4K(FILEA, PATHEX,OLAB, TRT_ROT, YRSIM, X)

!-----------------------------------------------------------------------
      FBTONS = FBIOM*0.01             !Biomass at forcing, t/ha
      PBIOMS = (BIOMAS*10.0)/1000.0   !Biomass at maturity, t/ha
!     FYIELD = YIELD/1000.0
!     PVEGWT = STOVER/1000.0
      Biomass_kg_ha = BIOMAS * 10.    !Convert from g/m2 to kg/ha

!     Change observed (FileA) dates to DAP
!       Forcing date to DAP
        CALL READA_Dates(X(6), YRSIM, IFORC)
        IF (IFORC .GT. 0 .AND. IPLTI .EQ. 'R' .AND. ISENS .EQ. 0) THEN
          DFORC = TIMDIF(YRPLT,IFORC)
        ELSE
          DFORC  = -99
        ENDIF
        OLAP(6) = 'FDAP  '
        CALL GetDesc(1,OLAP(6), DESCRIP(6))

!       Maturity date to DAP
        CALL READA_Dates(X(10), YRSIM, IMAT)
        IF (IMAT .GT. 0 .AND. IPLTI .EQ. 'R' .AND. ISENS .EQ. 0) THEN
          DMAT = TIMDIF(YRPLT,IMAT)
        ELSE
          DMAT  = -99
        ENDIF
        OLAP(10) = 'MDAP  '
        CALL GetDesc(1,OLAP(10), DESCRIP(10))

!       Harvest date to DAP
        CALL READA_Dates(X(11), YRSIM, IHARV)
        IF (IHARV .GT. 0 .AND. IPLTI .EQ. 'R' .AND. ISENS .EQ. 0) THEN
          DHARV = TIMDIF(YRPLT,IHARV)
        ELSE
          DHARV  = -99
        ENDIF
        OLAP(11) = 'HDAP  '
        CALL GetDesc(1,OLAP(11), DESCRIP(11))

!       Root Initiation date to DAP PART1
        CALL READA_Dates(X(1), RIDAP, IROOT)   !IROOT Rooting
        IF (IROOT .GT. 0 .AND. IPLTI .EQ. 'R' .AND. ISENS .EQ. 0) THEN
          DROOT = IROOT
          
        ELSE
          DROOT  = -99
        ENDIF

        CALL GetDesc(1,OLAP(1), DESCRIP(1))
        OLAP(1) = 'RIDAP '

!       First new leaf date to DAP PART1

       CALL READA_Dates(X(2), LEDAP, IFNL)  !IFNL first new leaf
        IF (IFNL .GT. 0 .AND. IPLTI .EQ. 'R' .AND. ISENS .EQ. 0) THEN
          DFNL = IFNL
          
        ELSE
          DFNL  = -99
        ENDIF
        
        CALL GetDesc(1,OLAP(2), DESCRIP(2))
          OLAP(2) = 'LEDAP '


!       Leaf cicle 1 date to DAP PART1

       CALL READA_Dates(X(3), L1DAP, ILC1)  !IFNL first new leaf
        IF (ILC1 .GT. 0 .AND. IPLTI .EQ. 'R' .AND. ISENS .EQ. 0) THEN
          DLC1 = ILC1
          
        ELSE
          DLC1  = -99
        ENDIF
        
        CALL GetDesc(1,OLAP(3), DESCRIP(3))
          OLAP(3) = 'L1DAP '

!       Leaf cicle 2 date to DAP PART1

       CALL READA_Dates(X(4), L2DAP, ILC2)  !IFNL first new leaf
        IF (ILC2 .GT. 0 .AND. IPLTI .EQ. 'R' .AND. ISENS .EQ. 0) THEN
          DLC2 = ILC2
          
        ELSE
          DLC2  = -99
        ENDIF
        
        CALL GetDesc(1,OLAP(4), DESCRIP(4))
          OLAP(4) = 'L2DAP '

!       Leaf cicle 3 date to DAP PART1

       CALL READA_Dates(X(5), L3DAP, ILC3)  !IFNL first new leaf
        IF (ILC3 .GT. 0 .AND. IPLTI .EQ. 'R' .AND. ISENS .EQ. 0) THEN
          DLC3 = ILC3
          
        ELSE
          DLC3  = -99
        ENDIF
        
        CALL GetDesc(1,OLAP(5), DESCRIP(5))
        OLAP(5) = 'L3DAP '


!       Reproductive Open Heart date to DAP PART1

       CALL READA_Dates(X(7), OHDAF, IROH)  !IFNL first new leaf
        IF (IROH .GT. 0 .AND. IPLTI .EQ. 'R' .AND. ISENS .EQ. 0) THEN
          DROH = IROH
          
        ELSE
          DROH  = -99
        ENDIF
        
        CALL GetDesc(1,OLAP(7), DESCRIP(7))
        OLAP(7) = 'OHDAF '


!       Reproductive early anthesys date to DAP PART1

       CALL READA_Dates(X(8), EADAF, IREA)  !IFNL first new leaf
        IF (IREA .GT. 0 .AND. IPLTI .EQ. 'R' .AND. ISENS .EQ. 0) THEN
          DREA = IREA
          
        ELSE
          DREA  = -99
        ENDIF
        
        CALL GetDesc(1,OLAP(8), DESCRIP(8))
        OLAP(8) = 'EADAF '


!       Reproductive last anthesys date to DAP PART1

       CALL READA_Dates(X(9), LADAF, IRLA)  !IFNL first new leaf
        IF (IRLA .GT. 0 .AND. IPLTI .EQ. 'R' .AND. ISENS .EQ. 0) THEN
          DRLA = IRLA
          
        ELSE
          DRLA  = -99
        ENDIF
        
        CALL GetDesc(1,OLAP(9), DESCRIP(9))
        OLAP(9) = 'LADAF '

!       Root Initiation date to DAP PART2
        IF (YRPLT .GT. 0) THEN
          RIDAP = TIMDIF (YRPLT,EDATE12)
          

          IF (RIDAP .LE. 0) THEN
            RIDAP = -99
          ENDIF
        ELSE
          RIDAP = -99
        ENDIF

!       First new leaf date to DAP PART2

       IF (YRPLT .GT. 0) THEN
          LEDAP = TIMDIF (YRPLT,EDATE13)
          IF (LEDAP .LE. 0) THEN
            LEDAP = -99
          ENDIF
        ELSE
          LEDAP = -99
        ENDIF

!       Leaf cycle 1 date to DAP PART2

       IF (EDATE1 .GT. 0 .AND. EDATE1 .LE. ISDATE) THEN
          L1DAP = TIMDIF (YRPLT,EDATE1)
          IF (L1DAP .LE. 0) THEN
            L1DAP = -99
          ENDIF
        ELSE
          L1DAP = -99
        ENDIF

!       Leaf cycle 2 date to DAP PART2

       IF (EDATE2 .GT. 0 .AND. EDATE2 .LE. ISDATE) THEN
          L2DAP = TIMDIF (YRPLT,EDATE2)
          IF (L2DAP .LE. 0) THEN
            L2DAP = -99
          ENDIF
        ELSE
          L2DAP = -99
        ENDIF

!       Leaf cycle 3 date to DAP PART2

       IF (EDATE3 .GT. 0 .AND. EDATE3 .LE. ISDATE) THEN
          L3DAP = TIMDIF (YRPLT,EDATE3)
          IF (L3DAP .LE. 0) THEN
            L3DAP = -99
          ENDIF
        ELSE
          L3DAP = -99
        ENDIF

!       Change simulated dates to DAP
!       isdate is forcing date
      
        IF (YRPLT .GT. 0) THEN
          DFR1 = TIMDIF (YRPLT,ISDATE)
          IF (DFR1 .LE. 0) THEN
            DFR1 = -99
          ENDIF
        ELSE
          DFR1 = -99
        ENDIF
  
        IF (YRPLT .GT. 0) THEN
          DNR3 = TIMDIF (YRPLT,MDATE)
          IF (DNR3 .LE. 0) THEN
            DNR3 = -99
        ELSE
          DNR3 = -99
        ENDIF
      
        IF (YRPLT .GT. 0) THEN
          HDAP = TIMDIF (YRPLT,YRDOY)
          IF (HDAP .LE. 0) THEN
            HDAP = -99
          ENDIF
        ELSE
          HDAP = -99
        ENDIF
      
        IF (YRPLT .GT. 0) THEN
          DNR7 = TIMDIF (YRPLT,PMDATE)
          IF (DNR7 .LE. 0)  THEN
            DNR7 = -99
          ENDIF
        ELSE
          DNR7 = -99
        ENDIF
      ENDIF
      
      !       Reproductive Open Heart date to DAP PART2

       IF (YRPLT .GT. 0) THEN
          OHDAF = (TIMDIF (YRPLT,EDATE5)) - DFR1
          IF (OHDAF .LE. 0) THEN
            OHDAF = -99
          ENDIF
        ELSE
          OHDAF = -99
        ENDIF

!       Reproductive Early Anthesys date to DAP PART2

       IF (YRPLT .GT. 0) THEN
          EADAF = (TIMDIF (YRPLT,EDATE6)) - DFR1
          IF (EADAF .LE. 0) THEN
            EADAF = -99
          ENDIF
        ELSE
          EADAF = -99
        ENDIF

!       Reproductive Last Anthesys date to DAP PART2

       IF (YRPLT .GT. 0) THEN
          LADAF = (TIMDIF (YRPLT,EDATE7)) - DFR1
          IF (LADAF .LE. 0) THEN
            LADAF = -99
          ENDIF
        ELSE
          LADAF = -99
        ENDIF

      WRITE(Simulated(1),'(I8)') RIDAP;         WRITE(Measured(1),'(A8)') TRIM(X(1))    !RIDAP
      WRITE(Simulated(2),'(I8)') LEDAP;         WRITE(Measured(2),'(A8)') TRIM(X(2))    !LEDAP
      WRITE(Simulated(3),'(I8)') L1DAP;         WRITE(Measured(3),'(A8)') TRIM(X(3))   !DLC1    !L1DAP
      WRITE(Simulated(4),'(I8)') L2DAP;         WRITE(Measured(4),'(A8)') TRIM(X(4))    !DLC2    !L2DAP
      WRITE(Simulated(5),'(I8)') L3DAP;         WRITE(Measured(5),'(A8)') TRIM(X(5))    !DLC3    !L3DAP
      WRITE(Simulated(6),'(I8)') DFR1;          WRITE(Measured(6),'(I8)') DFORC   !FDAT
      WRITE(Simulated(7),'(I8)') OHDAF;         WRITE(Measured(7),'(A8)') TRIM(X(7))    !DROH    !OHDAF
      WRITE(Simulated(8),'(I8)') EADAF;         WRITE(Measured(8),'(A8)') TRIM(X(8))    !DREA    !EADAF
      WRITE(Simulated(9),'(I8)') LADAF;         WRITE(Measured(9),'(A8)') TRIM(X(9))    !DRLA    !LADAF
      WRITE(Simulated(10),'(I8)') DNR7;         WRITE(Measured(10),'(I8)') DMAT   !MDAT  
      WRITE(Simulated(11),'(I8)') HDAP;         WRITE(Measured(11),'(I8)') DHARV  !HDAT   
      WRITE(Simulated(12),'(F8.1)') RIBIO;      WRITE(Measured(12),'(A8)') TRIM(X(12))  !RIBIO
      WRITE(Simulated(13),'(F8.2)') RILAI;      WRITE(Measured(13),'(A8)') TRIM(X(13))  !RILAI
      WRITE(Simulated(14),'(F8.1)') LEBIO;      WRITE(Measured(14),'(A8)') TRIM(X(14))  !LEBIO
      WRITE(Simulated(15),'(F8.2)') LELAI;      WRITE(Measured(15),'(A8)') TRIM(X(15))  !LELAI
      WRITE(Simulated(16),'(F8.1)') L1BIO;      WRITE(Measured(16),'(A8)') TRIM(X(16))  !L1BIO
      WRITE(Simulated(17),'(F8.2)') L1LAI;      WRITE(Measured(17),'(A8)') TRIM(X(17))  !L1LAI
      WRITE(Simulated(18),'(F8.1)') L1LN;       WRITE(Measured(18),'(A8)') TRIM(X(18))  !L1LN
      WRITE(Simulated(19),'(F8.1)') L2BIO;      WRITE(Measured(19),'(A8)') TRIM(X(19))  !L2BIO
      WRITE(Simulated(20),'(F8.2)') L2LAI;      WRITE(Measured(20),'(A8)') TRIM(X(20))  !L2LAI
      WRITE(Simulated(21),'(F8.1)') L2LN;       WRITE(Measured(21),'(A8)') TRIM(X(21))  !L2LN
      WRITE(Simulated(22),'(F8.1)') L3BIO;      WRITE(Measured(22),'(A8)') TRIM(X(22))  !L3BIO
      WRITE(Simulated(23),'(F8.2)') L3LAI;      WRITE(Measured(23),'(A8)') TRIM(X(23))  !L3LAI
      WRITE(Simulated(24),'(F8.1)') L3LN;       WRITE(Measured(24),'(A8)') TRIM(X(24))  !L3LN
      WRITE(Simulated(25),'(F8.2)') YIELDFresh/1000.                                     
                                                WRITE(Measured(25),'(A8)') TRIM(X(25))  !FWAH 
      WRITE(Simulated(26),'(F8.2)') YIELD/1000.                                     
                                                WRITE(Measured(26),'(A8)') TRIM(X(26))  !YDWAH 
      WRITE(Simulated(27),'(F8.1)') FBTONS                                      
                                                WRITE(Measured(27),'(A8)') TRIM(X(27))  !BADMF     
      WRITE(Simulated(28),'(F8.1)') PBIOMS                                        
                                                WRITE(Measured(28),'(A8)') TRIM(X(28)) !BADMH    
      WRITE(Simulated(29),'(F8.2)')VWATM;       WRITE(Measured(29),'(A8)') TRIM(X(29)) !VWATM  
      WRITE(Simulated(30),'(F8.2)') MAXLAI                                          
                                                WRITE(Measured(30),'(A8)') TRIM(X(30))  !LAIX     
      WRITE(Simulated(31),'(F8.1)') LN;         WRITE(Measured(31),'(A8)') TRIM(X(31))   !L#SM 
      WRITE(Simulated(32),'(F8.3)') HI;         WRITE(Measured(32),'(A8)') TRIM(X(32))   !HIAM
      WRITE(Simulated(33),'(I8)') NINT(GPSM)                                      
                                                WRITE(Measured(33),'(A8)') TRIM(X(33))   !E#AM 
      WRITE(Simulated(34),'(F8.1)') GPP;        WRITE(Measured(34),'(A8)') TRIM(X(34))   !E#UM     
      WRITE(Simulated(35),'(F8.3)') EYEWT                                          
                                                WRITE(Measured(35),'(A8)') TRIM(X(35))   !EWUM    
      WRITE(Simulated(36),'(F8.1)')CNAM;        WRITE(Measured(36),'(A8)') TRIM(X(36))   !CNAM
  
!     These aren't calculated - remove from Overview output
!     WRITE(Simulated(3),'(I8)') -99 ;  WRITE(Measured(3),'(I8)') -99     !PDFT - not used
!     WRITE(Simulated(14),'(I8)') -99 ; WRITE(Measured(14),'(I8)') -99    !THAM
!     WRITE(Simulated(15),'(F8.1)') FRNAM
!                                       WRITE(Measured(15),'(A8)') X(15)  !FRNAM
!     WRITE(Simulated(18),'(F8.1)')FRNpctM
!                                       WRITE(Measured(18),'(A8)') X(18)  !FRN%M
!     WRITE(Simulated(16),'(F8.1)')VNAM;WRITE(Measured(16),'(A8)')X(16)   !VNAM

      ENDIF

!-------------------------------------------------------------------
!     Send information to OPSUM to generate SUMMARY.OUT file
!-------------------------------------------------------------------
!     Store Summary.out labels and values in arrays to send to
!     OPSUM routines for printing.  Integers are temporarily 
!     saved as real numbers for placement in real array.
      LABEL(1)  = 'ADAT'; VALUE(1)  = FLOAT(ISDATE)
      LABEL(2)  = 'MDAT'; VALUE(2)  = FLOAT(MDATE)
      LABEL(3)  = 'DWAP'; VALUE(3)  = SDRATE
      LABEL(4)  = 'CWAM'; VALUE(4)  = BIOMAS*10.
      LABEL(5)  = 'HWAM'; VALUE(5)  = YIELD
      LABEL(6)  = 'HWAH'; VALUE(6)  = HWAH
!     BWAH multiplied by 10.0 in OPSUM - divide by 10. here to preserve units. (?????)
      LABEL(7)  = 'BWAH'; VALUE(7)  = BWAH / 10. 
      LABEL(8)  = 'HWUM'; VALUE(8)  = EYEWT       !unit eye weight g/unit
      LABEL(9)  = 'H#AM'; VALUE(9)  = GPSM        !# eyes/m2 at maturity
      LABEL(10) = 'H#UM'; VALUE(10) = GPP         !# eyes/fruit at maturity
      LABEL(11) = 'NFXM'; VALUE(11) = 0.0         !WTNFX*10.
      LABEL(12) = 'NUCM'; VALUE(12) = WTNUP*10.
      LABEL(13) = 'CNAM'; VALUE(13) = WTNCAN*10.
      LABEL(14) = 'GNAM'; VALUE(14) = WTNGRN*10.
      LABEL(15) = 'PWAM'; VALUE(15) = CRWNWT*10.
      LABEL(16) = 'LAIX'; VALUE(16) = MAXLAI
      LABEL(17) = 'HIAM'; VALUE(17) = HI
      LABEL(18) = 'EDAT'; VALUE(18) = FLOAT(YREMRG)

!     Send labels and values to OPSUM
      CALL SUMVALS (SUMNUM, LABEL, VALUE) 

      CALL OPVIEW(CONTROL, 
     &    Biomass_kg_ha, ACOUNT, DESCRIP, IDETO, LN, 
     &    Measured, PlantStres, Simulated, STGDOY, 
     &    STNAME, WTNCAN, LAI, NINT(YIELD), YRPLT, ISTAGE)

!     Send Measured and Simulated datat to OPSUM
      CALL EvaluateDat (ACOUNT, Measured, Simulated, DESCRIP, OLAP) 

!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!***********************************************************************
      RETURN
      END SUBROUTINE Aloha2_OPHARV
!=======================================================================



