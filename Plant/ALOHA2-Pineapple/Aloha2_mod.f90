!=======================================================================
!  MODULE Aloha2_mod
!  03/28/2017 CHP Written
!  09/05/2020 JVJ Stages inclusion for Overview.    
!=======================================================================

      MODULE Aloha2_mod
!     Contains data definitions for Aloha2 Pineaple model
      USE ModuleDefs

!     Data construct for control variables
      TYPE Aloha2Cul_type
        CHARACTER*6 VARNO,ECONO
        CHARACTER*12 CULfile
        CHARACTER*16 VRNAME
        CHARACTER*80 CULpath
        REAL TC, P1, P2, P3, P4, P5, P6, P7, P8
        REAL G1, G2, G3, PHINT
      END TYPE Aloha2Cul_type

      TYPE Aloha2Spe_type
        CHARACTER*12 SPEfile
        CHARACTER*80 SPEpath
        REAL CONV, FDMC
        REAL TBASV, TOPTV, TTOPV
        REAL TBASR, TOPTR, TTOPR
        REAL LIFAC
        REAL RWEP, PORM, RWMX, RLWR
        REAL CMFC
        REAL, DIMENSION(10) :: CO2X, CO2Y
      End Type Aloha2Spe_type

      TYPE Aloha2Sow_type
        REAL PLTPOP, SDEPTH, SDWTPL, PLANTSIZE
        INTEGER NFORCING, NDOF, PMTYPE, ForcingYRDOY
      End Type Aloha2Sow_type

      Type (Aloha2Cul_type) CULTIVAR
      Type (Aloha2Spe_type) SPECIES
      Type (Aloha2Sow_type) PLANTING

      CHARACTER*10 STNAME(20)
!     REAL TEMPM

      DATA STNAME/    &               
        'Foliar C1 ', &   ! 1   
        'Foliar C2 ', &   ! 2   
        'Foliar C3 ', &   ! 3    
        'Forcing   ', &   ! 4   
        'Open Heart', &   ! 5   
        'EarlyAnthe', &   ! 6   
        'LastAnthes', &   ! 7   
        'PhMaturity', &   ! 8   
        'Fruit Harv', &   ! 9   
        '          ', &   !10
        'Planting  ', &   !11         
        'Root Init.', &   !12
        'Leaf Emerg', &   !13
        'Start Sim ', &   !14
        'End Sim   ', &   !15
        '          ', &   !16
        '          ', &   !17
        '          ', &   !18
        '          ', &   !19
        'Harvest   '/     !20



      CONTAINS

!=======================================================================
!  Aloha2_IPCROP, Subroutine
!  Read crop parameters from PIALO980.SPE
!-----------------------------------------------------------------------
!  Revision history
!
!  06/15/1994 PWW Original written
!  03/29/2017 CGO Revised for v4.6
!  10/29/2021 FO  Fixed unit conflict on species file for gfortran
!=======================================================================
      SUBROUTINE Aloha2_IPCROP (CONTROL)

      IMPLICIT    NONE
      SAVE

      CHARACTER*6 SECTION
      CHARACTER*8, PARAMETER :: ERRKEY = 'IP_ALOHA' 
      CHARACTER*200 FILESPE, TEXTLINE
      INTEGER ERR, FOUND, I, ISECT, L, LNUM, LUNIO, LUNSPE

      TYPE (ControlType) CONTROL
      LUNIO  = CONTROL % LUNIO


!-----------------------------------------------------------------------
!   Open species file, parse headers, and read data.
    FILESPE = TRIM(SPECIES%SPEpath) // TRIM(SPECIES%SPEfile)
    CALL GETLUN('FILEC', LUNSPE)
    OPEN (LUNSPE, FILE = FILESPE, ACTION = 'READ', IOSTAT=ERR)
    IF (ERR .NE. 0) THEN
       CALL ERROR(ERRKEY,ERR,FILESPE,0)
    END IF

    SECTION = "*SPECI"
    CALL FIND(LUNSPE,SECTION,LNUM,FOUND)

    L=0
    DO WHILE (.TRUE.)  
      CALL IGNORE (LUNSPE, LNUM, ISECT, TEXTLINE)
      
      SELECT CASE(ISECT)
      CASE(0); EXIT   !End of file 

!     -----------------------------------------
!     Read species file parameters
      CASE(1);        !Data record 
        L = L + 1
        SELECT CASE(L)
        CASE (1); READ(TEXTLINE,'(2F7.0)',IOSTAT=ERR)       &
        SPECIES % CONV,  SPECIES % FDMC
        CASE (2); READ(TEXTLINE,'(3F7.0)',IOSTAT=ERR)       &
        SPECIES % TBASV, SPECIES % TOPTV, SPECIES % TTOPV
        CASE (3); READ(TEXTLINE,'(3F7.0)',IOSTAT=ERR)       & 
        SPECIES % TBASR, SPECIES % TOPTR, SPECIES % TTOPR
        CASE (4); READ(TEXTLINE,'(1F7.0)',IOSTAT=ERR)       &
        SPECIES % LIFAC
        CASE (5); READ(TEXTLINE,'(4F7.0)',IOSTAT=ERR)       &
        SPECIES % RWEP, SPECIES % PORM, SPECIES % RWMX,     & 
        SPECIES % RLWR
        CASE (6); READ(TEXTLINE,'(2F7.0)',IOSTAT=ERR)       &
        SPECIES % CMFC
        CASE (7:16); I = L - 6
                  READ(TEXTLINE,'(2F7.0)',IOSTAT=ERR)       &
                  SPECIES % CO2X(I), SPECIES % CO2Y(I)
        END SELECT !L = Data line selection
        IF (ERR .NE. 0) GOTO 100
!     -----------------------------------------

      CASE(2); EXIT   !End of section 
      END SELECT !ISECT
    ENDDO
    CLOSE(LUNSPE)

    CLOSE(LUNSPE)

    RETURN

100 CALL ERROR(ERRKEY,ERR,FILESPE,LNUM)   

!   -------------------------------------------
!   No ecotype file for now
!      READ(LUNIO,51,IOSTAT=ERR) FILEE, PATHEC; LNUM = LNUM + 1
!   51 FORMAT(15X,A12,1X,A80)
!      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEIO,LNUM)

      RETURN
      END SUBROUTINE Aloha2_IPCROP

!=======================================================================
!  Aloha2_IpPlant, Subroutine
!  Read crop parameters from PIALO980.SPE
!-----------------------------------------------------------------------
!  Revision history
!
!  06/15/1994 PWW Original written
!  03/29/2017 CHP Revised for v4.6
!=======================================================================
      SUBROUTINE Aloha2_IpPlant (CONTROL)

      IMPLICIT NONE
      SAVE

      CHARACTER*6, PARAMETER :: ERRKEY = 'IPPLNT'
      CHARACTER*5 ChemType
      CHARACTER*6 SECTION
      CHARACTER*12 FILEIO
      INTEGER LUNIO, ERR, LNUM, LINC, FOUND, ChemYRDOY
      TYPE (ControlType) CONTROL

      FILEIO = CONTROL % FILEIO
      
!     -----------------------------------------------------------------
!     Read input file name (ie. DSSAT45.INP) and path
      CALL GETLUN('FILEIO', LUNIO)
      OPEN (LUNIO, FILE = CONTROL % FILEIO, STATUS = 'OLD',       &
       IOSTAT=ERR)  
      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEIO,0)
      REWIND (LUNIO)
      READ(LUNIO,50,IOSTAT=ERR) SPECIES % SPEfile,                &
      SPECIES % SPEpath; LNUM = 7  !Species file, path
   50 FORMAT(//////,15X,A12,1X,A80)
      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEIO,LNUM)

!      READ(LUNIO,51,IOSTAT=ERR) FILEE, PATHER; LNUM = LNUM + 1  !Ecotype file name 
!      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEIO,LNUM)

      READ(LUNIO,51,IOSTAT=ERR) CULTIVAR % CULfile,               &
      CULTIVAR % CULpath; LNUM = LNUM + 2 
      IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEIO,LNUM)        !Cultivar file, path
   51 FORMAT(/,15X,A12,1X,A80)

!     -----------------------------------------------------------------
!     Read Planting Details Section
      SECTION = '*PLANT'
      CALL FIND(LUNIO, SECTION, LINC, FOUND) ; LNUM = LNUM + LINC
      IF (FOUND .EQ. 0) THEN
        CALL ERROR(SECTION, 42, FILEIO, LNUM)
      ELSE
        READ (LUNIO,70) PLANTING % PLTPOP, PLANTING % SDEPTH,     &
        PLANTING % SDWTPL !, PLANTING % NFORCING, & 
!                        PLANTING % PLANTSIZE, PLANTING % NDOF, PLANTING % PMTYPE
   70   FORMAT (24X,F6.0,24X,2F6.0,24X,I6,F6.0,2I6)
      ENDIF

!PLANTING DETAILS   
!        1         2         3         4         5         6         7         8         9        10        11        12
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
!  1989166     -99   2.6   2.6     S     R   53.    0.   5.0  595.  -99. -99.0 -99.0   0.0     2 551.2   459     0
!    yrplt   iemrg      pltpop  plme  plds       azir sdepth      sdage  atemp  plph    nforcing        ndof  
!                 plants                  rowspc           sdwtpl                   sprlap   plantsize       pmtype 

!     May want to remove PMTYPE and PLANTSIZE 
      PLANTING % PMTYPE = 0
      PLANTING % NFORCING = 0
      PLANTING % PLANTSIZE = 0.0
      PLANTING % PMTYPE = 0

!     -----------------------------------------------------------------
!     Read Chemical applications section to get forcing date
      SECTION = '*CHEMI'
      CALL FIND(LUNIO, SECTION, LINC, FOUND) ; LNUM = LNUM + LINC
      PLANTING % NFORCING = 0
      IF (FOUND .NE. 0) THEN
        DO WHILE (.TRUE.)
          READ (LUNIO,'(3X,I7,1X,A5)', IOSTAT=ERR) ChemYRDOY,ChemType
          IF (ERR .NE. 0) EXIT
          IF (ChemType == 'CH101' .OR. ChemTYPE == 'CH102') THEN
            PLANTING % NFORCING = 2 
            PLANTING % ForcingYRDOY = ChemYRDOY
            EXIT
          ENDIF
        ENDDO
      ENDIF

!*CHEMICALS          
!   1990350 CH102  1.00 AP006 -99.0   -99  -99     
                                
!     -----------------------------------------------------------------
 !     Read crop cultivar coefficients
       SECTION = '*CULTI'
       CALL FIND(LUNIO, SECTION, LINC, FOUND) ; LNUM = LNUM + LINC
       IF (FOUND .EQ. 0) THEN
           CALL ERROR(SECTION, 42, FILEIO, LNUM)
       ELSE
         READ (LUNIO,1800,IOSTAT=ERR)                              & 
            CULTIVAR % VARNO, CULTIVAR % VRNAME, CULTIVAR % ECONO, &
            CULTIVAR % TC, CULTIVAR % P1, CULTIVAR % P2,           &
            CULTIVAR % P3, CULTIVAR % P4, CULTIVAR % P5,           & 
            CULTIVAR % P6, CULTIVAR % P7, CULTIVAR % P8,           &  
            CULTIVAR % G1, CULTIVAR % G2, CULTIVAR % G3,           & 
            CULTIVAR % PHINT  
 1800    FORMAT (A6,1X,A16,1X,A6,1X,15F6.0)
!B0066 SC-F153          IB0001   60.0  629.  381. 2640.  400.  60.0  200.  14.0  40.0
         IF (ERR .NE. 0) CALL ERROR(ERRKEY,ERR,FILEIO,LNUM)
       ENDIF
      CLOSE(LUNIO)

      RETURN
      END SUBROUTINE Aloha2_IpPlant

!=======================================================================

!=======================================================================

      END MODULE Aloha2_mod
!======================================================================

