!=======================================================================
!
!     COPYRIGHT(C) 2017 ALL RIGHTS RESERVED BY
!     HAECHEON CHOI
!     DEPARTMENT OF MECHANICAL & AEROSPACE ENGINEERING
!     SEOUL NATIONAL UNIVERSITY
!     E-MAIL: CHOI@SNU.AC.KR
!     URL: HTTP://TFC.SNU.AC.KR
!
!     NONCOMMERCIAL USE WITH COPYRIGHTED MARK
!
!     Q&A: LICACODEFORUM@GMAIL.COM
!
!=======================================================================
!
!     2ND COPYRIGHT(C) 2018 ALL RIGHT RESERVED BY SANGJOON LEE
!     2018.02.28. MODIFIED FOR IN-LAB USAGE
!                 SANGJOON LEE FROM ENERGY AND ENVIRONMENTL FLOW LAB.
!                 (HTTP://EEFLOW.SNU.AC.KR)
!                 DEPARTMENT OF MECHANICAL & AEROSPACE ENGINEERING
!                 SEOUL NATIONAL UNIVERSITY
!     ONLY FOR IN-LAB USE. DO NOT DISTRIBUTE FOR COMMERCIAL USE.
!
!=======================================================================
!
!     TAILORED FOR RIBLET CHANNEL'S DRAG REDUCTION EVALUATION
!     2022.09.08. SANGJOON LEE FROM CFD LAB.
!     DEPARTMENT OF MECHANICAL ENGINEERING
!     UNIVERSITY OF CALIFORNIA, BERKELEY
!
!=======================================================================
      PROGRAM SOLVER
!=======================================================================
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE
      INTEGER*8   ::  I,J,K,N
      REAL*8      ::  OMP_GET_WTIME,TMP
      CHARACTER(LEN=200),DIMENSION(2) :: ARGS ! 1ST: GRID_DIR, 2ND: IBMPRE_DIR

      ! CALL PRINT_REAL_TIME()                               ! AT MISC_INIT LIBRARY

      N = COMMAND_ARGUMENT_COUNT()
      IF (N.NE.2) THEN
        WRITE (*,*) 'SOLVER: PROVIDE TWO ARGS: GRID_DIR & IBMPRE_DIR'
        STOP
      ENDIF

      CALL GET_COMMAND_ARGUMENT(1,ARGS(1))
      CALL GET_COMMAND_ARGUMENT(2,ARGS(2))

      TOTAL_TIME_B=OMP_GET_WTIME()                         ! INTRINSIC SUBROUTINE
      CALL READSETTINGS()                                  ! AT MISC_INIT LIBRARY
      CALL READGEOM(TRIM(ADJUSTL(ARGS(1))))                ! AT MISC_INIT LIBRARY
      CALL READBCS(TRIM(ADJUSTL(ARGS(2))))                 ! AT MISC_INIT LIBRARY
      CALL ALLOINIT(TRIM(ADJUSTL(ARGS(2))))                ! AT MISC_INIT LIBRARY
      CALL ALLO_ARRAY()                                    ! AT MISC_INIT LIBRARY

!=============== INITIALIZATION FROM PRE-PROCESSED DATA
      IF (IBMON .EQ. 1) THEN
        CALL IBMPREREAD(TRIM(ADJUSTL(ARGS(2))))            ! AT MISC_INIT LIBRARY
        IF (ILES .EQ. 1) THEN
            CALL NUTZEROREAD(TRIM(ADJUSTL(ARGS(2))))       ! AT MISC_INIT LIBRARY
        ENDIF
        IF (ICONJG .EQ. 1) THEN
          CALL CONJGREAD(TRIM(ADJUSTL(ARGS(2))))           ! AT MISC_INIT LIBRARY
        ENDIF
      ENDIF
      IF (ILES .EQ. 1) CALL SGSFILTERINIT()                ! AT SGS       LIBRARY

!=============== FIELD LOADING
      IF (IREAD .EQ. 1) THEN
        CALL PREFLD()                                      ! AT MISC_INIT LIBRARY
      ELSE
        CALL MAKEFLD()                                     ! AT MISC_INIT LIBRARY
      ENDIF

      TMP = EPS_PTR
      EPS_PTR = 1.D-5
      CALL ADDPERTURB()               ! AT MISC_INIT LIBRARY
      EPS_PTR = TMP

      CALL WRITEFIELD()
      IF (ICH .EQ. 1) CALL MEANPG()                        ! AT SLV_MMTM  LIBRARY

!=============== MISCELLANEOUS SETTINGS
      CALL DTTIMEINIT()                                    ! AT MISC_INIT LIBRARY
      CALL RK3COEFINIT()                                   ! AT MISC_INIT LIBRARY
      CALL LHSINIT()                                       ! AT SLV_MMTM  LIBRARY
      CALL POISINIT()                                      ! AT POISS     LIBRARY
      CALL FTRFILES(1)
      CFLMAX = CFLFAC
      NV = 101
      NAV = 3001
      QVOL = 0.

!=============== TIME-DEPENDENT CALCULATION (MAIN SOLVER)
      DO M = 1, NTST
        TIME_BEGIN=OMP_GET_WTIME()                         ! INTRINSIC SUBROUTINE
        CALL STEPINIT()                                    ! AT SLV_AUXI  LIBRARY
        SUBDT = 0.D0
!%%%%%%%%%%%%%%%%      RK3 MAIN SOLVER      %%%%%%%%%%%%%%%%
        DO MSUB = 1, 3
          CALL SUBSTEPINIT(MSUB)                           ! AT SLV_AUXI  LIBRARY
          IF (ICH .EQ. 1) THEN
            CALL QVOLCALC(QVOL(1))                         ! AT SLV_AUXI  LIBRARY
          ENDIF

          SGSTIME_B(MSUB)=OMP_GET_WTIME()                  ! INTRINSIC SUBROUTINE
          IF (ILES .EQ. 1) THEN
            CALL SGSCALC()                                 ! AT SGS       LIBRARY
            ! IF (IHTRANS .EQ. 1) CALL SGSCALC_T()           ! AT SGS       LIBRARY
          ENDIF
          SGSTIME_E(MSUB)=OMP_GET_WTIME()                  ! INTRINSIC SUBROUTINE

          ! ! IF (IMOVINGON .EQ. 1) CALL FINDFORCING()         ! AT IBM_BODY  LIBRARY
          ! !                                                  ! (GLOBAL LIBRARY)
          RHSNLHSTIME_B(MSUB)=OMP_GET_WTIME()              ! INTRINSIC SUBROUTINE
          CALL RHSNLHS()                                   ! AT SLV_MMTM  LIBRARY
          RHSNLHSTIME_E(MSUB)=OMP_GET_WTIME()              ! INTRINSIC SUBROUTINE

          ALLOCATE(DIVGSUM(0:N1,0:N2,0:N3))
          ALLOCATE(PHI(0:N1,0:N2,0:N3))
          DIVGSUM = 0.
          PHI = 0.

          POISSTIME_B(MSUB)=OMP_GET_WTIME()                ! INTRINSIC SUBROUTINE
          CALL DIVGS(DIVGSUM)                              ! AT SLV_CONT  LIBRARY
          CALL POISSON(PHI,DIVGSUM)                        ! AT POISS     LIBRARY
          POISSTIME_E(MSUB)=OMP_GET_WTIME()                ! INTRINSIC SUBROUTINE

          IF (ICH .EQ. 1) THEN
            CALL QVOLCALC(QVOL(2))                         ! AT SLV_AUXI  LIBRARY
          ENDIF

          ! IF (ICH .EQ. 1) THEN
          !   CALL QVOLCORR()                                ! AT SLV_AUXI  LIBRARY
          ! ENDIF
          
          CALL UCALC(PHI)

          IF (ICH .EQ. 1) THEN
            CALL QVOLCALC(QVOL(0))                         ! AT SLV_AUXI  LIBRARY
          ENDIF

          CALL PCALC(PHI,DIVGSUM)

          IF (ICH .EQ. 1) CALL MEANPG()                    ! AT SLV_AUXI  LIBRARY
          IF (IBMON .EQ. 1) CALL LAGFORCE

          ! IF (IHTRANS .EQ. 1) THEN
          !   IF (ICH .EQ. 1) CALL TVOLCALC(TVOL(1))         ! AT SLV_AUXI LIBRARY
          !   CALL RHSNLHS_T()                               ! AT SLV_ENGY  LIBRARY
          !   IF (ICH .EQ. 1) CALL TVOLCALC(TVOL(2))         ! AT SLV_AUXI LIBRARY
          !   IF (ICH .EQ. 1) CALL TVOLCORR()                ! AT SLV_AUXI LIBRARY
          !   CALL TCALC()                                   ! AT SLV_AUXI LIBRARY
          !   IF (ICH .EQ. 1) CALL TVOLCALC(TVOL(0))         ! AT SLV_AUXI LIBRARY
          ! ENDIF

          DEALLOCATE(DIVGSUM,PHI)

        ENDDO
!%%%%%%%%%%%%%%%%  END OF RK3 MAIN SOLVER  %%%%%%%%%%%%%%%%

        TIME = TIME + DT

        CALL CONVERGENCE_CHECK()
        IF ((IBMON .EQ. 1) .AND. (MASSON .EQ. 1)) CALL MASSCHECK()

        IF (MOD(NTIME,NPIN) .EQ. 0) THEN
          IHIST = IHIST + 1
          CALL WRITEHISTORY()
        ENDIF

        IF ((NTRACE .GT. 0) .AND. (MOD(M,NTR) .EQ. 0)) CALL TRACER()
        ! IF (IBMON .EQ. 1) 
        CALL DRAGLIFT()
        ! IF (ICH .EQ. 1) 
        CALL RIBFORCE()
        ! IF (IHTRANS .EQ. 1) CALL NUSSELT()

        IF (MOD(NTIME,NPRINT) .EQ. 0) CALL WRITEFIELD()
        IF (IAVG .EQ. 1) THEN
          IF ((TIME .GE. 500.) .AND. (TIME-DT .LT. 500.)) NPRIAVG = 1
          CALL FIELD_AVG()
          NPRIAVG = 0
        ENDIF

        TIME_END=OMP_GET_WTIME()
        CALL WRITEFTRTIME()

        IF (NTIME .EQ. 1) I = -1

        IF (QVOL(2).LT..9*XL*YL*ZL) THEN
          TMP = TIME
        ELSE
          IF ((INT((TIME-TMP)/5.) .EQ. I+1).AND.(INT((TIME-TMP)/5.) .LT. 1)) THEN
            I = I+1
            IF (EPS_PTR .NE. 0.) CALL ADDPERTURBSINE()
          ENDIF
        ENDIF
        IF (I.EQ.0) CFLFAC = 3.

        IF (TIME.GE.3000.) EXIT

      ENDDO

! !=============== FINISH AND CLEAR-OFF
      CALL FTRFILES(0)

      IF (MOD(NTIME,NPRINT) .NE. 0) CALL WRITEFIELD()
      IF (IAVG .EQ. 1) THEN
        NPRIAVG = 1
        CALL FIELD_AVG()
        NPRIAVG = 0
      ENDIF

      ! CALL CPU_TIME(TOTAL_TIME_E)
      ! CALL TOTAL_TIME()

      ! CALL PRINT_REAL_TIME()
      CALL DEALLO_ARRAY()

      STOP
      END PROGRAM SOLVER
!=======================================================================
!=======================================================================
