!> @brief USER INITIALISES THE WATER DEPTH H
      SUBROUTINE USER_CONDIN_H

      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D
      USE TPXO
      USE OKADA

      USE DECLARATIONS_SPECIAL

      IMPLICIT NONE


      INTEGER I

      DOUBLE PRECISION ETA0


!-----------------------------------------------------------------------
!
!     TUANZHOU L2-B1 STABILISED VERTICAL EROSION
!
!     DUAL WATER LEVEL INITIAL CONDITION
!
!     FIXED MOTHER MESH
!     INITIAL BREACH WIDTH = 20 m
!
!     Y < 0:
!          LAKE SIDE
!          FREE SURFACE = 31.5 m
!
!     Y >= 0:
!          POLDER SIDE
!          FREE SURFACE = 29.5 m
!
!     INITIAL WATER DEPTH:
!          H = MAX(ETA0 - ZF, 0)
!
!     NO BED CHANGE IS APPLIED DURING INITIALISATION.
!
!-----------------------------------------------------------------------


      DO I = 1,NPOIN


        IF(Y(I).LT.0.D0) THEN

          ETA0 = 31.5D0

        ELSE

          ETA0 = 29.5D0

        ENDIF


        H%R(I)=ETA0-ZF%R(I)


        IF(H%R(I).LT.0.D0) THEN

          H%R(I)=0.D0

        ENDIF


      ENDDO



!-----------------------------------------------------------------------
!
!     INITIAL VELOCITY
!
!     U = 0
!     V = 0
!
!-----------------------------------------------------------------------

      CALL OS('X=0     ',X=U)

      CALL OS('X=0     ',X=V)



      RETURN

      END


!=======================================================================
!=======================================================================


!> @brief USER TIME-STEP INTERFACE FOR DYNAMIC BOTTOM CONTROL
      SUBROUTINE USER_UTIMP_TELEMAC2D
     &(LTL,ATL,GRADEBL,GRAPRDL,LISDEBL,LISPRDL)

      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D

      USE DECLARATIONS_SPECIAL

      IMPLICIT NONE


      DOUBLE PRECISION, INTENT(IN) :: ATL

      INTEGER, INTENT(IN) :: LTL
      INTEGER, INTENT(IN) :: GRADEBL
      INTEGER, INTENT(IN) :: GRAPRDL
      INTEGER, INTENT(IN) :: LISDEBL
      INTEGER, INTENT(IN) :: LISPRDL


!-----------------------------------------------------------------------
!
!     TUANZHOU L2-B1 STABILISED VERTICAL EROSION
!
!     OFFICIAL TELEMAC USER TIME-STEP INTERFACE.
!
!     THIS ROUTINE DOES NOT CONTAIN THE EROSION LAW.
!     IT ONLY CALLS THE PROJECT-LEVEL DYNAMIC-BOTTOM ROUTINE.
!
!-----------------------------------------------------------------------


      CALL USER_DYNAMIC_BOTTOM(LTL,ATL)


      RETURN

      END


!=======================================================================
!=======================================================================


!> @brief PROJECT-LEVEL STABILISED VERTICAL EROSION MODEL
      SUBROUTINE USER_DYNAMIC_BOTTOM(LTL,ATL)

      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D

      USE DECLARATIONS_SPECIAL

      IMPLICIT NONE


      DOUBLE PRECISION, INTENT(IN) :: ATL

      INTEGER, INTENT(IN) :: LTL


      INTEGER I
      INTEGER NACTIVE
      INTEGER NERODED
      INTEGER NWET_ELIGIBLE
      INTEGER NLIMITED
      INTEGER N_WET_MIN
      INTEGER, ALLOCATABLE, SAVE :: WET_COUNT(:)
      LOGICAL, SAVE :: WET_COUNT_INITIALISED = .FALSE.

      DOUBLE PRECISION RHO_W
      DOUBLE PRECISION GRAV_B0
      DOUBLE PRECISION K_STRICKLER

      DOUBLE PRECISION HMIN_EROSION
      DOUBLE PRECISION TAU_C
      DOUBLE PRECISION K_D
      DOUBLE PRECISION POROSITY
      DOUBLE PRECISION Z_MIN
      DOUBLE PRECISION DZ_MAX_STEP

      DOUBLE PRECISION X_MIN_EROSION
      DOUBLE PRECISION X_MAX_EROSION
      DOUBLE PRECISION Y_MIN_EROSION
      DOUBLE PRECISION Y_MAX_EROSION

      DOUBLE PRECISION HI
      DOUBLE PRECISION VEL2
      DOUBLE PRECISION CF_LOCAL
      DOUBLE PRECISION TAUB
      DOUBLE PRECISION EROSION_RATE
      DOUBLE PRECISION DZ
      DOUBLE PRECISION DZ_RAW
      DOUBLE PRECISION Z_OLD
      DOUBLE PRECISION Z_NEW

      DOUBLE PRECISION MAX_TAUB
      DOUBLE PRECISION MAX_EROSION_RATE
      DOUBLE PRECISION MAX_DZ
      DOUBLE PRECISION MAX_DZ_RAW
      DOUBLE PRECISION TOTAL_DZ

      LOGICAL, SAVE :: MESSAGE_INITIAL = .FALSE.
      LOGICAL, SAVE :: MESSAGE_1       = .FALSE.
      LOGICAL, SAVE :: MESSAGE_10      = .FALSE.
      LOGICAL, SAVE :: MESSAGE_30      = .FALSE.
      LOGICAL, SAVE :: MESSAGE_60      = .FALSE.


!-----------------------------------------------------------------------
!
!     TUANZHOU L2-B1
!
!     PURPOSE:
!
!       DEVELOP AND VERIFY THE MINIMUM PHYSICAL VERTICAL-EROSION
!       LOOP:
!
!         H,U,V
!           ->
!         BED SHEAR STRESS
!           ->
!         EXCESS-SHEAR EROSION RATE
!           ->
!         BOTTOM LOWERING
!           ->
!         HYDRODYNAMIC FEEDBACK
!
!
!     IMPORTANT:
!
!       THIS IS A MODEL-DEVELOPMENT TEST.
!
!       K_D AND TAU_C BELOW ARE TEST PARAMETERS ONLY.
!       THEY ARE NOT CALIBRATED TUANZHOU SOIL PARAMETERS.
!
!
!     NO PRESCRIBED BREACH WIDENING IS USED.
!
!     NO SLOPE-STABILITY FAILURE IS USED.
!
!     NO AVALANCHE OR TERRAIN SMOOTHING IS USED.
!
!     B1 PER-TIME-STEP DZ LIMITER IS USED.
!
!     B1 WET/DRY TEMPORAL RELAXATION IS USED.
!
!     NO MASS-COMPENSATION MODEL IS USED.
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!     CONSTANTS, B0 TEST PARAMETERS, AND B1 CONTROLS
!
!-----------------------------------------------------------------------

      RHO_W       = 1000.D0

      GRAV_B0     = 9.81D0

      K_STRICKLER = 40.D0


!     MINIMUM WATER DEPTH REQUIRED FOR EROSION

      HMIN_EROSION = 0.01D0


!     TEST CRITICAL SHEAR STRESS [Pa]

      TAU_C       = 5.D0


!     TEST ERODIBILITY COEFFICIENT [m/(Pa s)]

      K_D         = 1.D-5

!     DAZZI (2019) BED POROSITY, INITIAL TEST VALUE
      POROSITY    = 0.4D0


!     MINIMUM ALLOWED BOTTOM ELEVATION [m]

      Z_MIN       = 27.D0

!     B1: REQUIRE 5 CONSECUTIVE WET STEPS BEFORE EROSION.
!     WITH DT=0.1 s, THIS CORRESPONDS TO 0.5 s.

      N_WET_MIN    = 5

!     B1: MAXIMUM BOTTOM LOWERING PER TIME STEP [m].
!     THIS IS ABOVE THE VALIDATED B0 BASELINE DZ AND SHOULD
!     REMAIN INACTIVE UNDER NORMAL B0-LIKE CONDITIONS.

      DZ_MAX_STEP  = 1.D-3


!-----------------------------------------------------------------------
!
!     FIXED EROSION-ACTIVATION ZONE
!
!     THIS ZONE ONLY DEFINES WHERE EROSION IS ALLOWED TO OPERATE.
!
!     IT DOES NOT PRESCRIBE THE BREACH WIDTH.
!
!     WHETHER EROSION ACTUALLY OCCURS IS DETERMINED BY:
!
!                 TAUB > TAU_C
!
!-----------------------------------------------------------------------

      X_MIN_EROSION = -113.D0

      X_MAX_EROSION =  113.D0

      Y_MIN_EROSION =  -25.D0

      Y_MAX_EROSION =   35.D0


!-----------------------------------------------------------------------
!
!     INITIALISE B1 CONSECUTIVE-WET COUNTER
!
!-----------------------------------------------------------------------

      IF(.NOT.WET_COUNT_INITIALISED) THEN
        ALLOCATE(WET_COUNT(NPOIN))
        WET_COUNT = 0
        WET_COUNT_INITIALISED = .TRUE.
      ENDIF


!-----------------------------------------------------------------------
!
!     INITIAL DIAGNOSTIC MESSAGE
!
!-----------------------------------------------------------------------

      IF(.NOT.MESSAGE_INITIAL) THEN

        WRITE(LU,*)
        WRITE(LU,*) '==============================================='
        WRITE(LU,*) 'L2-B1 STABILISED VERTICAL EROSION IS ACTIVE'
        WRITE(LU,*) 'USER_UTIMP_TELEMAC2D HOOK IS ACTIVE'
        WRITE(LU,*) 'NO PRESCRIBED BREACH WIDENING IS USED'
        WRITE(LU,*) 'NO BUILT-IN BREACH MODULE IS USED'
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*) 'RHO_W              = ',RHO_W
        WRITE(LU,*) 'GRAVITY             = ',GRAV_B0
        WRITE(LU,*) 'STRICKLER K         = ',K_STRICKLER
        WRITE(LU,*) 'HMIN_EROSION        = ',HMIN_EROSION
        WRITE(LU,*) 'TEST TAU_C          = ',TAU_C
        WRITE(LU,*) 'TEST K_D            = ',K_D
        WRITE(LU,*) 'POROSITY            = ',POROSITY
        WRITE(LU,*) 'Z_MIN               = ',Z_MIN
        WRITE(LU,*) 'N_WET_MIN           = ',N_WET_MIN
        WRITE(LU,*) 'DZ_MAX_STEP [m]      = ',DZ_MAX_STEP
        WRITE(LU,*) 'X EROSION RANGE     = ',
     &               X_MIN_EROSION,X_MAX_EROSION
        WRITE(LU,*) 'Y EROSION RANGE     = ',
     &               Y_MIN_EROSION,Y_MAX_EROSION
        WRITE(LU,*) '==============================================='
        WRITE(LU,*)

        MESSAGE_INITIAL = .TRUE.

      ENDIF


!-----------------------------------------------------------------------
!
!     INITIALISE STEP DIAGNOSTICS
!
!-----------------------------------------------------------------------

      NACTIVE          = 0

      NERODED          = 0
      NWET_ELIGIBLE     = 0
      NLIMITED          = 0

      MAX_TAUB         = 0.D0

      MAX_EROSION_RATE = 0.D0

      MAX_DZ           = 0.D0
      MAX_DZ_RAW       = 0.D0

      TOTAL_DZ         = 0.D0


!-----------------------------------------------------------------------
!
!     B1 STABILISED VERTICAL-EROSION LOOP
!
!     TELEMAC STRICKLER FRICTION LAW:
!
!       CF = 2*g / (K^2 * H^(1/3))
!
!     BED SHEAR STRESS:
!
!       TAUB = RHO_W * 0.5 * CF * (U^2 + V^2)
!
!     WHICH IS EQUIVALENT TO:
!
!       TAUB =
!         RHO_W * g * (U^2 + V^2)
!         / (K^2 * H^(1/3))
!
!
!     EXCESS-SHEAR EROSION LAW:
!
!       E = K_D * (TAUB - TAU_C),  IF TAUB > TAU_C
!
!       E = 0,                     OTHERWISE
!
!
!     BOTTOM UPDATE:
!
!       ZF_NEW = ZF_OLD - E * DT / (1 - POROSITY)
!       E = K_D * MAX(TAUB - TAU_C, 0)
!
!     WITH:
!
!       ZF_NEW >= Z_MIN
!
!-----------------------------------------------------------------------

      DO I = 1,NPOIN

!       CHECK FIXED EROSION-ACTIVATION ZONE

        IF((X(I).GE.X_MIN_EROSION).AND.
     &     (X(I).LE.X_MAX_EROSION).AND.
     &     (Y(I).GE.Y_MIN_EROSION).AND.
     &     (Y(I).LE.Y_MAX_EROSION)) THEN

          HI = H%R(I)

!         B1 WET/DRY TEMPORAL RELAXATION

          IF(HI.GT.HMIN_EROSION) THEN
            NACTIVE = NACTIVE + 1
            IF(WET_COUNT(I).LT.N_WET_MIN) THEN
              WET_COUNT(I) = WET_COUNT(I) + 1
            ENDIF
          ELSE
            WET_COUNT(I) = 0
          ENDIF

          IF(WET_COUNT(I).GE.N_WET_MIN) THEN

            NWET_ELIGIBLE = NWET_ELIGIBLE + 1

            VEL2 = U%R(I)*U%R(I) + V%R(I)*V%R(I)

!           LOCAL FRICTION COEFFICIENT CONSISTENT WITH
!           TELEMAC STRICKLER LAW

            CF_LOCAL =
     &        2.D0*GRAV_B0 /
     &        (K_STRICKLER*K_STRICKLER*HI**(1.D0/3.D0))

!           BED SHEAR STRESS [Pa]

            TAUB = RHO_W*0.5D0*CF_LOCAL*VEL2

            IF(TAUB.GT.MAX_TAUB) MAX_TAUB = TAUB

!           EXCESS-SHEAR EROSION

            IF((TAUB.GT.TAU_C).AND.
     &         (ZF%R(I).GT.Z_MIN)) THEN

              EROSION_RATE = K_D*(TAUB-TAU_C)/(1.D0-POROSITY)

!             RAW BOTTOM CHANGE BEFORE B1 LIMITER

              DZ_RAW = EROSION_RATE*DT

              IF(DZ_RAW.GT.MAX_DZ_RAW) MAX_DZ_RAW = DZ_RAW

!             B1 PER-TIME-STEP TERRAIN-CHANGE LIMITER

              DZ = MIN(DZ_RAW,DZ_MAX_STEP)

              IF(DZ_RAW.GT.DZ_MAX_STEP) NLIMITED = NLIMITED + 1

              Z_OLD = ZF%R(I)
              Z_NEW = Z_OLD-DZ

!             NEVER LOWER THE BOTTOM BELOW Z_MIN

              IF(Z_NEW.LT.Z_MIN) Z_NEW = Z_MIN

!             ACTUAL APPLIED BOTTOM CHANGE

              DZ = Z_OLD-Z_NEW

              IF(DZ.GT.0.D0) THEN
                ZF%R(I) = Z_NEW
                NERODED = NERODED + 1
                TOTAL_DZ = TOTAL_DZ + DZ

                IF(EROSION_RATE.GT.MAX_EROSION_RATE) THEN
                  MAX_EROSION_RATE = EROSION_RATE
                ENDIF

                IF(DZ.GT.MAX_DZ) MAX_DZ = DZ
              ENDIF

            ENDIF
          ENDIF

        ELSE

          WET_COUNT(I) = 0

        ENDIF

      ENDDO



!-----------------------------------------------------------------------
!
!     DIAGNOSTIC OUTPUT
!
!     THESE MESSAGES ARE FOR THE FIRST B1 DEVELOPMENT TEST.
!
!-----------------------------------------------------------------------


      IF((ATL.GE.1.D0).AND.(.NOT.MESSAGE_1)) THEN

        WRITE(LU,*)
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*) 'L2-B1 CHECK AT ABOUT 1 s'
        WRITE(LU,*) 'ATL                = ',ATL
        WRITE(LU,*) 'ACTIVE WET NODES   = ',NACTIVE
        WRITE(LU,*) 'WET-ELIGIBLE NODES = ',NWET_ELIGIBLE
        WRITE(LU,*) 'ERODED NODES       = ',NERODED
        WRITE(LU,*) 'DZ-LIMITED NODES   = ',NLIMITED
        WRITE(LU,*) 'MAX TAUB [Pa]      = ',MAX_TAUB
        WRITE(LU,*) 'MAX E [m/s]        = ',MAX_EROSION_RATE
        WRITE(LU,*) 'MAX RAW DZ [m]     = ',MAX_DZ_RAW
        WRITE(LU,*) 'MAX DZ/STEP [m]    = ',MAX_DZ
        WRITE(LU,*) 'SUM DZ/STEP [m]    = ',TOTAL_DZ
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*)

        MESSAGE_1 = .TRUE.

      ENDIF


      IF((ATL.GE.10.D0).AND.(.NOT.MESSAGE_10)) THEN

        WRITE(LU,*)
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*) 'L2-B1 CHECK AT ABOUT 10 s'
        WRITE(LU,*) 'ATL                = ',ATL
        WRITE(LU,*) 'ACTIVE WET NODES   = ',NACTIVE
        WRITE(LU,*) 'WET-ELIGIBLE NODES = ',NWET_ELIGIBLE
        WRITE(LU,*) 'ERODED NODES       = ',NERODED
        WRITE(LU,*) 'DZ-LIMITED NODES   = ',NLIMITED
        WRITE(LU,*) 'MAX TAUB [Pa]      = ',MAX_TAUB
        WRITE(LU,*) 'MAX E [m/s]        = ',MAX_EROSION_RATE
        WRITE(LU,*) 'MAX RAW DZ [m]     = ',MAX_DZ_RAW
        WRITE(LU,*) 'MAX DZ/STEP [m]    = ',MAX_DZ
        WRITE(LU,*) 'SUM DZ/STEP [m]    = ',TOTAL_DZ
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*)

        MESSAGE_10 = .TRUE.

      ENDIF


      IF((ATL.GE.30.D0).AND.(.NOT.MESSAGE_30)) THEN

        WRITE(LU,*)
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*) 'L2-B1 CHECK AT ABOUT 30 s'
        WRITE(LU,*) 'ATL                = ',ATL
        WRITE(LU,*) 'ACTIVE WET NODES   = ',NACTIVE
        WRITE(LU,*) 'WET-ELIGIBLE NODES = ',NWET_ELIGIBLE
        WRITE(LU,*) 'ERODED NODES       = ',NERODED
        WRITE(LU,*) 'DZ-LIMITED NODES   = ',NLIMITED
        WRITE(LU,*) 'MAX TAUB [Pa]      = ',MAX_TAUB
        WRITE(LU,*) 'MAX E [m/s]        = ',MAX_EROSION_RATE
        WRITE(LU,*) 'MAX RAW DZ [m]     = ',MAX_DZ_RAW
        WRITE(LU,*) 'MAX DZ/STEP [m]    = ',MAX_DZ
        WRITE(LU,*) 'SUM DZ/STEP [m]    = ',TOTAL_DZ
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*)

        MESSAGE_30 = .TRUE.

      ENDIF


      IF((ATL.GE.60.D0).AND.(.NOT.MESSAGE_60)) THEN

        WRITE(LU,*)
        WRITE(LU,*) '==============================================='
        WRITE(LU,*) 'L2-B1 CHECK AT ABOUT 60 s'
        WRITE(LU,*) 'ATL                = ',ATL
        WRITE(LU,*) 'ACTIVE WET NODES   = ',NACTIVE
        WRITE(LU,*) 'WET-ELIGIBLE NODES = ',NWET_ELIGIBLE
        WRITE(LU,*) 'ERODED NODES       = ',NERODED
        WRITE(LU,*) 'DZ-LIMITED NODES   = ',NLIMITED
        WRITE(LU,*) 'MAX TAUB [Pa]      = ',MAX_TAUB
        WRITE(LU,*) 'MAX E [m/s]        = ',MAX_EROSION_RATE
        WRITE(LU,*) 'MAX RAW DZ [m]     = ',MAX_DZ_RAW
        WRITE(LU,*) 'MAX DZ/STEP [m]    = ',MAX_DZ
        WRITE(LU,*) 'SUM DZ/STEP [m]    = ',TOTAL_DZ
        WRITE(LU,*) '==============================================='
        WRITE(LU,*)

        MESSAGE_60 = .TRUE.

      ENDIF


      RETURN

      END