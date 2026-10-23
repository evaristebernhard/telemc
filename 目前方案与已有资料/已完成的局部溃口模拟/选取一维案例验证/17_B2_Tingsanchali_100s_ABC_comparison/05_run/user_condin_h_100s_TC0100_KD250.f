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
!     TINGSANCHALI & CHINNARASRI (2001) RUN B-2
!
!     EROSION REFERENCE RUN
!
!     FULL EXPERIMENTAL FLUME:
!
!       -17.50 <= X <= 17.50 m
!         0.00 <= Y <=  1.00 m
!
!     DAM GEOMETRY:
!
!       UPSTREAM TOE       X = -2.40 m
!       CREST FRONT        X =  0.00 m
!       CREST BACK         X =  0.30 m
!       DOWNSTREAM TOE     X =  4.30 m
!       DAM HEIGHT             0.80 m
!
!     RUN B-2 INITIAL CONDITION:
!
!       X <= 0:
!         INITIAL FREE SURFACE = 0.85 m
!
!       X > 0:
!         INITIALLY DRY
!
!     INITIAL VELOCITY:
!
!       U = 0
!       V = 0
!
!     NO BED CHANGE IS APPLIED DURING INITIALISATION.
!
!-----------------------------------------------------------------------


      DO I = 1,NPOIN

        IF(X(I).LE.0.D0) THEN

          ETA0 = 0.85D0

          H%R(I) = ETA0-ZF%R(I)

          IF(H%R(I).LT.0.D0) THEN
            H%R(I) = 0.D0
          ENDIF

        ELSE

          H%R(I) = 0.D0

        ENDIF

      ENDDO


!-----------------------------------------------------------------------
!
!     INITIAL VELOCITY
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
!     TINGSANCHALI B-2 EROSION REFERENCE RUN
!
!     OFFICIAL TELEMAC USER TIME-STEP INTERFACE.
!
!     THIS ROUTINE ONLY CALLS THE PROJECT-LEVEL
!     DYNAMIC-BOTTOM ROUTINE.
!
!-----------------------------------------------------------------------


      CALL USER_DYNAMIC_BOTTOM(LTL,ATL)


      RETURN

      END


!=======================================================================
!=======================================================================


!> @brief B-2 EXCESS-SHEAR EROSION MODEL
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
      DOUBLE PRECISION GRAV_B2
      DOUBLE PRECISION K_STRICKLER

      DOUBLE PRECISION HMIN_EROSION

      DOUBLE PRECISION POROSITY
      DOUBLE PRECISION SOLID_FRACTION

      DOUBLE PRECISION TAU_C
      DOUBLE PRECISION K_D

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

      DOUBLE PRECISION SOLID_EROSION_RATE
      DOUBLE PRECISION BED_LOWERING_RATE

      DOUBLE PRECISION DZ
      DOUBLE PRECISION DZ_RAW

      DOUBLE PRECISION Z_OLD
      DOUBLE PRECISION Z_NEW

      DOUBLE PRECISION MAX_TAUB
      DOUBLE PRECISION MAX_SOLID_EROSION_RATE
      DOUBLE PRECISION MAX_BED_LOWERING_RATE
      DOUBLE PRECISION MAX_DZ
      DOUBLE PRECISION MAX_DZ_RAW
      DOUBLE PRECISION TOTAL_DZ

      LOGICAL, SAVE :: MESSAGE_INITIAL = .FALSE.
      LOGICAL, SAVE :: MESSAGE_1       = .FALSE.
      LOGICAL, SAVE :: MESSAGE_5       = .FALSE.
      LOGICAL, SAVE :: MESSAGE_10      = .FALSE.
      LOGICAL, SAVE :: MESSAGE_20      = .FALSE.
      LOGICAL, SAVE :: MESSAGE_30      = .FALSE.


!-----------------------------------------------------------------------
!
!     B-2 EROSION REFERENCE MODEL
!
!     PURPOSE:
!
!       ESTABLISH THE FORMAL EXCESS-SHEAR EROSION MODEL USED
!       FROM THE B-2 EXPERIMENTAL BENCHMARK ONWARD.
!
!
!     PHYSICAL CHAIN:
!
!       H,U,V
!         ->
!       BED SHEAR STRESS
!         ->
!       SOLID-MATERIAL EROSION RATE
!         ->
!       POROSITY-CORRECTED BED LOWERING
!         ->
!       ZF UPDATE
!         ->
!       HYDRODYNAMIC FEEDBACK
!
!
!     EXCESS-SHEAR LAW:
!
!       E_S = K_D * (TAUB - TAU_C)
!
!       FOR TAUB > TAU_C.
!
!       E_S = 0 OTHERWISE.
!
!
!     E_S IS DEFINED AS THE VOLUMETRIC REMOVAL RATE
!     OF SOLID MATERIAL PER UNIT BED AREA.
!
!
!     SOLID-VOLUME CONSERVATION:
!
!       (1-P) * (-DZB/DT) = E_S
!
!     THEREFORE:
!
!       -DZB/DT = E_S / (1-P)
!
!
!     IMPORTANT:
!
!       TAU_C AND K_D BELOW DEFINE A REFERENCE RUN ONLY.
!
!       THEY ARE NOT YET THE FINAL CALIBRATED B-2 PARAMETERS.
!
!       THE REFERENCE RUN IS USED TO VERIFY THE IMPLEMENTATION
!       BEFORE K_D - TAU_C SENSITIVITY / CALIBRATION.
!
!
!     NO SLOPE-STABILITY FAILURE IS USED.
!
!     NO AVALANCHE OR TERRAIN SMOOTHING IS USED.
!
!     NO PRESCRIBED BREACH WIDENING IS USED.
!
!     NO BUILT-IN TELEMAC BREACH MODULE IS USED.
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!     PHYSICAL CONSTANTS AND HYDRAULIC PARAMETERS
!
!-----------------------------------------------------------------------

      RHO_W       = 1000.D0

      GRAV_B2     = 9.81D0


!     RUN B-2:
!     MANNING n = 0.018
!
!     STRICKLER K = 1/n = 55.55555556

      K_STRICKLER = 55.55555556D0


!-----------------------------------------------------------------------
!
!     EROSION PARAMETERS
!
!-----------------------------------------------------------------------


!     MINIMUM WATER DEPTH REQUIRED FOR EROSION [m]

      HMIN_EROSION = 0.01D0


!     MATERIAL POROSITY [-]
!
!     REFERENCE VALUE FOR THE SAND-DAM MATERIAL.

      POROSITY = 0.40D0

      SOLID_FRACTION = 1.D0-POROSITY


!     REFERENCE CRITICAL SHEAR STRESS [Pa]
!
!     THIS VALUE IS USED FOR THE REFERENCE RUN ONLY.
!     TAU_C WILL BE EXAMINED IN THE FOLLOWING
!     SENSITIVITY / CALIBRATION STAGE.

      TAU_C = 0.1D0


!     REFERENCE ERODIBILITY COEFFICIENT [m/(Pa s)]
!
!     5.0E-4 m/(Pa s)
!       =
!     500 cm3/(N s)
!
!     THIS VALUE IS USED FOR THE REFERENCE RUN ONLY.
!     IT IS NOT THE FINAL CALIBRATED VALUE.
      K_D = 2.5D-4
!-----------------------------------------------------------------------
!
!     BED-ELEVATION CONTROL
!
!-----------------------------------------------------------------------


!     FIXED FLUME-BED ELEVATION [m]
!
!     THE ERODIBLE DAM CANNOT BE LOWERED BELOW
!     THE ORIGINAL FLUME BED.

      Z_MIN = 0.D0


!     REQUIRE 0.5 s OF CONSECUTIVE WETTING BEFORE EROSION.
!
!     DT = 0.01 s
!
!     50 * 0.01 = 0.50 s

      N_WET_MIN = 50


!     MAXIMUM ALLOWED BOTTOM LOWERING PER TIME STEP [m]
!
!     THIS IS A NUMERICAL SAFETY LIMITER.
!
!     IT MUST NOT CONTROL THE CALIBRATED EROSION RESPONSE.
!     THE NUMBER OF LIMITER ACTIVATIONS IS THEREFORE
!     EXPLICITLY REPORTED.

      DZ_MAX_STEP = 1.D-3


!-----------------------------------------------------------------------
!
!     EROSION-ACTIVATION REGION
!
!     ONLY THE INITIAL EXPERIMENTAL DAM BODY IS ERODIBLE.
!
!     UPSTREAM TOE   : X = -2.40 m
!     DOWNSTREAM TOE : X =  4.30 m
!
!     FULL FLUME WIDTH:
!
!       0 <= Y <= 1 m
!
!     THE FIXED FLUME BED OUTSIDE THE DAM IS NOT ERODIBLE.
!
!-----------------------------------------------------------------------

      X_MIN_EROSION = -2.40D0

      X_MAX_EROSION =  4.30D0

      Y_MIN_EROSION =  0.D0

      Y_MAX_EROSION =  1.D0


!-----------------------------------------------------------------------
!
!     INITIALISE CONSECUTIVE-WET COUNTER
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
        WRITE(LU,*) 'B-2 EROSION REFERENCE MODEL IS ACTIVE'
        WRITE(LU,*) 'USER_UTIMP_TELEMAC2D HOOK IS ACTIVE'
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*) 'K_D CALIBRATION RUN - NOT FINAL CALIBRATION'
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*) 'RHO_W               = ',RHO_W
        WRITE(LU,*) 'GRAVITY              = ',GRAV_B2
        WRITE(LU,*) 'STRICKLER K          = ',K_STRICKLER
        WRITE(LU,*) 'HMIN_EROSION         = ',HMIN_EROSION
        WRITE(LU,*) 'POROSITY             = ',POROSITY
        WRITE(LU,*) 'SOLID FRACTION       = ',SOLID_FRACTION
        WRITE(LU,*) 'REFERENCE TAU_C [Pa] = ',TAU_C
        WRITE(LU,*) 'REFERENCE K_D         = ',K_D
        WRITE(LU,*) 'Z_MIN                = ',Z_MIN
        WRITE(LU,*) 'N_WET_MIN            = ',N_WET_MIN
        WRITE(LU,*) 'DZ_MAX_STEP [m]      = ',DZ_MAX_STEP
        WRITE(LU,*) 'X EROSION RANGE      = ',
     &               X_MIN_EROSION,X_MAX_EROSION
        WRITE(LU,*) 'Y EROSION RANGE      = ',
     &               Y_MIN_EROSION,Y_MAX_EROSION
        WRITE(LU,*) 'NO SLOPE FAILURE'
        WRITE(LU,*) 'NO AVALANCHE'
        WRITE(LU,*) 'NO PRESCRIBED BREACH WIDENING'
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
      NWET_ELIGIBLE    = 0
      NERODED          = 0
      NLIMITED         = 0

      MAX_TAUB                 = 0.D0
      MAX_SOLID_EROSION_RATE   = 0.D0
      MAX_BED_LOWERING_RATE    = 0.D0
      MAX_DZ                   = 0.D0
      MAX_DZ_RAW               = 0.D0
      TOTAL_DZ                 = 0.D0


!-----------------------------------------------------------------------
!
!     B-2 EROSION LOOP
!
!
!     TELEMAC STRICKLER FRICTION LAW:
!
!       CF =
!         2*g /
!         (K^2 * H^(1/3))
!
!
!     BED SHEAR STRESS:
!
!       TAUB =
!         RHO_W * 0.5 * CF * (U^2 + V^2)
!
!
!     SOLID-MATERIAL EROSION RATE:
!
!       E_S =
!         K_D * (TAUB - TAU_C)
!
!
!     POROSITY-CORRECTED GEOMETRIC BED-LOWERING RATE:
!
!       E_Z =
!         E_S / (1-P)
!
!
!     RAW BED CHANGE:
!
!       DZ_RAW =
!         E_Z * DT
!
!
!     NUMERICAL SAFETY LIMIT:
!
!       DZ =
!         MIN(DZ_RAW,DZ_MAX_STEP)
!
!
!     BOTTOM UPDATE:
!
!       ZF_NEW =
!         MAX(Z_MIN,ZF_OLD-DZ)
!
!-----------------------------------------------------------------------

      DO I = 1,NPOIN


!       CHECK WHETHER NODE BELONGS TO INITIAL DAM REGION

        IF((X(I).GE.X_MIN_EROSION).AND.
     &     (X(I).LE.X_MAX_EROSION).AND.
     &     (Y(I).GE.Y_MIN_EROSION).AND.
     &     (Y(I).LE.Y_MAX_EROSION)) THEN


          HI = H%R(I)


!-----------------------------------------------------------------------
!
!         WET/DRY TEMPORAL RELAXATION
!
!-----------------------------------------------------------------------

          IF(HI.GT.HMIN_EROSION) THEN

            NACTIVE = NACTIVE + 1

            IF(WET_COUNT(I).LT.N_WET_MIN) THEN
              WET_COUNT(I) = WET_COUNT(I) + 1
            ENDIF

          ELSE

            WET_COUNT(I) = 0

          ENDIF


!-----------------------------------------------------------------------
!
!         EROSION IS ALLOWED ONLY AFTER 0.5 s
!         OF CONSECUTIVE WETTING.
!
!-----------------------------------------------------------------------

          IF(WET_COUNT(I).GE.N_WET_MIN) THEN

            NWET_ELIGIBLE = NWET_ELIGIBLE + 1


            VEL2 =
     &        U%R(I)*U%R(I) +
     &        V%R(I)*V%R(I)


!           LOCAL FRICTION COEFFICIENT CONSISTENT WITH
!           TELEMAC STRICKLER LAW

            CF_LOCAL =
     &        2.D0*GRAV_B2 /
     &        (K_STRICKLER*K_STRICKLER*
     &         HI**(1.D0/3.D0))


!           BED SHEAR STRESS [Pa]

            TAUB =
     &        RHO_W*0.5D0*CF_LOCAL*VEL2


            IF(TAUB.GT.MAX_TAUB) THEN
              MAX_TAUB = TAUB
            ENDIF


!-----------------------------------------------------------------------
!
!           EXCESS-SHEAR EROSION
!
!           ZF > Z_MIN PREVENTS EROSION OF THE
!           FIXED FLUME FLOOR.
!
!-----------------------------------------------------------------------

            IF((TAUB.GT.TAU_C).AND.
     &         (ZF%R(I).GT.Z_MIN)) THEN


!             SOLID-MATERIAL VOLUMETRIC EROSION RATE
!
!             UNIT:
!
!               m3/(m2 s) = m/s

              SOLID_EROSION_RATE =
     &          K_D*(TAUB-TAU_C)


!             CONVERT SOLID-MATERIAL REMOVAL TO
!             GEOMETRIC BED-LOWERING RATE

              BED_LOWERING_RATE =
     &          SOLID_EROSION_RATE /
     &          SOLID_FRACTION


              IF(SOLID_EROSION_RATE.GT.
     &           MAX_SOLID_EROSION_RATE) THEN

                MAX_SOLID_EROSION_RATE =
     &            SOLID_EROSION_RATE

              ENDIF


              IF(BED_LOWERING_RATE.GT.
     &           MAX_BED_LOWERING_RATE) THEN

                MAX_BED_LOWERING_RATE =
     &            BED_LOWERING_RATE

              ENDIF


!             RAW BOTTOM CHANGE BEFORE LIMITER

              DZ_RAW =
     &          BED_LOWERING_RATE*DT


              IF(DZ_RAW.GT.MAX_DZ_RAW) THEN
                MAX_DZ_RAW = DZ_RAW
              ENDIF


!             PER-TIME-STEP TERRAIN-CHANGE LIMITER

              DZ =
     &          MIN(DZ_RAW,DZ_MAX_STEP)


              IF(DZ_RAW.GT.DZ_MAX_STEP) THEN
                NLIMITED = NLIMITED + 1
              ENDIF


              Z_OLD = ZF%R(I)

              Z_NEW = Z_OLD-DZ


!             NEVER LOWER THE DAM BELOW
!             THE ORIGINAL FLUME FLOOR

              IF(Z_NEW.LT.Z_MIN) THEN
                Z_NEW = Z_MIN
              ENDIF


!             ACTUAL APPLIED BED CHANGE

              DZ = Z_OLD-Z_NEW


              IF(DZ.GT.0.D0) THEN

                ZF%R(I) = Z_NEW

                NERODED = NERODED + 1

                TOTAL_DZ =
     &            TOTAL_DZ+DZ


                IF(DZ.GT.MAX_DZ) THEN
                  MAX_DZ = DZ
                ENDIF

              ENDIF


            ENDIF

          ENDIF


        ELSE


!         NODE OUTSIDE INITIAL DAM REGION:
!         RESET TEMPORAL WETTING COUNTER.

          WET_COUNT(I) = 0


        ENDIF


      ENDDO


!-----------------------------------------------------------------------
!
!     DIAGNOSTIC OUTPUT
!
!     REFERENCE TIMES:
!
!       1, 5, 10, 20, 30 s
!
!-----------------------------------------------------------------------


      IF((ATL.GE.1.D0).AND.(.NOT.MESSAGE_1)) THEN

        WRITE(LU,*)
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*) 'B-2 EROSION CHECK AT ABOUT 1 s'
        WRITE(LU,*) 'ATL                 = ',ATL
        WRITE(LU,*) 'ACTIVE WET NODES    = ',NACTIVE
        WRITE(LU,*) 'WET-ELIGIBLE NODES  = ',NWET_ELIGIBLE
        WRITE(LU,*) 'ERODED NODES        = ',NERODED
        WRITE(LU,*) 'DZ-LIMITED NODES    = ',NLIMITED
        WRITE(LU,*) 'MAX TAUB [Pa]       = ',MAX_TAUB
        WRITE(LU,*) 'MAX SOLID E [m/s]   = ',
     &               MAX_SOLID_EROSION_RATE
        WRITE(LU,*) 'MAX BED RATE [m/s]  = ',
     &               MAX_BED_LOWERING_RATE
        WRITE(LU,*) 'MAX RAW DZ [m]      = ',MAX_DZ_RAW
        WRITE(LU,*) 'MAX DZ/STEP [m]     = ',MAX_DZ
        WRITE(LU,*) 'SUM DZ/STEP [m]     = ',TOTAL_DZ
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*)

        MESSAGE_1 = .TRUE.

      ENDIF


      IF((ATL.GE.5.D0).AND.(.NOT.MESSAGE_5)) THEN

        WRITE(LU,*)
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*) 'B-2 EROSION CHECK AT ABOUT 5 s'
        WRITE(LU,*) 'ATL                 = ',ATL
        WRITE(LU,*) 'ACTIVE WET NODES    = ',NACTIVE
        WRITE(LU,*) 'WET-ELIGIBLE NODES  = ',NWET_ELIGIBLE
        WRITE(LU,*) 'ERODED NODES        = ',NERODED
        WRITE(LU,*) 'DZ-LIMITED NODES    = ',NLIMITED
        WRITE(LU,*) 'MAX TAUB [Pa]       = ',MAX_TAUB
        WRITE(LU,*) 'MAX SOLID E [m/s]   = ',
     &               MAX_SOLID_EROSION_RATE
        WRITE(LU,*) 'MAX BED RATE [m/s]  = ',
     &               MAX_BED_LOWERING_RATE
        WRITE(LU,*) 'MAX RAW DZ [m]      = ',MAX_DZ_RAW
        WRITE(LU,*) 'MAX DZ/STEP [m]     = ',MAX_DZ
        WRITE(LU,*) 'SUM DZ/STEP [m]     = ',TOTAL_DZ
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*)

        MESSAGE_5 = .TRUE.

      ENDIF


      IF((ATL.GE.10.D0).AND.(.NOT.MESSAGE_10)) THEN

        WRITE(LU,*)
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*) 'B-2 EROSION CHECK AT ABOUT 10 s'
        WRITE(LU,*) 'ATL                 = ',ATL
        WRITE(LU,*) 'ACTIVE WET NODES    = ',NACTIVE
        WRITE(LU,*) 'WET-ELIGIBLE NODES  = ',NWET_ELIGIBLE
        WRITE(LU,*) 'ERODED NODES        = ',NERODED
        WRITE(LU,*) 'DZ-LIMITED NODES    = ',NLIMITED
        WRITE(LU,*) 'MAX TAUB [Pa]       = ',MAX_TAUB
        WRITE(LU,*) 'MAX SOLID E [m/s]   = ',
     &               MAX_SOLID_EROSION_RATE
        WRITE(LU,*) 'MAX BED RATE [m/s]  = ',
     &               MAX_BED_LOWERING_RATE
        WRITE(LU,*) 'MAX RAW DZ [m]      = ',MAX_DZ_RAW
        WRITE(LU,*) 'MAX DZ/STEP [m]     = ',MAX_DZ
        WRITE(LU,*) 'SUM DZ/STEP [m]     = ',TOTAL_DZ
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*)

        MESSAGE_10 = .TRUE.

      ENDIF


      IF((ATL.GE.20.D0).AND.(.NOT.MESSAGE_20)) THEN

        WRITE(LU,*)
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*) 'B-2 EROSION CHECK AT ABOUT 20 s'
        WRITE(LU,*) 'ATL                 = ',ATL
        WRITE(LU,*) 'ACTIVE WET NODES    = ',NACTIVE
        WRITE(LU,*) 'WET-ELIGIBLE NODES  = ',NWET_ELIGIBLE
        WRITE(LU,*) 'ERODED NODES        = ',NERODED
        WRITE(LU,*) 'DZ-LIMITED NODES    = ',NLIMITED
        WRITE(LU,*) 'MAX TAUB [Pa]       = ',MAX_TAUB
        WRITE(LU,*) 'MAX SOLID E [m/s]   = ',
     &               MAX_SOLID_EROSION_RATE
        WRITE(LU,*) 'MAX BED RATE [m/s]  = ',
     &               MAX_BED_LOWERING_RATE
        WRITE(LU,*) 'MAX RAW DZ [m]      = ',MAX_DZ_RAW
        WRITE(LU,*) 'MAX DZ/STEP [m]     = ',MAX_DZ
        WRITE(LU,*) 'SUM DZ/STEP [m]     = ',TOTAL_DZ
        WRITE(LU,*) '-----------------------------------------------'
        WRITE(LU,*)

        MESSAGE_20 = .TRUE.

      ENDIF


      IF((ATL.GE.30.D0).AND.(.NOT.MESSAGE_30)) THEN

        WRITE(LU,*)
        WRITE(LU,*) '==============================================='
        WRITE(LU,*) 'B-2 EROSION CHECK AT ABOUT 30 s'
        WRITE(LU,*) 'ATL                 = ',ATL
        WRITE(LU,*) 'ACTIVE WET NODES    = ',NACTIVE
        WRITE(LU,*) 'WET-ELIGIBLE NODES  = ',NWET_ELIGIBLE
        WRITE(LU,*) 'ERODED NODES        = ',NERODED
        WRITE(LU,*) 'DZ-LIMITED NODES    = ',NLIMITED
        WRITE(LU,*) 'MAX TAUB [Pa]       = ',MAX_TAUB
        WRITE(LU,*) 'MAX SOLID E [m/s]   = ',
     &               MAX_SOLID_EROSION_RATE
        WRITE(LU,*) 'MAX BED RATE [m/s]  = ',
     &               MAX_BED_LOWERING_RATE
        WRITE(LU,*) 'MAX RAW DZ [m]      = ',MAX_DZ_RAW
        WRITE(LU,*) 'MAX DZ/STEP [m]     = ',MAX_DZ
        WRITE(LU,*) 'SUM DZ/STEP [m]     = ',TOTAL_DZ
        WRITE(LU,*) '==============================================='
        WRITE(LU,*)

        MESSAGE_30 = .TRUE.

      ENDIF


      RETURN

      END