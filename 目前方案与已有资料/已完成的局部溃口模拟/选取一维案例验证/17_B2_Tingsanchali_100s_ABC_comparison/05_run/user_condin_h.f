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
!     B2 - TINGSANCHALI & CHINNARASRI (2001), RUN B-2
!
!     CORRECTED FULL-FLUME FIXED-BED HYDRODYNAMIC STARTUP
!
!     PHYSICAL FLUME:
!
!          X = -17.50 TO 17.50 m
!
!          Y =   0.00 TO  1.00 m
!
!
!     DAM:
!
!          UPSTREAM TOE   = -2.40 m
!
!          CREST FRONT    =  0.00 m
!
!          CREST BACK     =  0.30 m
!
!          DOWNSTREAM TOE =  4.30 m
!
!          CREST ELEVATION = 0.80 m
!
!
!     RUN B-2 INITIAL CONDITION:
!
!          INITIAL HEAD ABOVE CREST = 0.05 m
!
!          INITIAL UPSTREAM FREE SURFACE = 0.85 m
!
!
!     COMPUTATIONAL INITIAL CONDITION:
!
!          X <= 0:
!
!              H = MAX(0.85 - ZF, 0)
!
!          X > 0:
!
!              H = 0
!
!
!     THEREFORE THE FULL UPSTREAM PART OF THE 35 m FLUME IS
!     INITIALISED AS THE RESERVOIR, WHILE THE DOWNSTREAM PART
!     IS INITIALLY DRY.
!
!
!     INITIAL VELOCITY:
!
!          U = 0
!
!          V = 0
!
!
!     IMPORTANT:
!
!          THIS FILE IS FOR THE FIXED-BED HYDRODYNAMIC TEST ONLY.
!
!          NO EROSION IS ACTIVE.
!
!          NO DYNAMIC-BOTTOM ROUTINE IS CALLED.
!
!          NO USER_UTIMP_TELEMAC2D IS DEFINED.
!
!-----------------------------------------------------------------------


      ETA0 = 0.85D0


      DO I = 1,NPOIN


        IF(X(I).LE.0.D0) THEN

          H%R(I)=ETA0-ZF%R(I)

          IF(H%R(I).LT.0.D0) THEN

            H%R(I)=0.D0

          ENDIF


        ELSE

          H%R(I)=0.D0

        ENDIF


      ENDDO


!-----------------------------------------------------------------------
!
!     INITIAL VELOCITY
!
!-----------------------------------------------------------------------

      CALL OS('X=0     ',X=U)

      CALL OS('X=0     ',X=V)


!-----------------------------------------------------------------------
!
!     DIAGNOSTIC MESSAGE
!
!-----------------------------------------------------------------------

      WRITE(LU,*)
      WRITE(LU,*) '==============================================='
      WRITE(LU,*) 'B2 FULL-FLUME FIXED-BED INITIAL CONDITION'
      WRITE(LU,*) 'FLUME X                = -17.50 TO 17.50 m'
      WRITE(LU,*) 'UPSTREAM FREE SURFACE  = ',ETA0
      WRITE(LU,*) 'RESERVOIR REGION       = X <= 0'
      WRITE(LU,*) 'DOWNSTREAM INITIAL     = DRY'
      WRITE(LU,*) 'INITIAL U,V            = 0'
      WRITE(LU,*) 'NO EROSION IS ACTIVE'
      WRITE(LU,*) '==============================================='
      WRITE(LU,*)


      RETURN

      END