!                       ******************
                        SUBROUTINE POISSON
!                       ******************
     & (DUPDX,DVPDY,RM,LM,LPM,LQM,MNITERS,NITS,TOLSOR,TIEP,P,PP)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!***********************************************************************
! NAVIER-STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) COMPUTING REGULARIZATION AND POISSON EQUATION
!
!history  Sergio Castiblanco
!+        26/01/2021
!+        Translation for original Matlab implementation
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| DUPDX     |-->| D(UP)/DX VECTOR                                     |
!| DVPDY     |-->| D(VP)/DY VECTOR                                     |
!| RM        |-->| REGULARIZATION MATRIX                               |
!| LM        |-->| POISSON SYSTEM OF EQUATIONS MATRIX (L)              |
!| LPM       |-->| P PART OF L, FOR SOR SOLVER                         |
!| LQM       |-->| Q PART OF L, FOR SOR SOLVER                         |
!| MNITERS   |-->| MAXIMUM NUMBER OF ITERATIONS FOR SOR SOLVER         |
!| NITS      |<->| NUMBER OF ITERATIONS TAKEN BY THE SOLVER            |
!| TOLSOR    |-->| TOLERANCE FOR THE RESIDUAL OF SOR SOLVER            |
!| TIEP      |-->| VALUE TO TIE PRESURE, IN THE BOTTOM-LEFT CORNER     |
!| P         |<--| PRESURE VECTOR                                      |
!| PP        |<--| PRESURE FOR PLOT VECTOR                             |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE
      USE DECLARATIONS_NUMERICAL, ONLY: DX,DY,NX,NY,DT,DEBUG
      USE DECLARATIONS_PHYSIC
!
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(IN) :: DUPDX, DVPDY
      DOUBLE PRECISION, DIMENSION(NX*NY,NX*NY), INTENT(IN) :: RM
      TYPE(CSC_OBJ), INTENT(IN) :: LM, LPM, LQM
      DOUBLE PRECISION, INTENT(IN) :: TOLSOR, TIEP
      INTEGER, INTENT(IN) :: MNITERS
      INTEGER, INTENT(OUT) :: NITS
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(INOUT) :: P
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(OUT) :: PP
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER I
      DOUBLE PRECISION, DIMENSION(NX*NY) :: RHS, RHSR
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!  COMPUTING RIGHT HAND SIDE
!
      IF(DY.GE.DX) THEN
        DO I = 1,NX*NY
          RHS(I) = ((DX**2)*RHO/DT)*(DUPDX(I) + DVPDY(I))
        ENDDO
      ELSE
        DO I = 1,NX*NY
          RHS(I) = ((DY**2)*RHO/DT)*(DUPDX(I) + DVPDY(I))
        ENDDO
      ENDIF
!
!  REGULARIZATION OF RHS
!
      RHSR = MATMUL(RM,RHS)
!
!  SOLVING POISSON'S EQUATION
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO SOR SOLVER'
      CALL CSC_SOR(LM,LPM,LQM,RHSR,P,MNITERS,NITS,TOLSOR,NX*NY)
      IF(DEBUG) WRITE(*,*) 'OUT SOR SOLVER, NITS: ', NITS
!
!  COMPUTING PLOTTING PURPOSE PRESURE
!
      PP = P - (P(1)-TIEP)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE POISSON



