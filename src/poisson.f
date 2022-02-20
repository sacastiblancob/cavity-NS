!                       ******************
                        SUBROUTINE POISSON
!                       ******************
     & (LM,DUPDX,DVPDY,UL,SOLP,PCP,KSSP,LUVP,DLM,LPM,LQM,MNITERS,NITS,
     &     TOLSOR,TIEP,P,PP,FRESP)
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
!+        20/02/2022
!+        Update: adding new solvers options
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| LM        |-->| POISSON SYSTEM OF EQUATIONS MATRIX (L)              |
!| DUPDX     |-->| D(UP)/DX VECTOR                                     |
!| DVPDY     |-->| D(VP)/DY VECTOR                                     |
!| UL        |-->| SINGULAR VECTOR FOR REGULARIZATION PROCEDURE        |
!| SOLP      |-->| SOLVER FOR POISSON STEP                             |
!| PCP       |-->| PRECONDITIONING OPTION FOR POISSON STEP             |
!| KSSP      |-->| KRYLOV SUBSPACE SIZE FOR POISSON STEP, IF GMRES     |
!| LUVP      |-->| PACKED-LU VALUES IF SSOR OR ILU(0) PRECON.          |
!| DLM       |-->| DIAGONAL ENTRIES OF LM FOR (ABS)DIAG. PRECON.       |
!| LPM       |-->| P PART OF L, IF SOR SOLVER                          |
!| LQM       |-->| Q PART OF L, IF SOR SOLVER                          |
!| MNITERS   |-->| MAXIMUM NUMBER OF ITERATIONS FOR POISSON SOLVER     |
!| NITS      |<--| NUMBER OF ITERATIONS TAKEN BY THE SOLVER            |
!| TOLSOR    |-->| TOLERANCE FOR THE RESIDUAL OF POISSON SOLVER        |
!| TIEP      |-->| VALUE TO TIE PRESURE, IN THE BOTTOM-LEFT CORNER     |
!| P         |<--| PRESURE VECTOR                                      |
!| PP        |<--| PRESURE FOR PLOT VECTOR                             |
!| FRESP     |<--| FINAL RESIDUAL OBTAINED BY SOLVER                   |
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
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(IN) :: UL
      TYPE(CSC_OBJ), INTENT(IN) :: LM, LPM, LQM
      DOUBLE PRECISION, INTENT(IN) :: TOLSOR, TIEP
      DOUBLE PRECISION, INTENT(OUT) :: FRESP
      INTEGER, INTENT(IN) :: SOLP,PCP,KSSP,MNITERS
      INTEGER, INTENT(OUT) :: NITS
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(INOUT) :: P
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(OUT) :: PP
      DOUBLE PRECISION, DIMENSION(LM%NC), INTENT(IN) :: DLM
      DOUBLE PRECISION, DIMENSION(LM%NZ), INTENT(IN) :: LUVP
!
!  IN SUBROUTINE VARIABLES
!
      DOUBLE PRECISION, DIMENSION(NX*NY) :: RHS
      DOUBLE PRECISION :: RHTU
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!  COMPUTING RIGHT HAND SIDE
!
      IF(DY.GE.DX) THEN
        RHS = ((DX**2)*RHO/DT)*(DUPDX + DVPDY)
      ELSE
        RHS = ((DY**2)*RHO/DT)*(DUPDX + DVPDY)
      ENDIF
!
!  REGULARIZATION OF RHS
!
      RHTU = DOT_PRODUCT(RHS,UL)
      RHS = RHS - RHTU*UL
!
!  SOLVING POISSON'S EQUATION
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO SOR SOLVER'
!
      IF(SOLP.EQ.1)THEN        !SOR SOLVER      
        CALL CSC_SOR(LM,LPM,LQM,RHS,P,MNITERS,NITS,TOLSOR,NX*NY)
!      
      ELSE IF(SOLP.EQ.2)THEN   !BICGSTAB SOLVER
        CALL CSC_BICGSTAB(LM,RHS,P,MNITERS,TOLSOR,LUVP,DLM,PCP,NITS,
     &                     FRESP)
!
      ELSE IF(SOLP.EQ.3)THEN   !GMRES SOLVER
        CALL CSC_GMRES(LM,RHS,P,KSSP,MNITERS,TOLSOR,LUVP,DLM,PCP,NITS,
     &                  FRESP)
!
      ELSE
        WRITE(*,*) "INVALID OPTION FOR POISSON SOLVER: ABORTING!!!"
        STOP 0
      ENDIF
!
      IF(DEBUG) WRITE(*,*) 'OUT SOR SOLVER, NITS: ', NITS
!
!  COMPUTING PLOTTING PURPOSE PRESURE
!
      PP = P - (P(1)-TIEP)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE POISSON
