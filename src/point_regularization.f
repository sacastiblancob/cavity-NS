!                   *******************************
                    SUBROUTINE POINT_REGULARIZATION
!                   *******************************
     & (SING,TOLSING,MNITERM,W,ISUSERW,LM,LPM,LQM,UL,VL,RM)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!***********************************************************************
! NAVIER-STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) COMPUTING REGULARIZATION MATRIX AND PREPARING EVERYTING FOR
!            SOLVING THE PRESURE THROUGH REGULARIZATION AND SOR METHOD
!
!history  Sergio Castiblanco
!+        25/01/2021
!+        Translation for original Matlab implementation
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| SING      |-->| NUMBER OF ITERATIONS FOR INVERSE POWER METHOD       |
!| TOLSING   |-->| TOLERANCE FOR INVERSE POWER METHOD SOLVER           |
!| MNITERM   |-->| MAXIMUM NUMBER OF ITERATIONS FOR IPM SOLVER         |
!| W         |-->| SOR OVER RELAXATION COEFFICIENT                     |
!| ISUSERW   |-->| LOGICAL FOR SOR OVER RELAXATION COEFFICIENT         |
!| ML        |-->| LAPLACIAN MATRIX                                    |
!| MLP       |<->| P MATRIX FOR SOR SOLVER                             |
!| MLQ       |<->| Q MATRIX FOR SOR SOLVER                             |
!| UL        |<->| LHS SINGULAR VECTOR, EIGENVECTOR OF L'              |
!| UV        |<--| RHS SINGULAR VECTOR, EIGENVECTOR OF L               |
!| RM        |<--| REGULARIZATION MATRIX = (I - UL*UL')                |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE CSC_STORAGE
      USE DECLARATIONS_NUMERICAL, ONLY: NX,NY,DEBUG
!
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! INTENT VARIABLES
!
      INTEGER, INTENT(IN) :: SING, MNITERM
      DOUBLE PRECISION, INTENT(IN) :: TOLSING, W
      LOGICAL, INTENT(IN) :: ISUSERW
      TYPE(CSC_OBJ), INTENT(IN) :: LM
      TYPE(CSC_OBJ), INTENT(INOUT) :: LPM, LQM
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(INOUT) :: UL
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(OUT) :: VL
      DOUBLE PRECISION, DIMENSION(NX*NY,NX*NY), INTENT(OUT) :: RM
!
! IN SUBROUTINE VARIABLES
!
      INTEGER :: I,PIT
      !TRANSPOSE OF L
      TYPE(CSC_OBJ) :: LTM
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
! TRANSPOSING L
!
      CALL CSC_TRANS(LM,LTM,'LTM   ')




!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE POINT_REGULARIZATION
