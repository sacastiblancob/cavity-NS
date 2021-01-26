!                     *****************
                      SUBROUTINE CSC_CG
!                     *****************
     & (MA,BV,XV,MAXNITER,NITER,TOL,NE)
!
!***********************************************************************
! CSC_TOOLS - CONJUGATE GRADIENT SOLVER
!***********************************************************************
!
!brief    1) SOLVES SYSTEM Ax = b WITH CONJUGATE GRADIENT METHOD
!
!history  Sergio Castiblanco
!+        18/01/2021
!+        Translation for original Matlab implementation
!
!+
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| MA        |-->| MATRIX A IN CSC STORAGE                             |
!| BV        |-->| RIGHT HAND SIDE VECTOR B                            |
!| XV        |<->| INPUT WITH FIRST VALUES, AND OUTPUT WITH SOLUTION   |
!| MAXNITER  |-->| MAXIMUM NUMBER OF ITERATIONS                        |
!| NITER     |<--| NUMBER OF ITERATIONS                                |
!| TOL       |-->| SOLVER TOLERANCE OF THE RESIDUAL                    |
!| NE        |-->| NUMBER OF ELEMENTS OF BV, XV AND MA ROWS            |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_CG => CSC_CG
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE(CSC_OBJ), INTENT(IN) :: MA
      INTEGER, INTENT(IN) :: NE
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NE) :: BV
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NE) :: XV
      INTEGER, INTENT(IN) :: MAXNITER
      INTEGER, INTENT(OUT) :: NITER
      DOUBLE PRECISION, INTENT(IN) :: TOL
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: IT
      DOUBLE PRECISION, DIMENSION(NE) :: DUMV,RO,D,AD,R
      DOUBLE PRECISION :: ALF, BET
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  CHECKING DIMENSIONS
!
      IF(NE.NE.MA%NC) THEN
        WRITE(*,*) '!!!CSC_CG_ERROR!!!: DIMENSION OF B ',NE,
     &    ' DOES NOT AGREE WITH DIMENSIONS OF MATRIX',MA%NAM,': ',MA%NC
        STOP 1
      ENDIF
!
!  COMPUTING CG
!
      !INITIAL RO
      RO = 0D0
      CALL CSC_MMATVEC(MA,XV,DUMV,NE)
      RO = BV - DUMV
      !D
      D = RO
      DO IT = 1,MAXNITER
        !AD
        CALL CSC_MMATVEC(MA,D,AD,NE)
        !ALF
        ALF = (DOT_PRODUCT(RO,RO))/(DOT_PRODUCT(D,AD))
        !UPDATING XV
        XV = XV + (ALF*D)
        !CHECKING CONVERGENCE
        CALL CSC_MMATVEC(MA,XV,DUMV,NE)
        IF((NORM2(BV - DUMV)).LT.TOL) THEN
          EXIT
        ENDIF
        !R
        R = RO - (ALF*AD)
        !BET
        BET = (DOT_PRODUCT(R,R))/(DOT_PRODUCT(RO,RO))
        !UPDATING D
        D = R + (BET*D)
        !UPDATING RO
        RO = R
      ENDDO
!
!  NITER
      NITER = IT - 1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IF(NITER.EQ.MAXNITER) THEN
        WRITE(*,*) 'CSC_CG_WARNING!!: MAXNITER REACHED!!!'
      ENDIF
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_CG












