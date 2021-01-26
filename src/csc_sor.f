!                     ******************
                      SUBROUTINE CSC_SOR
!                     ******************
     & (MA,MP,MQ,BV,XV,MAXNITER,NITER,TOL,NE)
!
!***********************************************************************
! CSC_TOOLS - SUCCESIVE OVER-RELAXATION SOLVER
!***********************************************************************
!
!brief    1) SOLVES SYSTEM Ax = b WITH SUCCESIVE OVER-RELAXATION
!
!history  Sergio Castiblanco
!+        26/01/2021
!+        Translation for original Matlab implementation
!
!!
!!
!! P and Q:
!! Classic iterative methods use the next recurrence equation for solve 
!! the system:
!! 
!!               P * x(t+1) = Q * x(t) + b
!!
!!       For SOR:
!!            P = ((1/w)*D + L)
!!            Q = ((1/w - 1)*D - U); where: L = lower part of A
!!                                          U = upper part of A
!!                                          D = diagonal of A
!!
!! For the particular case when w = 1, we are in Gauss-Seidel method
!!
!!      Sergio A. Castiblanco B. - Advanced Numerical Methods
!!      Pontificia Universidad Javeriana - Bogota
!!
!+
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| MA        |-->| MATRIX A IN CSC STORAGE                             |
!| MP        |-->| MATRIX P IN CSC STORAGE                             |
!| MQ        |-->| MATRIX Q IN CSC STORAGE                             |
!| BV        |-->| RIGHT HAND SIDE VECTOR B                            |
!| XV        |<->| INPUT WITH FIRST VALUES, AND OUTPUT WITH SOLUTION   |
!| MAXNITER  |-->| MAXIMUM NUMBER OF ITERATIONS                        |
!| NITER     |<--| NUMBER OF ITERATIONS                                |
!| TOL       |-->| SOLVER TOLERANCE OF THE RESIDUAL                    |
!| NE        |-->| NUMBER OF ELEMENTS OF BV, XV AND MA ROWS            |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_SOR => CSC_SOR
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE(CSC_OBJ), INTENT(IN) :: MA,MP,MQ
      INTEGER, INTENT(IN) :: NE
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NE) :: BV
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NE) :: XV
      INTEGER, INTENT(IN) :: MAXNITER
      INTEGER, INTENT(OUT) :: NITER
      DOUBLE PRECISION, INTENT(IN) :: TOL
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: IT,J,I
      DOUBLE PRECISION, DIMENSION(NE) :: DUMV,DXV
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DO IT = 1,MAXNITER
        CALL CSC_MMATVEC(MQ,XV,DXV,SIZE(XV))
        XV = DXV + BV
        DO J = 1,MA%NC
          DO I = MP%C(J),MP%C(J+1)-MP%C(1)
            IF(J.EQ.MP%R(I)) THEN
              XV(J) = XV(J)/MP%V(I)
            ELSE
              XV(MP%R(I)) = XV(MP%R(I)) - MP%V(I)*XV(J)
            ENDIF
          ENDDO
        ENDDO
        CALL CSC_MMATVEC(MA,XV,DUMV,SIZE(XV))
        IF(NORM2(BV-DUMV).LT.TOL) THEN
          EXIT
        ENDIF
      ENDDO
!
      NITER = IT-1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_SOR




