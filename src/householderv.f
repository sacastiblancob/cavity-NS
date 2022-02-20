!                   ***********************
                    SUBROUTINE HOUSEHOLDERV
!                   ***********************
     & (VX,BETA,V,N)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) COMPUTES THE ENTRIES OF SSOR PRECONDITIONER AND RETURN THEM
!            IN PQV, CSC STORAGE 
!
!history  Sergio Castiblanco
!+        17/02/2022
!+        Translation for original Matlab implementation
!
!+
! This function computes the vector v for constructing Householder
! Projector related to vector x, i.e
!           P = I - beta*v*v'
! Then,          Px = ||x||*e_1
!
! Algorithm 5.1.1 (Householder Vector) Matrix Computations
! (Golub & Van Loan)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| X        |-->| VECTOR FOR COMPUTING HOUSEHOLDER MIRROR              |
!| BETA     |<--| BETA COEFFICIENT                                     |
!| V        |<->| HOUSEHOLDER MIRROR VECTOR                            |
!| N        |-->| DIMENSSIONS OF VECTORS X AND V                       |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_HOUSEHOLDERV => HOUSEHOLDERV
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      INTEGER, INTENT(IN) :: N
      DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: VX
      DOUBLE PRECISION, DIMENSION(N), INTENT(INOUT) :: V
      DOUBLE PRECISION, INTENT(OUT) :: BETA
!
!  IN SUBROUTINE VARIABLES
!
      DOUBLE PRECISION, DIMENSION(N) :: X
      DOUBLE PRECISION :: SIG, MU, NX
      INTEGER :: I
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! NORMALIZING TO AVOID OVERFLOW
      NX = NORM2(VX)
      IF(NX.EQ.0.0D0)THEN
        WRITE(*,*) "ERROR HOUSEHOLDERV, NORM(X)==0"
        STOP 0
      ENDIF
      X = VX/NX
!
      SIG = DOT_PRODUCT(X(2:N),X(2:N))
      V = 0.D0
      V(1) = 1.D0
      DO I=2,N
        V(I)=X(I)
      ENDDO
      IF(N.EQ.1)THEN
        BETA = 2
        V = 1.0D0
      ELSE
        IF(SIG.EQ.0)THEN
          BETA = 0
        ELSE
          MU = SQRT(X(1)**2 + SIG)
          IF(X(1).LE.0)THEN
            V(1) = X(1) - MU
          ELSE
            V(1) = -SIG/(X(1)+MU)
          ENDIF
          BETA = (2*V(1)**2)/(SIG+V(1)**2)
          DO I=1,N
            V(I)=V(I)/V(1)
          ENDDO
        ENDIF
      ENDIF      
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE HOUSEHOLDERV
