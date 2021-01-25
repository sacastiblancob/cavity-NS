!                     **********************
                      SUBROUTINE CSC_MMATVEC
!                     **********************
     & (MA,BV,CV,NE)
!
!***********************************************************************
! 2D-DIFFUSION SOLVER SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) COMPUTES THE MATRIX VECTOR MULTIPLICATION Ab = c
!
!history  Sergio Castiblanco
!+        18/01/2021
!+        Translation for original Matlab implementation
!
!+
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| MA        |-->| MATRIX A IN CSC STORAGE                             |
!| BV        |-->| VECTOR FOR WHICH A IT'S BEEN MULTIPLIED             |
!| CV        |<->| SOLUTION OF THE MULTIPLICATION                      |
!| NE        |-->| NUMBER OF ELEMENTS OF BV, CV AND MA COLUMNS         |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_MMATVEC => CSC_MMATVEC
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE(CSC_OBJ), INTENT(IN) :: MA
      INTEGER, INTENT(IN) :: NE
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NE) :: BV
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(NE) :: CV
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: I,J
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!  CHECKING DIMENSIONS
!
      IF(NE.NE.MA%NC) THEN
        WRITE(*,*) '!!!MMATVEC_ERROR!!!: DIMENSION OF B ',NE,
     &    ' DOES NOT AGREE WITH DIMENSIONS OF MATRIX',MA%NAM,': ',MA%NC
        STOP 1
      ENDIF
!
!  COMPUTING THE MULTIPLICATION
!
      CV = 0D0
      DO J = 1,MA%NC
        DO I = MA%C(J),(MA%C(J+1)-1)
          CV(MA%R(I)) = CV(MA%R(I)) + MA%V(I)*BV(J)
        ENDDO
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_MMATVEC

