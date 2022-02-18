!                   ************************
                    SUBROUTINE CSC_SOLPACKLU
!                   ************************
     & (N,M,NZ,LUV,LUR,LUC,VB,VX)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) THIS SUBROUTINE SOLVES THE SYSTEM Ax=b WHERE MATRIX A IS LU
!            DECOMPOSITION, i.e., A=LU AND STORED IN CSC-PACKED STORAGE 
!
!history  Sergio Castiblanco
!+        16/02/2022
!+        Translation for original Matlab implementation
!
!-----------------------------------------------------------------------
! CSC-PACKED STORAGE IS THE NEXT:
! As L is Unitary Lower Triangular, i.e., its diagonal are all ones,
! there is no need to store its diagonal, insteadi, in the upper part
! and diagonal the components of matrix U can be stored. Using just one
! matrix to store both L and U. This is done in CSC storage, and it is
! useful for preconditioners as SSOR or ILU, as those have same
! structure as original matrix.
!
!               solve -> Ly = b (forward substitution)
!               solve -> Ux = y (backward substitution)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| N        |-->| ROWS OF MATRIX LU                                    |
!| M        |-->| COLUMNS OF MATRIX LU                                 |
!| NZ       |-->| NUMBERS OF NONZERO ELEMNTS OF MATRIX LU              |
!| LUV      |-->| VALUES OF MATRIX LU IN CSC-PACKED STORAGE            |
!| LUR      |-->| ROW POINTERS OF MATRIX LU IN CSC-PACKED STORAGE      |
!| LUC      |-->| COLUMN POINTERS OF MATRIX LU IN CSC-PACKED STORAGE   |
!| VB       |-->| RIGHT HAND SIDE VECTOR b                             |
!| VX       |<->| SOLUTION VECTOR FOR THE SYSTEM                       |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_SOLPACKLU => CSC_SOLPACKLU
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      INTEGER, INTENT(IN) :: N,M,NZ
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NZ) :: LUV,LUR
      DOUBLE PRECISION, INTENT(IN), DIMENSION(M+1) :: LUC
      DOUBLE PRECISION, INTENT(IN), DIMENSION(N) :: VB
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(N) :: VX
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: J, I
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! REPLACING VALUES OF VX WITH THOSE OF VB
!
      DO I=1,N
        VX(I)=VB(I)
      ENDDO
!
! FORWARD SUBSTITUTION
!
      DO J=1,M
        DO I=LUC(J),(LUC(J+1)-LUC(1))
          IF(J.LT.LUR(I))THEN
            VX(LUR(I)) = VX(LUR(I)) - LUV(I)*VX(J)
          ENDDO
        ENDDO
      ENDDO
!
! BACKWARD SUBSTITUTION
!
      DO J=M,1,-1
        DO I = (LUC(J+1)-LUC(1)),LUC(J),-1
          IF(J.EQ.LUR(I))THEN
            VX(J) = VX(J)/LUV(I)
          ELSE IF(J.GT.LUR(I))THEN
            VX(LUR(I)) = VX(LUR(I)) - LUV(I)*VX(J)
          ENDIF
        ENDDO
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_SOLPACKLU
