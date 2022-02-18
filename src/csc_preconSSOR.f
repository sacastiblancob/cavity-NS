!                   *************************
                    SUBROUTINE CSC_PRECONSSOR
!                   *************************
     & (PQV,MA,W)
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
!
! This function takes the matrix A in CSC storage, and prepares the
! matrices for solve the system Ax=b with SSOR preconditioner
!
! P and Q:
! The SSOR preconditioner is defined with the next matrix:
! 
!               Mssor = (D - w*E) D^(-1) (D - wF)
!
! In which the system can be solved throug the matrices:
!
!               P = (D- w*E)D^(-1) = (I - w*E*D^(-1))
!               Q = (D - w*F)
!
! Then, for finding z = Mssor^(-1)*x
!
!               solve Py = x;  -> Forward substitution
!               solve Qz = y;  -> Backward substitution
!
! For the particular case when w = 1, we are in Symmetric Gauss-Seidel
! preconditioner.
! Note that diagonal of matrix P = 1. Then this decomposition on P and Q
! can be stored in packed format,i.e., use the same structure as A, but
! storing in the lower part -w*E*D^(-1), and in the upper part Q. Which
! is wath PQV stores.
!
! Entries:
!     MA : Matrix A in CSC stotage
!     W : SOR Over-Relaxation coefficient
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| PQV      |<->| VECTOR OF VALUES OF SSOOR PRECON IN CSC              |
!| MA       |-->| MATRIX A IN CSC                                      |
!| W        |-->| SOR COEFFICIENT, IF W==1 SYMMETRIC GAUSS-SEIDEL      |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_PRECONSSOR => CSC_PRECONSSOR
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE(CSC_OBJ), INTENT(IN) :: MA
      DOUBLE PRECISION, DIMENSION(MA%NZ), INTENT(INOUT) :: PQV
      DOUBLE PRECISION, INTENT(IN) :: W
!
!  IN SUBROUTINE VARIABLES
!
      DOUBLE PRECISION, DIMENSION(MA%NC) :: D
      INTEGER ::M,I,J,PQP
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      M = SIZE(MA%C)-1
!
! GETTING DIAGONAL OF A
      D = 0.0D0
      DO J=1,M
        DO I=MA%C(J),MA%C(J+1)-MA%C(1)
          IF(MA%R(I).EQ.J)THEN
            D(J) = MA%V(I)
          ENDIF
        ENDDO
      ENDDO
!
      PQV = 0.0D0
!
! GETTING SSOR DECOMPOSITION
      PQP=1
      DO J=1,M
        DO I=MA%C(J),MA%C(J+1)-MA%C(1)
          IF(MA%R(I).LT.J)THEN
            PQV(PQP) = W*MA%V(I)
          ELSEIF(MA%R(I).GT.J)THEN
            PQV(PQP) = W*MA%V(I)/D(J)
          ELSE
            PQV(PQP) = MA%V(I)
          ENDIF
          PQP = PQP + 1
        ENDDO
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_PRECONSSOR
