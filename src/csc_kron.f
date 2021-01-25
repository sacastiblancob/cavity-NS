!                     *******************
                      SUBROUTINE CSC_KRON
!                     *******************
     & (MA,MB,MC,NAMC)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) COMPUTES THE KRONECKER PRODUCT BETWEEN TO MATRICES IN CSC
!
!history  Sergio Castiblanco
!+        17/01/2021
!+        Translation for original Matlab implementation
!
!+
!!
!! This function computes the kronecker product between A and B with 
!! A and B stored in CSC
!!
!! Example
!!          A               B
!!      |1  1  0|        |1 1 0|
!!      |2  2  2|        |0 1 0|
!!      |0  3  3|        |0 1 1|
!!
!! Output
!!
!!                     |1 1   1 1        |
!!                     |  1     1        |
!!          C          |  1 1   1 1      |
!!      |1B  1B  0 |   |2 2   2 2   2 2  |
!!      |2B  2B  2B| = |  2     2     2  |
!!      |0   3B  3B|   |  2 2   2 2   2 2|
!!                     |      3 3   3 3  |
!!                     |        3     3  |
!!                     |        3 3   3 3|
!!
!! The output, as the entries, is stored in CSC
!!
!!
!!      Sergio A. Castiblanco B. - Advanced Numerical Methods
!!      Pontificia Universidad Javeriana - Bogota
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| MA       |-->| MATRIX A IN CSC                                      |
!| MB       |-->| MATRIX B IN CSC                                      |
!| MC       |<->| MATRIX C IN CSC                                      |
!| NAMC     |-->| NAME OF MATRIX C                                     |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_KRON => CSC_KRON
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE(CSC_OBJ), INTENT(IN) :: MA,MB
      TYPE(CSC_OBJ), INTENT(INOUT) :: MC
      CHARACTER(LEN=6), INTENT(IN) :: NAMC
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: AJ, BJ, P, AF, BF, IA, IB, J
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
! ALLOCATING MC
      CALL ALL_CSC(MC, MA%NR*MB%NR, MA%NC*MB%NC, MA%NZ*MB%NZ,NAMC)
!
! WRITTING MC
      MC%C(1) = 1
      AJ = 1
      BJ = 1
      P = 1
      DO J=1,MC%NC
        AF = MA%C(AJ+1) - MA%C(AJ)
        BF = MB%C(BJ+1) - MB%C(BJ)
        MC%C(J+1) = MC%C(J) + AF*BF
        DO IA = MA%C(AJ),(MA%C(AJ+1)-1)
          DO IB = MB%C(BJ),(MB%C(BJ+1)-1)
            MC%R(P) = MB%NR*(MA%R(IA)-1) + MB%R(IB)
            MC%V(P) = MA%V(IA)*MB%V(IB)
            P = P+1
          ENDDO
        ENDDO
        BJ = BJ+1
        IF(BJ.EQ.(MB%NC+1)) THEN
          AJ = AJ + 1
          BJ = 1
        ENDIF
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_KRON

