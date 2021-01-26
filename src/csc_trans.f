!                     ********************
                      SUBROUTINE CSC_TRANS
!                     ********************
     & (MA,MT,NAMT)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) COMPUTES THE TRANSPOSE OF MATRIX A IN CSC
!
!history  Sergio Castiblanco
!+        25/01/2021
!+        Translation for original Matlab implementation
!
!+
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| MA       |-->| MATRIX A IN CSC                                      |
!| MT       |<->| TRANSPOSE OF MATRIX A IN CSC                         |
!| NAMT     |-->| NAME OF MATRIX MT                                     |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_TRANS => CSC_TRANS
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE(CSC_OBJ), INTENT(IN) :: MA
      TYPE(CSC_OBJ), INTENT(INOUT) :: MT
      CHARACTER(LEN=6), INTENT(IN) :: NAMT
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: I,J,K
      INTEGER, DIMENSION(MA%NZ) :: COLS
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  STORING COLUMNS INDICES OF A IN THE VECTOR COLS
      DO J = 1,MA%NC
        DO I = MA%C(J),MA%C(J+1)-1
          COLS(I) = J
        ENDDO
      ENDDO
!
! ALLOCATING MT
!
      CALL ALL_CSC(MT, MA%NC, MA%NR, MA%NZ,NAMT)
!
! FILLING MT
!
      MT%C(1) = 1
      K = 1
      DO J = 1,MT%NC
        DO I = 1,MA%NZ
          IF(MA%R(I).EQ.J) THEN
            MT%R(K) = COLS(I)
            MT%V(K) = MA%V(I)
            MT%C(J+1) = MT%C(J+1)+1
            K = K + 1
          ENDIF
        ENDDO
      ENDDO
!
! ACCUMULATION OVER MT%C
!
      DO J = 2,MT%NC+1
        MT%C(J) = MT%C(J) + MT%C(J-1)
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_TRANS



