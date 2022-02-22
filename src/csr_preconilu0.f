!                   *************************
                    SUBROUTINE CSR_PRECONILU0
!                   *************************
     & (LUV,MA)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) COMPUTES THE INCOMPLETE-LU(0) DECOMPOSITION OF MA AND
!            RETURN THE VALUES IN LUV, CSR STORAGE
!
!history  Sergio Castiblanco
!+        21/02/2022
!+        Translation for original Matlab implementation
!
!+
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| LUV      |<->| VECTOR OF VALUES OF LU MATRIX IN CSC                 |
!| MA       |-->| MATRIX A IN CSR                                      |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSR_PRECONILU0 => CSR_PRECONILU0
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE(CSC_OBJ), INTENT(IN) :: MA
      DOUBLE PRECISION, DIMENSION(MA%NZ), INTENT(INOUT) :: LUV
!
!  IN SUBROUTINE VARIABLES
!
      DOUBLE PRECISION, DIMENSION(MA%NZ) :: AV
      INTEGER, DIMENSION(MA%NZ) :: AC
      INTEGER, DIMENSION(MA%NC+1) :: AR
      INTEGER, DIMENSION(MA%NC) :: DIA, POINT
      INTEGER :: N,I,J,K,W,V
!
! COMPUTING ILU0 FACTORIZATION (MITTAL & AL-KURDI, 2003)
!
      AV = MA%V
      AR = MA%C
      AC = MA%R
      LUV = AV
!
      N = SIZE(AR)-1
      DIA = 0
      DO I=1,N
        DO J=AR(I),(AR(I+1)-AR(1))
          IF(AC(J).EQ.I)THEN
            DIA(I)=J
          ENDIF
        ENDDO
      ENDDO
      POINT = 0
!
      DO I=2,N
        DO V=AR(I)+1,(AR(I+1)-AR(1))
          POINT(AC(V))=V
        ENDDO
        DO V=AR(I),DIA(I)-1
          J=AC(V)
          LUV(V) = LUV(V)/LUV(DIA(J))
          DO W = (DIA(J)+1),(AR(J+1)-AR(1))
            K = POINT(AC(W))
            IF(K.GT.0)THEN
              LUV(K) = LUV(K) - LUV(V)*LUV(W)
            !ELSE   !IF REVISED, IF NOT (DIAGONAL MAY HAVE ZEROS)
            !  LUV(DIA(I)) = LUV(DIA(I)) - LUV(V)*LUV(W)
            ENDIF
          ENDDO
        ENDDO
        DO V=AR(I)+1,(AR(I+1)-AR(1))
          POINT(AC(V))=0
        ENDDO
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSR_PRECONILU0
