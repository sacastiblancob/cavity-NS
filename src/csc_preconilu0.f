!                   *************************
                    SUBROUTINE CSC_PRECONILU0
!                   *************************
     & (LUV,MA)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) COMPUTES THE INCOMPLETE-LU(0) DECOMPOSITION OF MA AND
!            RETURN THE VALUES IN LUV, CSC STORAGE
!
!history  Sergio Castiblanco
!+        17/02/2022
!+        Translation for original Matlab implementation
!
!+
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| LUV      |<->| VECTOR OF VALUES OF LU MATRIX IN CSC                 |
!| MA       |-->| MATRIX A IN CSC                                      |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_PRECONILU0 => CSC_PRECONILU0
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
      INTEGER, DIMENSION(MA%NZ) :: AC, LUC
      INTEGER, DIMENSION(MA%NR+1) :: AR
      INTEGER, DIMENSION(MA%NC) :: DIA, POINT
      INTEGER, DIMENSION(MA%NC+1) :: LUR
      INTEGER :: N,NZ,MH,I,IAA,IAB,J,JP,K,W,V
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! COMPUTING A TRANSPOSE, AND STORE RESULT IN CSR FORMAT IN AV,AC,AR
      AV = 0.D0
      AC = 0
      AR = 0
!
      N = SIZE(MA%C)-1
      NZ = MA%C(N+1)-MA%C(1)
!
      MH = MAXVAL(MA%R)+1
!
      DO I=1,NZ
        J = MA%R(I)+2;
        IF(J.LE.MH)THEN
          AR(J)=AR(J)+1
        ENDIF
      ENDDO
!
      AR(1)=MA%C(1)
      AR(2)=MA%C(1)
!
      DO I=3,MH
        AR(I)=AR(I)+AR(I-1)
      ENDDO
!
      DO I=1,N
        IAA = MA%C(I)
        IAB = MA%C(I+1)-MA%C(1)
        IF(IAB.LT.IAA)THEN
          EXIT
        ENDIF
        DO JP=IAA,IAB
          J = MA%R(JP)+1
          K = AR(J)
          AC(K) = I
          AV(K) = MA%V(JP)
          AR(J) = K+1
        ENDDO
      ENDDO
!
! COMPUTING ILU0 FACTORIZATION (MITTAL & AL-KURDI, 2003)
!
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
            ELSE   !IF REVISED, IF NOT (DIAGONAL MAY HAVE ZEROS)
              LUV(DIA(I)) = LUV(DIA(I)) - LUV(V)*LUV(W)
            ENDIF
          ENDDO
        ENDDO
        DO V=AR(I)+1,(AR(I+1)-AR(1))
          POINT(AC(V))=0
        ENDDO
      ENDDO
!
! NOW ILU0 IS STORED IN LUV, BUT IN CSR FORMAT, TRANSPOSING IT
!
      N = SIZE(AR)-1
      NZ = AR(N+1)-1
!
      LUC = 0
      LUR = 0
      AV = LUV
      LUV = 0D0
!
      MH = MAXVAL(AC)+1
!
      DO I=1,NZ
        J = AC(I)+2;
        IF(J.LE.MH)THEN
          LUR(J)=LUR(J)+1
        ENDIF
      ENDDO
!
      LUR(1)=AR(1)
      LUR(2)=AR(1)
!
      DO I=3,MH
        LUR(I)=LUR(I)+LUR(I-1)
      ENDDO
!
      DO I=1,N
        IAA = AR(I)
        IAB = AR(I+1)-AR(1)
        IF(IAB.LT.IAA)THEN
          EXIT
        ENDIF
        DO JP=IAA,IAB
          J = AC(JP)+1
          K = LUR(J)
          LUC(K) = I
          LUV(K) = AV(JP)
          LUR(J) = K+1
        ENDDO
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_PRECONILU0
