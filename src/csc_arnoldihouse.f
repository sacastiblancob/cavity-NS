!                  ***************************
                   SUBROUTINE CSC_ARNOLDIHOUSE
!                  ***************************
     & (MA,VV,M,LUV,VD,PC,WH,VB)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) COMPUTES THE ARNOLDI METHOD USING HOUSEHOLDER PROJECTIONS 
!
!history  Sergio Castiblanco
!+        17/02/2022
!+        Translation for original Matlab implementation
!
!+
!
! Entries:
!   - MA : matrix A of dimensions n*n in CSC storage
!   - VX : vector for compute Krylov subspace = ...
!             span{v,Av,(A^2)v,...,(A^m-1)v}
!   - M : dimension of Krylov subspace
!   - LUV  : values of LU decomposition matrix (for SSOR or ILU(0)) in
!              CSC_packed storage, i.e., if pc==3 or pc==4
!   - VD : Vector with values of D or abs(D) if pc==1 or pc==2
!   - PC  : preconditioning type
!       - 0 : No preconditioning
!       - 1 : Diagonal preconditioning
!       - 2 : Absolute diagonal preconditioning
!       - 3 : Symmetric SOR (Symmetric Gauss-Seidel -> SSOR, w=1) or ILU(0)
!
! Output
!   - WH : Matrix with Hessemberg matrix in its upper+1 part, and 
!          W vectors in its low part
!   - VB : beta coefficients for computing Householder Projectors
!           P = I - beta*w*w'
! 
! WH:  (n+1)x(m+1)
!    | h_10   h_11   h_12    .    .    .       h_1m  |
!    | 1.0    h_21   h_22    .    .    .       h_2m  |
!    | w_21   1.0    h_32    .    .    .       h_3m  |
!    | w_31   w_32   1.0     .    .    .       h_4m  |
!    |  .      .      .      .                  .    |
!    |  .      .      .           .             .    |
!    |  .      .      .                h_mm-1  h_mm  |
!    | w_m-11 w_m-12 w_m-13  .    .    1.0     h_m+1m|
!    | w_m1   w_m2   w_m3    .    .    w_mm+1  1.0   |
!    |  .      .      .      .    .     .      .     |
!    |  .      .      .      .    .     .      .     |
!    | w_n1   w_n2   w_n3    .    .     .     w_nm   |
!
! Note that w_ii=1.0, and shall not be stored, but that get computations
! quite difficult, instead the m w_ii are stored, which increase by m the
! real optimal storage.
!
! WH(1,1) = beta for GMRES residual computation
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| MA       |-->| MATRIX A IN CSC STORAGE                              |
!| VV       |-->| VECTOR FOR KRYLOV SUBSPACE K(MA,VV)_M                |
!| M        |-->| KRYLOV SUBSPACE SIZE                                 |
!| LUV      |-->| LU DECOMPOSITION VALUES (FOR SSOR OR ILU0)           |
!| VD       |-->| VECTOR WITH (ABS)DIAGONAL OF A (IF PC=1 OR PC=2)     |
!| PC       |-->| PRECONDITIONING OPTION                               |
!| WH       |<->| ORTONORMAL BASE OF KRYLOV (W) AND HESSENBERG MATRIX  |
!| VB       |<->| VECTOR WITH BETA COEFFICIENTS FOR HOUSEHOLDER P's    |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_ARNOLDIHOUSE => CSC_ARNOLDIHOUSE
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE(CSC_OBJ), INTENT(IN) :: MA
      INTEGER, INTENT(IN) :: M,PC
      DOUBLE PRECISION, DIMENSION(MA%NR), INTENT(IN) :: VV
      DOUBLE PRECISION, DIMENSION(MA%NC), INTENT(IN) :: VD
      DOUBLE PRECISION, DIMENSION(MA%NZ), INTENT(IN) :: LUV
      DOUBLE PRECISION, DIMENSION(MA%NR+1,M+1), INTENT(INOUT) :: WH
      DOUBLE PRECISION, DIMENSION(M), INTENT(INOUT) :: VB
!
!  IN SUBROUTINE VARIABLES
!
      DOUBLE PRECISION, DIMENSION(MA%NR) :: VZ, VZAUX, VQ
      DOUBLE PRECISION :: BQTW, BZTW
      INTEGER :: N,I,J,K
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!-----------------------------------------------------------------------
! NO PRECONDITIONING
!-----------------------------------------------------------------------
      IF(PC==0)THEN
! INIT
      N = SIZE(VV)
      IF(M.GT.N)THEN
        WRITE(*,*) "WARNING: SIZE OF KRYLOV SUSPACE GREATER THAN SIZE OF
     &   MATRIX A"
        STOP 0
      ENDIF
      VZ = VV
      WH = 0.0D0
      VB = 0.0D0
      CALL HOUSEHOLDERV(VZ,VB(1),WH(2:N+1,1),N)
!
!-----------------------------------------------------------------------
      IF(M.LT.N)THEN
      ! CASE WHEN M<N
!
      DO J=1:M+1
      ! COMPUTING H
      WH(1:J-1,J) = VZ(1:J-1)
      WH(J,J) = VZ(J) - VB(J)*(DOT_PRODUCT(VZ(J:N),WH(J+1:N+1,J)))
!
      IF(J.LE.M)THEN
!
      !COMPUTING V_J, BASED IN WH PREVIOUS RESULTS
      !V_J = P1*P2*...*PJ*ej ,i.e.,SUCCESIVE RIGHT OUTER PRODUCTS
        VV = 0.0D0
        VQ = 0.0D0
        VQ(J) = 1.0D0
        DO K=J,1,-1
          BQTW = VB(K)*(DOT_PRODUCT(VQ(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VV(I) = VQ(I) - WH(I+1,K)*BQTW
          ENDDO
          VQ = VV
        ENDDO
!
      !COMPUTING Z_J+1
      !Z_J+1 = PJ*PJ-1*PJ-2*...*P2*P1*(A*VV_J)
      !SUCCESIVE OUTER PRODUCTS BY LEFT SIDE
        CALL CSC_MMATVEC(MA,VV,VZAUX,N)
        DO K=1,J
          BZTW = VB(K)*(DOT_PRODUCT(VZAUX(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VZ(I) = VZAUX(I) - WH(I+1,K)*BZTW
          ENDDO
          VZAUZ = VV
        ENDDO
!
      !COMPUTING NEW HOUSEHOLDER PROJECTOR VECTOR AND NEW BETA
        CALL HOUSEHOLDERV(VZ(J+1:N),VB(J+1),WH(J+2:N+1,J+1),N-J)
      ENDIF
!      
      ENDDO
!-----------------------------------------------------------------------
      ELSE IF(M.EQ.N)THEN
      ! CASE WHEN M==N
!
      DO J=1:M+1
      ! COMPUTING H
      WH(1:J-1,J) = VZ(1:J-1)
!
      IF(J.LE.M)THEN
      WH(J,J) = VZ(J) - VB(J)*(DOT_PRODUCT(VZ(J:N),WH(J+1:N+1,J)))
!
      !COMPUTING V_J, BASED IN WH PREVIOUS RESULTS
      !V_J = P1*P2*...*PJ*ej ,i.e.,SUCCESIVE RIGHT OUTER PRODUCTS
        VV = 0.0D0
        VQ = 0.0D0
        VQ(J) = 1.0D0
        DO K=J,1,-1
          BQTW = VB(K)*(DOT_PRODUCT(VQ(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VV(I) = VQ(I) - WH(I+1,K)*BQTW
          ENDDO
          VQ = VV
        ENDDO
!
      !COMPUTING Z_J+1
      !Z_J+1 = PJ*PJ-1*PJ-2*...*P2*P1*(A*VV_J)
      !SUCCESIVE OUTER PRODUCTS BY LEFT SIDE
        CALL CSC_MMATVEC(MA,VV,VZAUX,N)
        DO K=1,J
          BZTW = VB(K)*(DOT_PRODUCT(VZAUX(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VZ(I) = VZAUX(I) - WH(I+1,K)*BZTW
          ENDDO
          VZAUX = VZ
        ENDDO
!
      !COMPUTING NEW HOUSEHOLDER PROJECTOR VECTOR AND NEW BETA
        CALL HOUSEHOLDERV(VZ(J+1:N),VB(J+1),WH(J+2:N+1,J+1),N-J)
      ENDIF
!      
      ENDDO
!
      ENDIF
!-----------------------------------------------------------------------
! DIAGONAL OR ABSOULTE-DIAGONAL PRECONDITIONING
!-----------------------------------------------------------------------
      ELSE IF((PC.EQ.1).OR.(PC.EQ.2))THEN
! INIT
      N = SIZE(VV)
      IF(M.GT.N)THEN
        WRITE(*,*) "WARNING: SIZE OF KRYLOV SUSPACE GREATER THAN SIZE OF
     &   MATRIX A"
        STOP 0
      ENDIF
      VZ = VV
      WH = 0.0D0
      VB = 0.0D0
      CALL HOUSEHOLDERV(VZ,VB(1),WH(2:N+1,1),N)
!
!-----------------------------------------------------------------------
      IF(M.LT.N)THEN
      ! CASE WHEN M<N
!
      DO J=1:M+1
      ! COMPUTING H
      WH(1:J-1,J) = VZ(1:J-1)
      WH(J,J) = VZ(J) - VB(J)*(DOT_PRODUCT(VZ(J:N),WH(J+1:N+1,J)))
!
      IF(J.LE.M)THEN
!
      !COMPUTING V_J, BASED IN WH PREVIOUS RESULTS
      !V_J = P1*P2*...*PJ*ej ,i.e.,SUCCESIVE RIGHT OUTER PRODUCTS
        VV = 0.0D0
        VQ = 0.0D0
        VQ(J) = 1.0D0
        DO K=J,1,-1
          BQTW = VB(K)*(DOT_PRODUCT(VQ(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VV(I) = VQ(I) - WH(I+1,K)*BQTW
          ENDDO
          VQ = VV
        ENDDO
!
      !COMPUTING Z_J+1
      !Z_J+1 = PJ*PJ-1*PJ-2*...*P2*P1*(A*VV_J)
      !SUCCESIVE OUTER PRODUCTS BY LEFT SIDE
        VZAUX = VV/VD    !RIGHT PRECONDITIONING
        CALL CSC_MMATVEC(MA,VZAUX,VZAUX,N)
        DO K=1,J
          BZTW = VB(K)*(DOT_PRODUCT(VZAUX(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VZ(I) = VZAUX(I) - WH(I+1,K)*BZTW
          ENDDO
          VZAUX = VZ
        ENDDO
!
      !COMPUTING NEW HOUSEHOLDER PROJECTOR VECTOR AND NEW BETA
        CALL HOUSEHOLDERV(VZ(J+1:N),VB(J+1),WH(J+2:N+1,J+1),N-J)
      ENDIF
!      
      ENDDO
!-----------------------------------------------------------------------
      ELSE IF(M.EQ.N)THEN
      ! CASE WHEN M==N
!
      DO J=1:M+1
      ! COMPUTING H
      WH(1:J-1,J) = VZ(1:J-1)
!
      IF(J.LE.M)THEN
      WH(J,J) = VZ(J) - VB(J)*(DOT_PRODUCT(VZ(J:N),WH(J+1:N+1,J)))
!
      !COMPUTING V_J, BASED IN WH PREVIOUS RESULTS
      !V_J = P1*P2*...*PJ*ej ,i.e.,SUCCESIVE RIGHT OUTER PRODUCTS
        VV = 0.0D0
        VQ = 0.0D0
        VQ(J) = 1.0D0
        DO K=J,1,-1
          BQTW = VB(K)*(DOT_PRODUCT(VQ(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VV(I) = VQ(I) - WH(I+1,K)*BQTW
          ENDDO
          VQ = VV
        ENDDO
!
      !COMPUTING Z_J+1
      !Z_J+1 = PJ*PJ-1*PJ-2*...*P2*P1*(A*VV_J)
      !SUCCESIVE OUTER PRODUCTS BY LEFT SIDE
        VZAUX = VV/VD    !RIGTH PRECONDITIONING
        CALL CSC_MMATVEC(MA,VZAUX,VZAUX,N)
        DO K=1,J
          BZTW = VB(K)*(DOT_PRODUCT(VZAUX(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VZ(I) = VZAUX(I) - WH(I+1,K)*BZTW
          ENDDO
          VZAUX = VZ
        ENDDO
!
      !COMPUTING NEW HOUSEHOLDER PROJECTOR VECTOR AND NEW BETA
        CALL HOUSEHOLDERV(VZ(J+1:N),VB(J+1),WH(J+2:N+1,J+1),N-J)
      ENDIF
!      
      ENDDO
!
      ENDIF
!-----------------------------------------------------------------------
! SSOR OR ILU(0) PRECONDITIONING
!-----------------------------------------------------------------------
      ELSE IF((PC.EQ.3).OR.(PC.EQ.4))THEN
! INIT
      N = SIZE(VV)
      IF(M.GT.N)THEN
        WRITE(*,*) "WARNING: SIZE OF KRYLOV SUSPACE GREATER THAN SIZE OF
     &   MATRIX A"
        STOP 0
      ENDIF
      VZ = VV
      WH = 0.0D0
      VB = 0.0D0
      CALL HOUSEHOLDERV(VZ,VB(1),WH(2:N+1,1),N)
!
!-----------------------------------------------------------------------
      IF(M.LT.N)THEN
      ! CASE WHEN M<N
!
      DO J=1:M+1
      ! COMPUTING H
      WH(1:J-1,J) = VZ(1:J-1)
      WH(J,J) = VZ(J) - VB(J)*(DOT_PRODUCT(VZ(J:N),WH(J+1:N+1,J)))
!
      IF(J.LE.M)THEN
!
      !COMPUTING V_J, BASED IN WH PREVIOUS RESULTS
      !V_J = P1*P2*...*PJ*ej ,i.e.,SUCCESIVE RIGHT OUTER PRODUCTS
        VV = 0.0D0
        VQ = 0.0D0
        VQ(J) = 1.0D0
        DO K=J,1,-1
          BQTW = VB(K)*(DOT_PRODUCT(VQ(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VV(I) = VQ(I) - WH(I+1,K)*BQTW
          ENDDO
          VQ = VV
        ENDDO
!
      !COMPUTING Z_J+1
      !Z_J+1 = PJ*PJ-1*PJ-2*...*P2*P1*(A*VV_J)
      !SUCCESIVE OUTER PRODUCTS BY LEFT SIDE
        CALL CSC_SOLPACKLU(N,N,LUV,MA%R,MA%C,VV,VZAUX)
        CALL CSC_MMATVEC(MA,VZAUX,VZAUX,N)
        DO K=1,J
          BZTW = VB(K)*(DOT_PRODUCT(VZAUX(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VZ(I) = VZAUX(I) - WH(I+1,K)*BZTW
          ENDDO
          VZAUX = VZ
        ENDDO
!
      !COMPUTING NEW HOUSEHOLDER PROJECTOR VECTOR AND NEW BETA
        CALL HOUSEHOLDERV(VZ(J+1:N),VB(J+1),WH(J+2:N+1,J+1),N-J)
      ENDIF
!      
      ENDDO
!-----------------------------------------------------------------------
      ELSE IF(M.EQ.N)THEN
      ! CASE WHEN M==N
!
      DO J=1:M+1
      ! COMPUTING H
      WH(1:J-1,J) = VZ(1:J-1)
!
      IF(J.LE.M)THEN
      WH(J,J) = VZ(J) - VB(J)*(DOT_PRODUCT(VZ(J:N),WH(J+1:N+1,J)))
!
      !COMPUTING V_J, BASED IN WH PREVIOUS RESULTS
      !V_J = P1*P2*...*PJ*ej ,i.e.,SUCCESIVE RIGHT OUTER PRODUCTS
        VV = 0.0D0
        VQ = 0.0D0
        VQ(J) = 1.0D0
        DO K=J,1,-1
          BQTW = VB(K)*(DOT_PRODUCT(VQ(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VV(I) = VQ(I) - WH(I+1,K)*BQTW
          ENDDO
          VQ = VV
        ENDDO
!
      !COMPUTING Z_J+1
      !Z_J+1 = PJ*PJ-1*PJ-2*...*P2*P1*(A*VV_J)
      !SUCCESIVE OUTER PRODUCTS BY LEFT SIDE
        CALL CSC_SOLPACKLU(N,N,LUV,MA%R,MA%C,VV,VZAUX)
        CALL CSC_MMATVEC(MA,VZAUX,VZAUX,N)
        DO K=1,J
          BZTW = VB(K)*(DOT_PRODUCT(VZAUX(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VZ(I) = VZAUX(I) - WH(I+1,K)*BZTW
          ENDDO
          VZAUX = VZ
        ENDDO
!
      !COMPUTING NEW HOUSEHOLDER PROJECTOR VECTOR AND NEW BETA
        CALL HOUSEHOLDERV(VZ(J+1:N),VB(J+1),WH(J+2:N+1,J+1),N-J)
      ENDIF
!      
      ENDDO
!
      ENDIF
!-----------------------------------------------------------------------
      ENDIF
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_ARNOLDIHOUSE
