!                     ********************
                      SUBROUTINE CSC_GMRES
!                     ********************
     & (MA,VB,VX,M,NITER,TOL,LUV,VD,PC,NT,RES)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) SOLVES SYSTEM Ax=b WITH PRECONDITIONED GENERALIZED MINIMAL
!            RESIDUAL (GMRES) METHOD (RESTARTED) 
!
!history  Sergio Castiblanco
!+        17/02/2022
!+        Translation for original Matlab implementation
!
!+
! This function solves the system Ax=b with the Restarted Generalized
! Minimal Residual (Restarted GMRES)
!
! Entries:
!     MA : Matrix A in CSC storage
!     VB : Right hand side vector
!     VX : Final solution (must have the first guest for the solution)
!     M : Krylov space dimension
!     NITER : Max. number of iterations
!     TOL : tolerance for the stop through the norm of the residual
!     LUV  : values of LU decomposition matrix (for SSOR or ILU(0)) in
!             CSC_packed storage, if not SSOR or ILU(0), it must be have
!             arbitrary values.
!     VD  : DIAGONAL VALUES OF A, FOR (ABS)DIAGONAL PRECONDITIONINGS
!     PC  : preconditioning type
!       - 0 : No preconditioning
!       - 1 : Diagonal preconditioning
!       - 2 : Absolute diagonal preconditioning
!       - 3 : Symmetric SOR (Symmetric Gauss-Seidel -> SSOR, w=1)
!       - 4 : ILU(0)
!
!      Sergio A. Castiblanco B.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| MA       |-->| MATRIX A IN CSC STORAGE                              |
!| VB       |-->| RIGTH HAND SIDE VECTOR                               |
!| VX       |<->| SOLUTION VECTOR, AND INTIAL GUESS                    |
!| M        |-->| KRYLOV SUBSPACE SIZE                                 |
!| NITER    |-->| MAXIMUM NUMBER OF ITERATIONS                         |
!| TOL      |-->| TOLERANCE FOR THE NORM OF RESIDUAL                   |
!| LUV      |-->| LU DECOMPOSITION VALUES (FOR SSOR OR ILU0)           |
!| VD       |-->| VECTOR WITH (ABS)DIAGONAL OF A (IF PC=1 OR PC=2)     |
!| PC       |-->| PRECONDITIONING OPTION                               |
!| NT       |<--| FINAL NUMBER OF OUTER ITERATIONS                     |
!| RES      |<--| FINAL NORM OF THE RESIDUAL                           |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_GMRES => CSC_GMRES
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE(CSC_OBJ), INTENT(IN) :: MA
      INTEGER, INTENT(IN) :: M,PC,NITER
      INTEGER, INTENT(OUT) :: NT
      DOUBLE PRECISION, INTENT(IN) :: TOL
      DOUBLE PRECISION, INTENT(OUT) :: RES
      DOUBLE PRECISION, DIMENSION(MA%NR), INTENT(IN) :: VB
      DOUBLE PRECISION, DIMENSION(MA%NC), INTENT(IN) :: VD
      DOUBLE PRECISION, DIMENSION(MA%NZ), INTENT(IN) :: LUV
      DOUBLE PRECISION, DIMENSION(MA%NC), INTENT(INOUT) :: VX
!
!  IN SUBROUTINE VARIABLES
!
      DOUBLE PRECISION, DIMENSION(MA%NR+1,M+1) :: WH
      DOUBLE PRECISION, DIMENSION(M+1) :: PBET
      DOUBLE PRECISION, DIMENSION(M+1) :: VG
      DOUBLE PRECISION :: D,S,C,HIJ,HI1J
      DOUBLE PRECISION, DIMENSION(MA%NR) :: VZ,VQ,VR
      DOUBLE PRECISION :: BQTW
      INTEGER :: N,NB,T,I,J,K
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! INIT
      WH = 0.0D0
      PBET = 0.0D0
      VG = 0.0D0
      N = MA%NC
      NB = SIZE(VB)
      IF(N.NE.NB)THEN
        WRITE(*,*) "ERROR GMRES, DIMENSSIONS DO NOT AGREE"
        STOP 0
      ENDIF
!-----------------------------------------------------------------------
! NO PRECONDITIONING
!-----------------------------------------------------------------------
      IF(PC==0)THEN
!
      DO T=1,NITER
!
        !RESIDUAL
        CALL CSC_MMATVEC(MA,VX,VR,N)
        VR = VB - VR
!
        !ARNOLDI METHOD WITH HOUSEHOLDER
        CALL CSC_ARNOLDIHOUSE(MA,VR,M,LUV,VD,PC,WH,PBET)
        !IF(T.LE.10) WRITE(*,*) "PBET",PBET
        !IF(T.LE.10) WRITE(*,*) "WH", Wh(1:(M+M),:)
!
        !ELIMINATION PROCESS FOR HESSENBERG MATRIX (GIVEN'S ROTATIONS)
        VG(1) = WH(1,1)
        DO I=1,M
          D = SQRT(WH(I,I+1)**2 + WH(I+1,I+1)**2)
          S = WH(I+1,I+1)/D
          C = WH(I,I+1)/D
          WH(I,I+1) = WH(I,I+1)*C + WH(I+1,I+1)*S
          WH(I+1,I+1) = 0.0D0
          DO J=I+1,M
            !UPDATING HESSENBERG MATRIX
            HIJ = WH(I,J+1)
            HI1J = WH(I+1,J+1)
            WH(I,J+1) = C*HIJ + S*HI1J
            WH(I+1,J+1) = -S*HIJ + C*HI1J
          ENDDO
          !UPDATING RHS VECTOR VG
          VG(I+1) = -S*VG(I)
          VG(I) = C*VG(I)
        ENDDO
!
        !BACKWARD SUBSTITUTION
        DO I=M,1,-1
          VG(I) = VG(I)/WH(I,I+1)
          VG(1:I-1) = VG(1:I-1) - WH(1:I-1,I+1)*VG(I)
        ENDDO
!
        !COMPUTING Z FOR X_M = X0 + Z
        VZ = 0.D0
        DO K = M,1,-1
          VZ(K) = VZ(K) + VG(K)
          VQ = VZ
          BQTW = PBET(K)*(DOT_PRODUCT(VQ(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VZ(I) = VQ(I) - WH(I+1,K)*BQTW
          ENDDO
        ENDDO
!
        !UPDATING THE SOLUTION
        VX = VX + VZ
!
        !TOLERANCE OVER RESIDUAL
        RES = ABS(VG(M+1))
        IF(RES.LT.TOL)THEN
          EXIT
        ENDIF
!
      ENDDO
      NT = T
!
!-----------------------------------------------------------------------
! DIAGONAL OR ABSOULTE-DIAGONAL PRECONDITIONING
!-----------------------------------------------------------------------
      ELSE IF((PC.EQ.1).OR.(PC.EQ.2))THEN
!
      DO T=1,NITER
!
        !RESIDUAL
        CALL CSC_MMATVEC(MA,VX,VR,N)
        VR = VB - VR
!
        !ARNOLDI METHOD WITH HOUSEHOLDER
        CALL CSC_ARNOLDIHOUSE(MA,VR,M,LUV,VD,PC,WH,PBET)
!
        !ELIMINATION PROCESS FOR HESSENBERG MATRIX (GIVEN'S ROTATIONS)
        VG(1) = WH(1,1)
        DO I=1,M
          D = SQRT(WH(I,I+1)**2 + WH(I+1,I+1)**2)
          S = WH(I+1,I+1)/D
          C = WH(I,I+1)/D
          WH(I,I+1) = WH(I,I+1)*C + WH(I+1,I+1)*S
          WH(I+1,I+1) = 0.0D0
          DO J=I+1,M
            !UPDATING HESSENBERG MATRIX
            HIJ = WH(I,J+1)
            HI1J = WH(I+1,J+1)
            WH(I,J+1) = C*HIJ + S*HI1J
            WH(I+1,J+1) = -S*HIJ + C*HI1J
          ENDDO
          !UPDATING RHS VECTOR VG
          VG(I+1) = -S*VG(I)
          VG(I) = C*VG(I)
        ENDDO
!
        !BACKWARD SUBSTITUTION
        DO I=M,1,-1
          VG(I) = VG(I)/WH(I,I+1)
          VG(1:I-1) = VG(1:I-1) - WH(1:I-1,I+1)*VG(I)
        ENDDO
!
        !COMPUTING Z FOR X_M = X0 + Z
        VZ = 0.D0
        DO K = M,1,-1
          VZ(K) = VZ(K) + VG(K)
          VQ = VZ
          BQTW = PBET(K)*(DOT_PRODUCT(VQ(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VZ(I) = VQ(I) - WH(I+1,K)*BQTW
          ENDDO
        ENDDO
!
        !RIGHT PRECONDITIONING
        VZ = VZ/VD
!
        !UPDATING THE SOLUTION
        VX = VX + VZ
!
        !TOLERANCE OVER RESIDUAL
        RES = ABS(VG(M+1))
        IF(RES.LT.TOL)THEN
          EXIT
        ENDIF
!
      ENDDO
      NT = T
!
!-----------------------------------------------------------------------
! SSOR OR ILU(0) PRECONDITIONING
!-----------------------------------------------------------------------
      ELSE IF((PC.EQ.3).OR.(PC.EQ.4))THEN
!
      DO T=1,NITER
!
        !RESIDUAL
        CALL CSC_MMATVEC(MA,VX,VR,N)
        VR = VB - VR
!
        !ARNOLDI METHOD WITH HOUSEHOLDER
        CALL CSC_ARNOLDIHOUSE(MA,VR,M,LUV,VD,PC,WH,PBET)
!
        !ELIMINATION PROCESS FOR HESSENBERG MATRIX (GIVEN'S ROTATIONS)
        VG(1) = WH(1,1)
        DO I=1,M
          D = SQRT(WH(I,I+1)**2 + WH(I+1,I+1)**2)
          S = WH(I+1,I+1)/D
          C = WH(I,I+1)/D
          WH(I,I+1) = WH(I,I+1)*C + WH(I+1,I+1)*S
          WH(I+1,I+1) = 0.0D0
          DO J=I+1,M
            !UPDATING HESSENBERG MATRIX
            HIJ = WH(I,J+1)
            HI1J = WH(I+1,J+1)
            WH(I,J+1) = C*HIJ + S*HI1J
            WH(I+1,J+1) = -S*HIJ + C*HI1J
          ENDDO
          !UPDATING RHS VECTOR VG
          VG(I+1) = -S*VG(I)
          VG(I) = C*VG(I)
        ENDDO
!
        !BACKWARD SUBSTITUTION
        DO I=M,1,-1
          VG(I) = VG(I)/WH(I,I+1)
          VG(1:I-1) = VG(1:I-1) - WH(1:I-1,I+1)*VG(I)
        ENDDO
!
        !COMPUTING Z FOR X_M = X0 + Z
        VZ = 0.D0
        DO K = M,1,-1
          VZ(K) = VZ(K) + VG(K)
          VQ = VZ
          BQTW = PBET(K)*(DOT_PRODUCT(VQ(K:N),WH(K+1:N+1,K)))
          DO I=K,N
            VZ(I) = VQ(I) - WH(I+1,K)*BQTW
          ENDDO
        ENDDO
!
        !RIGHT PRECONDITIONING
        CALL CSC_SOLPACKLU(MA%NR,MA%NC,MA%NZ,LUV,MA%R,MA%C,VZ,VZ)
!
        !UPDATING THE SOLUTION
        VX = VX + VZ
!
        !TOLERANCE OVER RESIDUAL
        RES = ABS(VG(M+1))
        IF(RES.LT.TOL)THEN
          EXIT
        ENDIF
!
      ENDDO
      NT = T
!
!-----------------------------------------------------------------------
      ENDIF
!-----------------------------------------------------------------------
      IF(NT.EQ.NITER)THEN
        WRITE(*,*) "WARNING: GMRES REACH MAX. ITERATIONS ", NT
        WRITE(*,*) "WARNING: RESIDUAL ", RES
      ENDIF
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_GMRES
