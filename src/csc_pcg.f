!                     ******************
                      SUBROUTINE CSC_PCG
!                     ******************
     & (MA,VB,VX,M,NITER,TOL,LUV,VD,PC,NT,RES)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) SOLVES SYSTEM Ax=b WITH PRECONDITIONED CONJUGATE GRADIENT
!             METHOD 
!
!history  Sergio Castiblanco
!+        17/02/2022
!+        Translation for original Matlab implementation
!
!+
! This function solves the system Ax=b with Preconditioned CG
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
      USE CSC_STORAGE, EX_CSC_PCG => CSC_PCG
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
      DOUBLE PRECISION, DIMENSION(MA%NR), INTENT(INOUT) :: VX
!
!  IN SUBROUTINE VARIABLES
!
      DOUBLE PRECISION :: RR,A,B,RZ
      DOUBLE PRECISION, DIMENSION(SIZE(VB)) :: VR,VP,VAP
      DOUBLE PRECISION, ALLOCATABLE :: VZ
      INTEGER :: N,T,NB
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! INIT
      N = MA%NC
      NB = SIZE(VB)
      IF(N.NE.NB)THEN
        WRITE(*,*) "ERROR CG, DIMENSSIONS DO NOT AGREE"
        STOP 0
      ENDIF
!-----------------------------------------------------------------------
! NO PRECONDITIONING
!-----------------------------------------------------------------------
      IF(PC==0)THEN
!
      !RESIDUAL AND P-VECTOR
      CALL CSC_MMATVEC(MA,VX,VR,N)
      VR = VB - VR
      VP = VR
!
      DO T=1,NITER
!
        RR = DOT_PRODUCT(VR,VR)
        CALL CSC_MMATVEC(MA,VP,VAP,N)
        A = RR/DOT_PRODUCT(VAP,VP)
        VX = VX + A*VP
        VR = VR - A*VAP
        RES = NORM2(VR)
        IF(RES.LT.TOL)THEN
          EXIT
        ENDIF
        B = DOT_PRODUCT(VR,VR)/RR
        VP = VR + B*VP
!
      ENDDO
      NT = T
!
!-----------------------------------------------------------------------
! DIAGONAL OR ABSOULTE-DIAGONAL PRECONDITIONING
!-----------------------------------------------------------------------
      ELSE IF((PC.EQ.1).OR.(PC.EQ.2))THEN
!
      !ALLOCATING VZ
      ALLOCATE(VZ,SIZE(VB))
!
      !RESIDUAL AND P-VECTOR
      CALL CSC_MMATVEC(MA,VX,VR,N)
      VR = VB - VR
      VZ = VR/VD
      VP = VZ
!
      DO T=1,NITER
!
        RZ = DOT_PRODUCT(VR,VZ)
        CALL CSC_MMATVEC(MA,VP,VAP,N)
        A = RZ/DOT_PRODUCT(VAP,VP)
        VX = VX + A*VP
        VR = VR - A*VAP
        RES = NORM2(VR)
        IF(RES.LT.TOL)THEN
          EXIT
        ENDIF
        VZ = VR/VD
        B = DOT_PRODUCT(VR,VZ)/RZ
        VP = VZ + B*VP
!
      ENDDO
      NT = T
!
!-----------------------------------------------------------------------
! SSOR OR ILU(0) PRECONDITIONING
!-----------------------------------------------------------------------
      ELSE IF((PC.EQ.3).OR.(PC.EQ.4))THEN
!
      !ALLOCATING VZ
      ALLOCATE(VZ,SIZE(VB))
!
      !RESIDUAL AND P-VECTOR
      CALL CSC_MMATVEC(MA,VX,VR,N)
      VR = VB - VR
      CALL CSC_SOLPACKLU(MA%NR,MA%NC,MA%NZ,LUV,MA%R,MA%C,VR,VZ)
      VP = VZ
!
      DO T=1,NITER
!
        RZ = DOT_PRODUCT(VR,VZ)
        CALL CSC_MATVEC(MA,VP,VAP,N)
        A = RZ/DOT_PRODUCT(VAP,VP)
        VX = VX + A*VP
        VR = VR - A*VAP
        RES = NORM2(VR)
        IF(RES.LT.TOL)THEN
          EXIT
        ENDIF
        CALL CSC_SOLPACKLU(MA%NR,MA%NC,MA%NZ,LUV,MA%R,MA%C,VR,VZ)
        B = DOT_PRODUCT(VR,VZ)/RZ
        VP = VZ + B*VP
!
      ENDDO
      NT = T
!
!-----------------------------------------------------------------------
      ENDIF
!
!-----------------------------------------------------------------------
      IF(NT.EQ.NITER)THEN
        WRITE(*,*) "CG WARNING: MAXIMUM ITERATIONS REACHED ", NT
        WRITE(*,*) "NORM OF RESIDUAL ", RES
      ENDIF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_PCG
