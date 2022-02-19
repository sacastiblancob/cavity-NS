!                   ***********************
                    SUBROUTINE CSC_BICGSTAB
!                   ***********************
     & (MA,VB,VX,NITER,TOL,LUV,VD,PC,NT,RES)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) SOLVES SYSTEM Ax=b WITH PRECONDITIONED BI-CONJUGATE
!            GRADIENT METHOD STABILIZED (BICGSTAB)
!
!history  Sergio Castiblanco
!+        18/02/2022
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
      USE CSC_STORAGE, EX_CSC_BICGSTAB => CSC_BICGSTAB
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE(CSC_OBJ), INTENT(IN) :: MA
      INTEGER, INTENT(IN) :: PC,NITER
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
      DOUBLE PRECISION :: A,B,W
      DOUBLE PRECISION, DIMENSION(MA%NR) :: VR,VP,VAP,VR0,VS,VAS
      DOUBLE PRECISION, ALLOCATABLE :: VPG(:),VSG(:)
      INTEGER :: N,T,NB,SVB
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! INIT
      N = MA%NC
      NB = SIZE(VB)
      IF(N.NE.NB)THEN
        WRITE(*,*) "ERROR BICGSTAB, DIMENSSIONS DO NOT AGREE"
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
!     HERE R0 IS TAKEN TO BE A*R0, IN ORDER TO AVOID DOT_PRODUCT=0
      CALL CSC_MMATVEC(MA,VR,VR0,N)
      DO T=1,NITER
!
        CALL CSC_MMATVEC(MA,VP,VAP,N)
        A = DOT_PRODUCT(VR,VR0)/DOT_PRODUCT(VAP,VR0)
        VS = VR - A*VAP
        CALL CSC_MMATVEC(MA,VS,VAS,N)
        W = DOT_PRODUCT(VAS,VS)/DOT_PRODUCT(VAS,VAS)
        VX = VX + A*VP + W*VS
        B = (A/W)*DOT_PRODUCT(VS - W*VAS,VR0)/DOT_PRODUCT(VR,VR0)
        VR = VS - W*VAS
        RES = NORM2(VR)
        IF(RES.LT.TOL)THEN
          EXIT
        ENDIF
        VP = VR + B*(VP - W*VAP)
!
      ENDDO
      NT = T
!
!-----------------------------------------------------------------------
! DIAGONAL OR ABSOULTE-DIAGONAL PRECONDITIONING
!-----------------------------------------------------------------------
      ELSE IF((PC.EQ.1).OR.(PC.EQ.2))THEN
!
      !ALLOCATING VPG, VSG
      SVB = SIZE(VB)
      ALLOCATE(VPG(SVB))
      ALLOCATE(VSG(SVB))
!
      !RESIDUAL AND P-VECTOR
      CALL CSC_MMATVEC(MA,VX,VR,N)
      VR = VB - VR
      VP = VR
!
!     HERE R0 IS TAKEN TO BE A*R0, IN ORDER TO AVOID DOT_PRODUCT=0
      CALL CSC_MMATVEC(MA,VR,VR0,N)
      DO T=1,NITER
!
        VPG = VP/VD
        CALL CSC_MMATVEC(MA,VPG,VAP,N)
        A = DOT_PRODUCT(VR,VR0)/DOT_PRODUCT(VAP,VR0)
        VS = VR - A*VAP
        VSG = VS/VD
        CALL CSC_MMATVEC(MA,VSG,VAS,N)
        W = DOT_PRODUCT(VAS,VS)/DOT_PRODUCT(VAS,VAS)
        VX = VX + A*VPG + W*VSG
        B = (A/W)*DOT_PRODUCT(VS - W*VAS,VR0)/DOT_PRODUCT(VR,VR0)
        VR = VS - W*VAS
        RES = NORM2(VR)
        IF(RES.LT.TOL)THEN
          EXIT
        ENDIF
        VP = VR + B*(VP - W*VAP)
!
      ENDDO
      NT = T
!
!-----------------------------------------------------------------------
! SSOR OR ILU(0) PRECONDITIONING
!-----------------------------------------------------------------------
      ELSE IF((PC.EQ.3).OR.(PC.EQ.4))THEN
!
      !ALLOCATING VPG, VSG
      SVB = SIZE(VB)
      ALLOCATE(VPG(SVB))
      ALLOCATE(VSG(SVB))
!
      !RESIDUAL AND P-VECTOR
      CALL CSC_MMATVEC(MA,VX,VR,N)
      VR = VB - VR
      VP = VR
!
!     HERE R0 IS TAKEN TO BE A*R0, IN ORDER TO AVOID DOT_PRODUCT=0
      CALL CSC_MMATVEC(MA,VR,VR0,N)
      DO T=1,NITER
!
        CALL CSC_SOLPACKLU(MA%NR,MA%NC,MA%NZ,LUV,MA%R,MA%C,VP,VPG)
        CALL CSC_MMATVEC(MA,VPG,VAP,N)
        A = DOT_PRODUCT(VR,VR0)/DOT_PRODUCT(VAP,VR0)
        VS = VR - A*VAP
        CALL CSC_SOLPACKLU(MA%NR,MA%NC,MA%NZ,LUV,MA%R,MA%C,VS,VSG)
        CALL CSC_MMATVEC(MA,VSG,VAS,N)
        W = DOT_PRODUCT(VAS,VS)/DOT_PRODUCT(VAS,VAS)
        VX = VX + A*VPG + W*VSG
        B = (A/W)*DOT_PRODUCT(VS - W*VAS,VR0)/DOT_PRODUCT(VR,VR0)
        VR = VS - W*VAS
        RES = NORM2(VR)
        IF(RES.LT.TOL)THEN
          EXIT
        ENDIF
        VP = VR + B*(VP - W*VAP)
!
      ENDDO
      NT = T
!
!-----------------------------------------------------------------------
      ENDIF
!
!-----------------------------------------------------------------------
      IF(NT.EQ.NITER)THEN
        WRITE(*,*) "BICGSTAB WARNING: MAXIMUM ITERATIONS REACHED ", NT
        WRITE(*,*) "NORM OF RESIDUAL ", RES
      ENDIF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_BICGSTAB
