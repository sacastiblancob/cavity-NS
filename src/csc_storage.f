!                     ******************
                      MODULE CSC_STORAGE
!                     ******************
!
!
!***********************************************************************
! STRUCTURE DECLARATION FOR CSC STRUCTURES
!***********************************************************************
!
!brief    STRUCTURE CSC FOR STORAGE VALUES, ROWS INDICES AND COLUMN
!         STARTS
!
!history  Sergio Castiblanco
!+        16/01/2021
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_CSC
      INTERFACE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      SUBROUTINE ALL_CSC(OBJ, N, M, NZ, NAM)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(INOUT) :: OBJ
        INTEGER, INTENT(IN)          :: N, M, NZ
        CHARACTER(LEN=6), INTENT(IN) :: NAM
      END SUBROUTINE
!
      SUBROUTINE LOG2INT(DIMI,LOGI,INTE)
        INTEGER, INTENT(IN) :: DIMI
        INTEGER, DIMENSION(DIMI), INTENT(OUT) :: INTE
        LOGICAL, DIMENSION(DIMI), INTENT(IN) :: LOGI
      END SUBROUTINE
!
      SUBROUTINE LOG2IND(DIMI,LOGI,INTE)
        INTEGER, INTENT(IN) :: DIMI
        INTEGER, ALLOCATABLE, INTENT(OUT) :: INTE(:)
        LOGICAL, DIMENSION(DIMI), INTENT(IN) :: LOGI
      END SUBROUTINE
! 
      SUBROUTINE CSC_DIAG(NCH,NRH,H2,D,MA,NAMMA)
        USE DECLARATIONS_CSC
        INTEGER, INTENT(IN)          :: NCH
        INTEGER, INTENT(IN)          :: NRH
        DOUBLE PRECISION, INTENT(IN), DIMENSION(NRH,NCH) :: H2
        INTEGER, INTENT(IN), DIMENSION(NCH)              :: D
        TYPE(CSC_OBJ), INTENT(INOUT)                     :: MA
        CHARACTER(LEN=6), INTENT(IN)                     :: NAMMA
      END SUBROUTINE CSC_DIAG
!
      SUBROUTINE CSC_DIAGA(MA,DA)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA
        DOUBLE PRECISION, DIMENSION(MA%NC), INTENT(INOUT) :: DA
      END SUBROUTINE CSC_DIAGA
!
      SUBROUTINE CSC_KRON(MA,MB,MC,NAMC)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA, MB
        TYPE(CSC_OBJ), INTENT(INOUT) :: MC
        CHARACTER(LEN=6), INTENT(IN) :: NAMC
      END SUBROUTINE CSC_KRON
!
      SUBROUTINE CSC_SUM(MA,MB,MC,NAMC)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA, MB
        TYPE(CSC_OBJ), INTENT(INOUT) :: MC
        CHARACTER(LEN=6), INTENT(IN) :: NAMC
      END SUBROUTINE CSC_SUM
!
      SUBROUTINE UNION(VA,VB,VC)
        INTEGER, INTENT(IN) :: VA(:), VB(:)
        INTEGER, ALLOCATABLE, INTENT(OUT) :: VC(:)
      END SUBROUTINE UNION
!
      SUBROUTINE CSC_MMATVEC(MA,BV,CV,NE)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA
        INTEGER, INTENT(IN) :: NE
        DOUBLE PRECISION, INTENT(IN), DIMENSION(NE) :: BV
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(NE) :: CV
      END SUBROUTINE CSC_MMATVEC
!
      SUBROUTINE CSC_PRECONSSOR(PQV,MA,W)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA
        DOUBLE PRECISION, DIMENSION(MA%NZ), INTENT(INOUT) :: PQV
        DOUBLE PRECISION, INTENT(IN) :: W
      END SUBROUTINE CSC_PRECONSSOR
!
      SUBROUTINE CSC_PRECONILU0(LUV,MA)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA
        DOUBLE PRECISION, DIMENSION(MA%NZ), INTENT(INOUT) :: LUV
      END SUBROUTINE CSC_PRECONILU0
!
      SUBROUTINE CSC_SOLPACKLU(N,M,NZ,LUV,LUR,LUC,VB,VX)
        INTEGER, INTENT(IN) :: N,M,NZ
        DOUBLE PRECISION, INTENT(IN), DIMENSION(NZ) :: LUV
        INTEGER, INTENT(IN), DIMENSION(NZ) :: LUR
        INTEGER, INTENT(IN), DIMENSION(M+1) :: LUC
        DOUBLE PRECISION, INTENT(IN), DIMENSION(N) :: VB
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(N) :: VX
      END SUBROUTINE CSC_SOLPACKLU
!
      SUBROUTINE HOUSEHOLDERV(VX,BETA,V,N)
        INTEGER, INTENT(IN) :: N
        DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: VX
        DOUBLE PRECISION, DIMENSION(N), INTENT(INOUT) :: V
        DOUBLE PRECISION, INTENT(OUT) :: BETA
      END SUBROUTINE HOUSEHOLDERV
!
      SUBROUTINE CSC_ARNOLDIHOUSE(MA,VVV,M,LUV,VD,PC,WH,VB)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA
        INTEGER, INTENT(IN) :: M,PC
        DOUBLE PRECISION, DIMENSION(MA%NR), INTENT(IN) :: VVV
        DOUBLE PRECISION, DIMENSION(MA%NC), INTENT(IN) :: VD
        DOUBLE PRECISION, DIMENSION(MA%NZ), INTENT(IN) :: LUV
        DOUBLE PRECISION, DIMENSION(MA%NR+1,M+1), INTENT(INOUT) :: WH
        DOUBLE PRECISION, DIMENSION(M), INTENT(INOUT) :: VB
      END SUBROUTINE CSC_ARNOLDIHOUSE
!
      SUBROUTINE CSC_CG(MA,BV,XV,MAXNITER,NITER,TOL,NE)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA
        INTEGER, INTENT(IN) :: NE
        DOUBLE PRECISION, INTENT(IN), DIMENSION(NE) :: BV
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NE) :: XV
        INTEGER, INTENT(IN) :: MAXNITER
        INTEGER, INTENT(OUT) :: NITER
        DOUBLE PRECISION, INTENT(IN) :: TOL
      END SUBROUTINE CSC_CG
!
      SUBROUTINE CSC_TRANS(MA,MT,NAMT)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA
        TYPE(CSC_OBJ), INTENT(INOUT) :: MT
        CHARACTER(LEN=6), INTENT(IN) :: NAMT
      END SUBROUTINE CSC_TRANS
!
      SUBROUTINE CSC_PRESOR(MA,MP,MQ,W,NAMP,NAMQ)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA
        TYPE(CSC_OBJ), INTENT(INOUT) :: MP,MQ
        DOUBLE PRECISION, INTENT(IN) :: W
        CHARACTER(LEN=6), INTENT(IN) :: NAMP,NAMQ
      END SUBROUTINE CSC_PRESOR
!
      SUBROUTINE CSC_SOR(MA,MP,MQ,BV,XV,MAXNITER,NITER,TOL,NE)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA,MP,MQ
        INTEGER, INTENT(IN) :: NE
        DOUBLE PRECISION, INTENT(IN), DIMENSION(NE) :: BV
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NE) :: XV
        INTEGER, INTENT(IN) :: MAXNITER
        INTEGER, INTENT(OUT) :: NITER
        DOUBLE PRECISION, INTENT(IN) :: TOL
      END SUBROUTINE CSC_SOR
!
      SUBROUTINE CSC_PCG(MA,VB,VX,NITER,TOL,LUV,VD,PC,NT,RES)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA
        INTEGER, INTENT(IN) :: PC,NITER
        INTEGER, INTENT(OUT) :: NT
        DOUBLE PRECISION, INTENT(IN) :: TOL
        DOUBLE PRECISION, INTENT(OUT) :: RES
        DOUBLE PRECISION, DIMENSION(MA%NR), INTENT(IN) :: VB
        DOUBLE PRECISION, DIMENSION(MA%NC), INTENT(IN) :: VD
        DOUBLE PRECISION, DIMENSION(MA%NZ), INTENT(IN) :: LUV
        DOUBLE PRECISION, DIMENSION(MA%NC), INTENT(INOUT) :: VX
      END SUBROUTINE CSC_PCG
!
      SUBROUTINE CSC_BICGSTAB(MA,VB,VX,NITER,TOL,LUV,VD,PC,NT,RES)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA
        INTEGER, INTENT(IN) :: PC,NITER
        INTEGER, INTENT(OUT) :: NT
        DOUBLE PRECISION, INTENT(IN) :: TOL
        DOUBLE PRECISION, INTENT(OUT) :: RES
        DOUBLE PRECISION, DIMENSION(MA%NR), INTENT(IN) :: VB
        DOUBLE PRECISION, DIMENSION(MA%NC), INTENT(IN) :: VD
        DOUBLE PRECISION, DIMENSION(MA%NZ), INTENT(IN) :: LUV
        DOUBLE PRECISION, DIMENSION(MA%NC), INTENT(INOUT) :: VX
      END SUBROUTINE CSC_BICGSTAB
!
      SUBROUTINE CSC_GMRES(MA,VB,VX,M,NITER,TOL,LUV,VD,PC,NT,RES)
        USE DECLARATIONS_CSC
        TYPE(CSC_OBJ), INTENT(IN) :: MA
        INTEGER, INTENT(IN) :: M,PC,NITER
        INTEGER, INTENT(OUT) :: NT
        DOUBLE PRECISION, INTENT(IN) :: TOL
        DOUBLE PRECISION, INTENT(OUT) :: RES
        DOUBLE PRECISION, DIMENSION(MA%NR), INTENT(IN) :: VB
        DOUBLE PRECISION, DIMENSION(MA%NC), INTENT(IN) :: VD
        DOUBLE PRECISION, DIMENSION(MA%NZ), INTENT(IN) :: LUV
        DOUBLE PRECISION, DIMENSION(MA%NC), INTENT(INOUT) :: VX
      END SUBROUTINE CSC_GMRES
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END INTERFACE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END MODULE CSC_STORAGE

