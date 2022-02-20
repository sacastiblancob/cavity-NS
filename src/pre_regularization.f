!                   *****************************
                    SUBROUTINE PRE_REGULARIZATION
!                   *****************************
     & (LM,SOLSING,PCS,KSSS,SING,TOLSING,MNITERM,LUVR,DLMR,UL,VL)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!***********************************************************************
! NAVIER-STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) COMPUTING RIGHT AND LEFT HAND SIDE SINGULAR VECTORS OF
!            MATRIX LM TO ACHIEVE REGULARIZATION 
!
!history  Sergio Castiblanco
!+        19/02/2022
!+        Translation for original Matlab implementation
!+        16/02/2022
!+        Previous version was awfully made POINT_REGULARIZATION
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| LM        |-->| LAPLACIAN MATRIX                                    |
!| SOLSING   |-->| SOLVER FOR SINGULAR VECTORS PROBLEM                 |
!| PCS       |-->| PRECONDITIONING FOR SINGULAR VECTOR PROBLEM         |
!| KSSS      |-->| KRYLOV SUBSPACE DIMENSION IF SOLVER==GMRES          |
!| SING      |-->| NUMBER OF ITERATIONS FOR INVERSE POWER METHOD       |
!| TOLSING   |-->| TOLERANCE FOR INVERSE POWER METHOD SOLVER           |
!| MNITERM   |-->| MAXIMUM NUMBER OF ITERATIONS FOR IPM SOLVER         |
!| LUVR      |-->| LU VALUES IF SSOR OR ILU(0) PRECON.                 |
!| DLMR      |-->| DIAGONAL ENTRIES OF LM IF (ABS)DIAG PRECON.         |
!| UL        |<->| LHS SINGULAR VECTOR, EIGENVECTOR OF L'              |
!| UV        |<--| RHS SINGULAR VECTOR, EIGENVECTOR OF L               |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE CSC_STORAGE
      USE DECLARATIONS_NUMERICAL, ONLY: NX,NY,DEBUG
!
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! INTENT VARIABLES
!
      INTEGER, INTENT(IN) :: SOLSING,PCS,KSSS,SING,MNITERM
      DOUBLE PRECISION, INTENT(IN) :: TOLSING
      TYPE(CSC_OBJ), INTENT(IN) :: LM
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(INOUT) :: UL
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(OUT) :: VL
      DOUBLE PRECISION, DIMENSION(LM%NC), INTENT(IN) :: DLMR
      DOUBLE PRECISION, DIMENSION(LM%NZ), INTENT(IN) :: LUVR
!
! IN SUBROUTINE VARIABLES
!
      !LOOP VARIABLES
      INTEGER :: I
      !TRANSPOSE OF L
      TYPE(CSC_OBJ) :: LTM
      !RIGHT HAND SIDE FOR INVERSE POWER METHOD
      DOUBLE PRECISION, DIMENSION(NX*NY) :: B
      !FINAL NUMBER OF ITERATIONS TAKEN BY THE SOLVER, AND RESIDUAL
      INTEGER :: NIT
      DOUBLE PRECISION :: RESS
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
! TRANSPOSING L
!
      CALL CSC_TRANS(LM,LTM,'LTM   ')
!
! COMPUTING RHS SINGULAR VECTOR VL
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING VL'
      VL = (1D0/SQRT(REAL(NX*NY,8)))
!
! GOING TO INVERSE POWER METHOD (COMPUTING LHS SINGULAR VECTOR UL)
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING UL'
      UL = VL
      !UL = 0.5D0
      B = 0.0D0
      RESS = 0.0D0
      NIT = 0
!
      IF(SOLSING.EQ.1)THEN          !CONJUGATE GRADIENT
!
      DO I = 1,SING
        CALL CSC_PCG(LTM,B,UL,MNITERM,TOLSING,LUVR,DLMR,PCS,NIT,RESS)
        UL = UL*(1D0/NORM2(UL))
      ENDDO 
!
      ELSE IF(SOLSING.EQ.2)THEN     !BICGSTAB
!
      DO I = 1,SING
        CALL CSC_BICGSTAB(LTM,B,UL,MNITERM,TOLSING,LUVR,DLMR,PCS,NIT,
     &                    RESS)
        UL = UL*(1D0/NORM2(UL))
      ENDDO 
!
      ELSE IF(SOLSING.EQ.3)THEN     !GMRES
!
      DO I = 1,SING
        CALL CSC_GMRES(LTM,B,UL,KSSS,MNITERM,TOLSING,LUVR,DLMR,PCS,
     &                 NIT,RESS)
        UL = UL*(1D0/NORM2(UL))
      ENDDO 
!
      ELSE
        WRITE(*,*) "NO VALID OPTION FOR SOLSING, ABORTING!!!"
        STOP 0
      ENDIF
      
!
! WRITING THE NUMBER OF ITERATIONS CARRIED OUT BY THE SOLVER
!
      WRITE(*,*) REPEAT('~',72)
      WRITE(*,*) 'NUMBER OF ITERATIONS FOR INVERSE POWER METHOD : ',SING
      WRITE(*,*) 'NUMBER OF MAX. ITERATIONS FOR SOLVER: ',MNITERM
      WRITE(*,*) 'TOLERANCE FOR THE RESIDUAL: ',TOLSING
      WRITE(*,*) 'NUMBER OT ITERATIONS CARRIED OUT BY THE SOLVER: ',NIT
      WRITE(*,*) 'INFINITE NORM OF LHS SINGULAR VECTOR (UL): ',
     &    MAXVAL(ABS(UL))
      WRITE(*,*) 'FINAL RESIDUAL OF SOLVER ',PCS,': ',RESS
      WRITE(*,*) REPEAT('~',72)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE PRE_REGULARIZATION
