!                   *******************************
                    SUBROUTINE POINT_REGULARIZATION
!                   *******************************
     & (SING,TOLSING,MNITERM,W,LM,LPM,LQM,UL,VL,RM)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!***********************************************************************
! NAVIER-STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) COMPUTING REGULARIZATION MATRIX AND PREPARING EVERYTING FOR
!            SOLVING THE PRESURE THROUGH REGULARIZATION AND SOR METHOD
!
!history  Sergio Castiblanco
!+        25/01/2021
!+        Translation for original Matlab implementation
!+        16/02/2022
!+        Removing varible ISUSERW, not needed anymore
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| SING      |-->| NUMBER OF ITERATIONS FOR INVERSE POWER METHOD       |
!| TOLSING   |-->| TOLERANCE FOR INVERSE POWER METHOD SOLVER           |
!| MNITERM   |-->| MAXIMUM NUMBER OF ITERATIONS FOR IPM SOLVER         |
!| W         |<->| SOR OVER RELAXATION COEFFICIENT                     |
!| ML        |-->| LAPLACIAN MATRIX                                    |
!| MLP       |<->| P MATRIX FOR SOR SOLVER                             |
!| MLQ       |<->| Q MATRIX FOR SOR SOLVER                             |
!| UL        |<->| LHS SINGULAR VECTOR, EIGENVECTOR OF L'              |
!| UV        |<--| RHS SINGULAR VECTOR, EIGENVECTOR OF L               |
!| RM        |<--| REGULARIZATION MATRIX = (I - UL*UL')                |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE CSC_STORAGE
      USE DECLARATIONS_NUMERICAL, ONLY: NX,NY,NTIM,DEBUG
!
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! INTENT VARIABLES
!
      INTEGER, INTENT(IN) :: SING, MNITERM
      DOUBLE PRECISION, INTENT(IN) :: TOLSING
      DOUBLE PRECISION, INTENT(INOUT) :: W
      TYPE(CSC_OBJ), INTENT(IN) :: LM
      TYPE(CSC_OBJ), INTENT(INOUT) :: LPM, LQM
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(INOUT) :: UL
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(OUT) :: VL
      DOUBLE PRECISION, DIMENSION(NX*NY,NX*NY), INTENT(OUT) :: RM
!
! IN SUBROUTINE VARIABLES
!
      !LOOP VARIABLES
      INTEGER :: I
      !TRANSPOSE OF L
      TYPE(CSC_OBJ) :: LTM
      !RIGHT HAND SIDE FOR INVERSE POWER METHOD
      DOUBLE PRECISION, DIMENSION(NX*NY) :: B
      !UL ROW AND UL COLUMN
      DOUBLE PRECISION, DIMENSION(NX*NY,1) :: ULC
      DOUBLE PRECISION, DIMENSION(1,NX*NY) :: ULR
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
! TRANSPOSING L
!
      CALL CSC_TRANS(LM,LTM,'LTM   ')
!
! COMPUTING OVER-RELAXATION COEFFICIENT
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING W'
      W = 13.4523570058092D0 * EXP(-0.206450260650164D0 * 
     &    (NX*NY)**(-0.434163866503769D0)) - 11.4497834449085D0
      WRITE(*,*) REPEAT('~',72)
      WRITE(*,*) 'SOR OVER-RELAXATION COEFFICIENT :',W
      WRITE(*,*) REPEAT('~',72)
!
! COMPUTING LP AND LQ FOR SOR SOLVER
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING LP AND LQ FOR SOR'
      CALL CSC_PRESOR(LM,LPM,LQM,W,'LPM   ','LQM   ')
!
! GOING TO INVERSE POWER METHOD
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING UL'
      UL = 0.5D0
      B = 0D0
      NTIM = 0
      DO I = 1,SING
        CALL CSC_CG(LTM,B,UL,MNITERM,NTIM,TOLSING,SIZE(UL))
        UL = UL*(1D0/NORM2(UL))
      ENDDO 
      !FILLING ULC AND ULR
      DO I = 1,NX*NY
        ULC(I,1) = UL(I)
        ULR(1,I) = UL(I)
      ENDDO
!
! COMPUTING RHS SINGULAR VECTOR VL
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING VL'
      VL = (1D0/SQRT(REAL(NX*NY,8)))
!
! COMPUTING REGULARIZATION MATRIX
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING REGULARIZATION MATRIX'
      RM = -1D0*MATMUL(ULC,ULR)
      DO I = 1,NX*NY
        RM(I,I) = 1D0 + RM(I,I)
      ENDDO
!
! WRITING THE NUMBER OF ITERATIONS CARRIED OUT BY THE SOLVER
!
      WRITE(*,*) REPEAT('~',72)
      WRITE(*,*) 'NUMBER OF ITERATIONS FOR INVERSE POWER METHOD : ',SING
      WRITE(*,*) 'NUMBER OF MAX. ITERATIONS FOR SOLVER: ',MNITERM
      WRITE(*,*) 'TOLERANCE FOR THE RESIDUAL: ',TOLSING
      WRITE(*,*) 'NUMBER OT ITERATIONS CARRIED OUT BY THE SOLVER: ',NTIM
      WRITE(*,*) 'INFINITE NORM OF LHS SINGULAR VECTOR (UL): ',
     &    MAXVAL(ABS(UL))
      WRITE(*,*) REPEAT('~',72)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE POINT_REGULARIZATION
