!                   **************************
                    SUBROUTINE PRECONDITIONING
!                   **************************
     & (LM,KM,PCD,WSSORD,PCS,WSSORS,SOLP,PCP,WSSORP,W,LPM,LQM)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!***********************************************************************
! NAVIER-STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) COMPUTING PRECONDITIONERS, IF ANY
!
!history  Sergio Castiblanco
!+        19/02/2022
!+        Translation for original Matlab implementation
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| LM        |-->| POISSON EQUATION MATRIX                             |
!| KM        |-->| DIFUSSION EQUATION MATRIX                           |
!| PCD       |-->| PRECONDITIONING FOR DIFUSSION STEP                  |
!| WSSORD    |-->| SSOR-PRECON. COEFFICIENT FOR DIFUSSION STEP         |
!| PCS       |-->| PRECONDITIONING FOR REGULARIZATION                  |
!| WSSORS    |-->| SSOR-PRECON. COEFFICIENT FOR REGULARIZATION         |
!| SOLP      |-->| SOLVER FOR POISSON STEP                             |
!| PCP       |-->| PRECONDITIONING FOR POISSON SOLVER                  |
!| WSSORP    |-->| SSOR-PRECON. COEFFICIENT FOR REGULARIZATION         |
!| W         |-->| SOR SOLVER COEFFICIENT FOR POISSON                  |
!| LPM       |<->| LOWER PART OF MATRIX LM IF SOR SOLVER               |
!| LQM       |<->| UPPER PART OF MATRIX LM IF SOR SOLVER               |
!| LUVP(mod) |<->| VALUES VECTOR IF SSOR OR ILU0 PRECON., POISSON      |
!| DLM(mod)  |<->| DIAGONAL ENTRIES OF LM IF (ABS)DIAG PRECON.         |
!| LUVR(mod) |<->| VALUES VECTOR IF SSOR OR ILU0 PRECON., REGULARIZATI.|
!| DLMR(mod) |<->| DIAGONAL ENTRIES OF LM IF (ABS)DIAG PRECON. REG.    |
!| LUVD(mod) |<->| VALUES VECTOR IF SSOR OR ILU0 PRECON., DIFUSSION    |
!| DKM(mod)  |<->| DIAGONAL ENTRIES OF KM IF (ABS)DIAG PRECON.         |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE CSC_STORAGE
      USE DECLARATIONS_NUMERICAL, ONLY: NX,NY,DEBUG,LUVP,DLM,LUVR,DLMR,
     &      LUVD,DKM
!
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! INTENT VARIABLES
!
      INTEGER, INTENT(IN) :: PCD,PCS,SOLP,PCP
      DOUBLE PRECISION, INTENT(IN) :: WSSORD,WSSORS,WSSORP
      !DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: LUVP(:),DLM(:)
      !DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: LUVD(:),DKM(:)
      !DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: LUVR(:),DLMR(:)
      DOUBLE PRECISION, INTENT(INOUT) :: W
      TYPE(CSC_OBJ), INTENT(IN) :: LM, KM
      TYPE(CSC_OBJ), INTENT(INOUT) :: LPM, LQM
!
! IN SUBROUTINE VARIABLES
!
      DOUBLE PRECISION, ALLOCATABLE :: ADLM(:),ADKM(:),ADLMR(:)
!
! CHEKING WSSOR
!
      IF(.NOT.((0.0D0.LT.WSSORD).AND.(WSSORD.LT.2.0D0)))THEN
        WRITE(*,*) "INVALID OPTION FOR WSSORD, IT E (0,2): ABORTING!!!"
        STOP 0
      ENDIF
      IF(.NOT.((0.0D0.LT.WSSORS).AND.(WSSORS.LT.2.0D0)))THEN
        WRITE(*,*) "INVALID OPTION FOR WSSORD, IT E (0,2): ABORTING!!!"
        STOP 0
      ENDIF
      IF(.NOT.((0.0D0.LT.WSSORP).AND.(WSSORP.LT.2.0D0)))THEN
        WRITE(*,*) "INVALID OPTION FOR WSSORD, IT E (0,2): ABORTING!!!"
        STOP 0
      ENDIF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! FOR DIFUSSION STEP
!-----------------------------------------------------------------------
!
      IF(DEBUG) WRITE(*,*) "DIFUSSION STEP PRECONDITIONING"
      IF(PCD.EQ.0)THEN         !NO PRECONDITIONING
        CONTINUE
      ELSE IF(PCD.EQ.1)THEN    !DIAGONAL PRECONDITIONING
        ALLOCATE(DKM(KM%NC))
        CALL CSC_DIAGA(KM,DKM)
      ELSE IF(PCD.EQ.2)THEN    !ABSOLUTE DIAGONAL PRECONDITIONING
        ALLOCATE(DKM(KM%NC))
        ALLOCATE(ADKM(KM%NC))
        CALL CSC_DIAGA(KM,ADKM)
        DKM = ABS(ADKM)
      ELSE IF(PCD.EQ.3)THEN    !SSOR PRECONDITIONING
        ALLOCATE(LUVD(KM%NZ))
        CALL CSC_PRECONSSOR(LUVD,KM,WSSORD)
      ELSE IF(PCD.EQ.4)THEN    !ILU(0) PRECONDITIONING (DIRECT LU)
        ALLOCATE(LUVD(KM%NZ))
        CALL CSC_PRECONILU0(LUVD,KM)
      ELSE
        WRITE(*,*) "DIFFUSION INVALID PRECONDITIONING OPTION: ABORTING!"
        STOP 0
      ENDIF
!
!-----------------------------------------------------------------------
! REGULARIZATION (INVERSE POWER METHOD)
!-----------------------------------------------------------------------
!
      IF(DEBUG) WRITE(*,*) "REGULARIZATION STEP PRECONDITIONING"
      IF(PCS.EQ.0)THEN           !NO PRECONDITIONING
        CONTINUE
      ELSE IF(PCS.EQ.1)THEN      !DIAGONAL PRECONDITIONING
        WRITE(*,*) "LM%NC", LM%NC
        ALLOCATE(DLMR(LM%NC))
        CALL CSC_DIAGA(LM,DLMR)
      ELSE IF(PCS.EQ.2)THEN      !ABSOLUTE DIAGONAL PRECONDITIONING
        ALLOCATE(DLMR(LM%NC))
        ALLOCATE(ADLMR(LM%NC))
        CALL CSC_DIAGA(LM,ADLMR)
        DLMR = ABS(ADLMR)
      ELSE IF(PCS.EQ.3)THEN      !SSOR PRECONDITIONING
        ALLOCATE(LUVR(LM%NZ))
        CALL CSC_PRECONSSOR(LUVR,LM,WSSORS)
      ELSE IF(PCS.EQ.4)THEN      !ILU0 PRECONDITIONING
        ALLOCATE(LUVR(LM%NZ))
        CALL CSC_PRECONILU0(LUVR,LM)
      ELSE
        WRITE(*,*) "REGULAR. INVALID PRECONDITIONING OPTION: ABORTING!"
        STOP 0
      ENDIF
!
!-----------------------------------------------------------------------
! POISSON EQUATION
!-----------------------------------------------------------------------
!
      IF(DEBUG) WRITE(*,*) "POISSON STEP PRECONDITIONING"
      IF(PCP.EQ.0)THEN           !NO PRECONDITIONING
        CONTINUE
      ELSE IF(PCP.EQ.1)THEN      !DIAGONAL PRECONDITIONING
        ALLOCATE(DLM(LM%NC))
        CALL CSC_DIAGA(LM,DLM)
      ELSE IF(PCP.EQ.2)THEN      !ABSOLUTE DIAGONAL PRECONDITIONING
        ALLOCATE(DLM(LM%NC))
        ALLOCATE(ADLM(LM%NC))
        CALL CSC_DIAGA(LM,ADLM)
        DLM = ABS(ADLM)
      ELSE IF(PCP.EQ.3)THEN      !SSOR PRECONDITIONING
        ALLOCATE(LUVP(LM%NZ))
        CALL CSC_PRECONSSOR(LUVP,LM,WSSORP)
      ELSE IF(PCP.EQ.4)THEN      !ILU(0) PRECONDITIONING
        ALLOCATE(LUVP(LM%NZ))
        CALL CSC_PRECONILU0(LUVP,LM)
      ELSE
        WRITE(*,*) "POISSON INVALID PRECONDITIONING OPTION: ABORTING!"
        STOP 0
      ENDIF
!
! IF SOR SOLVER
      IF(SOLP.EQ.1)THEN
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
      ENDIF
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE PRECONDITIONING
