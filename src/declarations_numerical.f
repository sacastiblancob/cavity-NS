!               *****************************
                MODULE DECLARATIONS_NUMERICAL
!               *****************************
!
!
!***********************************************************************
! NAVIER STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    DECLARATION OF NUMERICAL VARIABLES FOR NS SOLVER
!
!history  Sergio Castiblanco
!+        12/01/2021
!+        Translation for original Matlab implementation
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE CSC_STORAGE
      IMPLICIT NONE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     LOGICALS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! DEBUGGER OPTION
!
      LOGICAL :: DEBUG
!
! IS THE W OVER-RELAXATION GIVEN BY THE USER?
!
      LOGICAL :: ISUSERW
!
! IS THE BOUNDARY CONDITION GIVEN IN A FILE
!
      LOGICAL :: ISBOUND
!
! IS THE HOT START GIVEN IN A FILE
!
      LOGICAL :: ISSTART
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     INTEGERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! LOGICAL UNIT FOR WRITING OUPUTS
!      INTEGER :: LU=6
!
! 
! NUMBER OF NODES IN DIRECTIONS X AND Y
!
      INTEGER :: NX
      INTEGER :: NY
!
! GLOBAL INDEX
!
      INTEGER, ALLOCATABLE :: POS(:)
      INTEGER, ALLOCATABLE :: MPOS(:,:)
!
! BOUNDARY INDEXES
!
      INTEGER, ALLOCATABLE :: BOUND(:), UPBOUND(:), DOBOUND(:),
     &    RIBOUND(:), LEBOUND(:), UPBOUNDI(:), DOBOUNDI(:),
     &    RIBOUNDI(:), LEBOUNDI(:), BOUNDINT(:)
!
! NUMBER OF TIME STEPS
!
      INTEGER :: NT
!
! MAXIMUM NUMBER OF ITERATIONS FOR CG SOLVER (DIFFUSION STEP)
!
      INTEGER :: MNITERD
!
! MAXIMUM NUMBER OF ITERATIONS FOR THE SOLVER USED TO FIND MINIMUM 
! RELATED SINGULAR VECTOR
! (REGULARIZATION STEP)
!
      INTEGER :: MNITERM
!
! NUMBER OF ITERATIONS FOR INVERSE POWER METHOD
! (TO FIND MINIMUM RELATED SINGULAR VECTOR)
!
      INTEGER :: SING
!
! MAXIMUM NUMBER OF ITERATIONS FOR SOR SOLVER (POISSON EQUATION)
!
      INTEGER :: MNITERS
!
! NUMBER OF ITERATIONS TAKEN TO COMPUTE DIFFUSION STEP
!
      INTEGER :: NTIDX,NTIDY
!
! NUMBER OF ITERATIONS TAKEN TO COMPUTE MINIMUM RELATED SINGULAR VECTOR
!
      INTEGER :: NTIM
!
! NUMBER OF ITERATIONS TAKEN TO COMPUTE PRESURE THROUGH POISSON EQUATION
!
      INTEGER :: NITS
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     DOUBLE PRECISION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! DIMENSIONS OF THE DOMAIN
!
      DOUBLE PRECISION :: XMIN
      DOUBLE PRECISION :: XMAX
      DOUBLE PRECISION :: YMIN
      DOUBLE PRECISION :: YMAX
!
! DIFFERENTIALS
!
      DOUBLE PRECISION :: DX
      DOUBLE PRECISION :: DY
!
! TIME DISCRETIZATION
!
!     INITIAL TIME
      DOUBLE PRECISION :: TO
!     FINAL TIME
      DOUBLE PRECISION :: TF
!
! DESIRED COURANT NUMBER
!
      DOUBLE PRECISION :: CFL
!
! GRID
!
      DOUBLE PRECISION, ALLOCATABLE :: X(:), Y(:)
!
! VELOCITIES VECTORS
!
      ! INITIAL CONDITION AND PREVIOUS TIME STEP
      DOUBLE PRECISION, ALLOCATABLE :: UO(:),VO(:)
      ! VELOCITIES AFTER FIRST FRACTIONAL STEP
      DOUBLE PRECISION, ALLOCATABLE :: UP(:),VP(:)
      ! PRESURE VECTORS
      DOUBLE PRECISION, ALLOCATABLE :: P(:), PP(:)
      ! VELOCITIES AFTER SECOND FRACTIONAL STEP
      DOUBLE PRECISION, ALLOCATABLE :: UPP(:), VPP(:)
      ! FINAL RESULT FOR ACTUAL TIME STEP
      DOUBLE PRECISION, ALLOCATABLE :: U(:), V(:)
!
! DERIVATIVES OF VELOCITIES AND PRESURE
!
      !DERIVATIVES FOR NON-LINEAL ADVECTION STEP
      DOUBLE PRECISION, ALLOCATABLE :: DUDX(:),DUDY(:),DVDX(:),DVDY(:)
      !DERIVATIVES FOR RHS OF POISSON EQUATION
      DOUBLE PRECISION, ALLOCATABLE :: DUPDX(:),DVPDY(:)
      !DERIVATIVES OF PRESURE
      DOUBLE PRECISION, ALLOCATABLE :: DPDX(:),DPDY(:)
! TIME AND TIME STEP
!
      DOUBLE PRECISION, ALLOCATABLE  :: T(:)
      DOUBLE PRECISION :: DT
!
! DIFFUSION STABILITY PARAMETERS
!
      DOUBLE PRECISION :: SX
      DOUBLE PRECISION :: SY
!
! SINGULAR VECTORS OF L
!
      DOUBLE PRECISION, ALLOCATABLE :: UL(:),VL(:)
!
! REGULARIZATION MATRIX (I - UL*UL')
!
      DOUBLE PRECISION, ALLOCATABLE :: RM(:,:)
!
! ANGLES OF REGULARIZATION (NOT CORRECTED, NC, AND CORRECTED, C)
!
      DOUBLE PRECISION :: ANGNC, ANGC
!
! SOLVERS TOLERANCES FOR THE RESIDUAL
!
      !CONJUGATE GRADIENT METHOD (DIFFUSION STEP)
      DOUBLE PRECISION :: TOLCG
      !FOR FINDING MINIMUM SINGULAR VECTOR RELATED
      DOUBLE PRECISION :: TOLSING
      !SUCCESIVE OVER-RELAXATION METHOD (POISSON EQUATION)
      DOUBLE PRECISION :: TOLSOR
!
! OVER-RELAXATION COEFFICIENT FOR SOR SOLVER
!
      DOUBLE PRECISION :: W
!
! VALUE TO TIE BOTTOM-LEFT PRESURE
!
      DOUBLE PRECISION :: TIEP
!
! TIME INTERVAL FOR WRITING
!
      DOUBLE PRECISION :: WTIME
      DOUBLE PRECISION :: NTIME = 0D0
!
! BOUNDARY CONDITIONS MATRIX
!
      DOUBLE PRECISION, ALLOCATABLE :: BCOND(:,:)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! CSC OBJECTS (MATRICES)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! MATRIX TO SOLVE POISSON EQUATION
!
      TYPE(CSC_OBJ) :: LM, LPM, LQM
!
! MATRIX TO SOLVE DIFFUSION STEP
!
      TYPE(CSC_OBJ) :: KM
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! CHARACTERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! BOUNDARY CONDITIONS FILE PATH
!
      CHARACTER(LEN=250) :: BOUNDFILE
!
! HOT START FILE PATH
!
      CHARACTER(LEN=250) :: STARTFILE
!
!     ============================================
!
      END MODULE DECLARATIONS_NUMERICAL

