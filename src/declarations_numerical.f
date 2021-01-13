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
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     LOGICALS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! DEBUGGER OPTION
!
      LOGICAL :: DEBUG
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     INTEGERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! LOGICAL UNIT FOR WRITING OUPUTS
      INTEGER :: LU=6
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
!
! BOUNDARY INDEXES
!
      INTEGER, ALLOCATABLE :: BOUND(:), UPBOUND(:), DOBOUND(:),
     &    RIBOUND(:), LEBOUND(:), BOUNDINT(:)
!
! NUMBER OF TIME STEPS
!
      INTEGER :: NT
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
      DOUBLE PRECISION :: CFL = 0.25
!
! GRID
!
      DOUBLE PRECISION, ALLOCATABLE :: X(:), Y(:)
!
! STIFFNESS AND LAPLACE MATRICES
!
      DOUBLE PRECISION, ALLOCATABLE :: K(:,:), L(:,:)
!
! VELOCITIES VECTORS
!
      DOUBLE PRECISION, ALLOCATABLE :: UO(:),VO(:)
!
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
!     ============================================
!
      END MODULE DECLARATIONS_NUMERICAL

