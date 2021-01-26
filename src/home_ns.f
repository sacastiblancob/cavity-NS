!                    ***************
                     PROGRAM HOME_NS
!                    ***************
!
!
!***********************************************************************
! NAVIER STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) MAIN PROGRAM NAVIER STOKES SOLVER.
!
!history  Sergio Castiblanco
!+        12/01/2021
!+        Translation for original Matlab implementation
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_PHYSIC
      USE DECLARATIONS_NUMERICAL
!
      IMPLICIT NONE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER TI
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  READING USER INPUT FILE
!
      NAMELIST /NSCONF/ NU,RHO,NX,NY,XMIN,XMAX,YMIN,YMAX,CFL,TO,TF,
     &    TOLCG,MNITERD,TOLSING,SING,MNITERM,W,ISUSERW,TOLSOR,MNITERS,
     &    DEBUG
!
      OPEN(1010, FILE = "nsconf.nml", STATUS = 'OLD')
      READ(1010, NML = NSCONF)
      CLOSE(1010)
      IF(DEBUG) WRITE(*,*) 'EXIT READING USER ENTRIES'
!
!  WRITING HEADERS
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO WRITE_HEADERS'
      CALL WRITE_HEADERS
      IF(DEBUG) WRITE(*,*) 'EXIT WRITE_HEADERS'
!
!  ALLOCATING MEMORY AND SETTING INITIAL CONDITION
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO POINT_NS'
      CALL POINT_NS
      IF(DEBUG) WRITE(*,*) 'EXIT POINT_NS'
!
!  COMPUTING STIFFNESS DIFUSSION MATRIX
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO DIFFUSION_MATRIX'
      CALL DIFFUSION_MATRIX(DX,DY,DT,NU,NU,KM,SX,SY) 
      IF(DEBUG) WRITE(*,*) 'EXIT DIFFUSION_MATRIX'
!
!  COMPUTING LAPLACIAN OPERATOR MATRIX FOR SOLVE POISSON EQUATION
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO LAPLACIAN_MATRIX'
      CALL LAPLACIAN_MATRIX(DX,DY,LM)
      IF(DEBUG) WRITE(*,*) 'EXIT LAPLACIAN_MATRIX'
!
!  COMPUTING REGULARIZATION MATRIX, AND PREPARING ENTRIES FOR SOR SOLVER
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO POINT_REGULARIZATION'
      CALL POINT_REGULARIZATION(SING,TOLSING,MNITERM,W,ISUSERW,LM,LPM,
     &    LQM,UL,VL,RM)
      IF(DEBUG) WRITE(*,*) 'EXIT POINT_REGULARIZATION'
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  TIME LOOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO TI = 2,SIZE(T)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF(DEBUG) WRITE(*,*) 'TIME LOOP INIT ON TIME = ',T(TI),':',TI
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  FIRST FRACTIONAL STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF(DEBUG) WRITE(*,*) 'GOING INTO FIRST FRACTIONAL STEP'
      CALL DIVER2D(DUDX,DUDY,DVDX,DVDY,UO,VO,BOUND,UPBOUNDI,
     &     DOBOUNDI,LEBOUNDI,RIBOUNDI)
      IF(DEBUG) WRITE(*,*) 'EXIT FIRST FRACTIONAL STEP'
      
!
!  COMPUTING RIGHT HAND SIDE VECTOR
!
!     IF(DEBUG) WRITE(*,*) 'GOING INTO POINT_RHS'
!     CALL POINT_RHS(CO,CA,CS,CCS,CB,BOUND,UPBOUNDI,DOBOUNDI,LEBOUNDI,
!    &    RIBOUNDI,VX,VY)
!     IF(DEBUG) WRITE(*,*) 'EXIT FROM POINT_RHS'
!
!  SOLVING THE SYSTEM OF EQUATIONS WITH CONJUGATE GRADIENT
!
!     IF(DEBUG) WRITE(*,*) 'CALLING SOLVER'
!     CALL CSC_CG(KM,CS,CCS,MAXNITER,NITER,TOL,SIZE(CS))
!     IF(DEBUG) WRITE(*,*) 'EXIT SOLVER'
!
!  UPDATING VARIABLES AND WRITING OUTPUTS
!
!     IF(DEBUG) WRITE(*,*) 'CALLING UPDATE AND WRITING'
!     CALL UPDATE_AND_WRITE(X,Y,CA,ERRC,CO,CCS,CB,BOUND,T,TI,NTIME,
!    &     WTIME,NITER,TOL,WTI)
!     IF(DEBUG) WRITE(*,*) 'EXIT UPDATE AND WRITING'
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  END TIME LOOP
      ENDDO
      
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      STOP 0
      END
