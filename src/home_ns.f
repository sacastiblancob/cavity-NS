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
      STOP 0
      END
