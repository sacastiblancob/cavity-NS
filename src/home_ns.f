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
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  READING USER INPUT FILE
!
      NAMELIST /NSCONF/ NU,RHO,NX,NY,XMIN,XMAX,YMIN,YMAX,TO,TF,DEBUG
!
      OPEN(1010, FILE = "nsconf.nml", STATUS = 'OLD')
      READ(1010, NML = NSCONF)
      CLOSE(1010)
      IF(DEBUG) WRITE(LU,*) 'EXIT READING USER ENTRIES'
!
!  WRITING HEADERS
!
      IF(DEBUG) WRITE(LU,*) 'GOING INTO WRITE_HEADERS'
      CALL WRITE_HEADERS(LU)
      IF(DEBUG) WRITE(LU,*) 'EXIT WRITE_HEADERS'
!
!  ALLOCATING MEMORY AND SETTING INITIAL CONDITION
!
      IF(DEBUG) WRITE(LU,*) 'GOING INTO POINT_NS'
      CALL POINT_NS
      IF(DEBUG) WRITE(LU,*) 'EXIT POINT_NS'
!

      STOP 0
      END
