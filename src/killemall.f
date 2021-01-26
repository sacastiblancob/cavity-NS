!                      ********************
                       SUBROUTINE KILLEMALL
!                      ********************
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!***********************************************************************
! NAVIER-STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) DEALLOCATING.
!
!history  Sergio Castiblanco
!+        26/01/2021
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE CSC_STORAGE
      USE DECLARATIONS_PHYSIC
      USE DECLARATIONS_NUMERICAL
!
      IMPLICIT NONE
!    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ALLOCATING GRID
!
      IF(DEBUG) WRITE(*,*) 'DEALLOCATING GRID, X AND Y'
      DEALLOCATE(X)
      DEALLOCATE(Y)
!
! DEALLOCATE BOUNDARY POSITION VECTORS
!
      IF(DEBUG) WRITE(*,*) 'DEALLOCATING BOUNDARY POSITIONAL VECTORS'
      DEALLOCATE(POS)
      DEALLOCATE(MPOS)
!
!  EXTERNAL BOUNDARIES
      DEALLOCATE(BOUND)
      DEALLOCATE(UPBOUND)
      DEALLOCATE(DOBOUND)
      DEALLOCATE(LEBOUND)
      DEALLOCATE(RIBOUND)
!  INTERNAL BOUNDARIES
!
      DEALLOCATE(BOUNDINT)
      DEALLOCATE(UPBOUNDI)
      DEALLOCATE(DOBOUNDI)
      DEALLOCATE(RIBOUNDI)
      DEALLOCATE(LEBOUNDI)
!
! DEALLOCATE VELOCITIES VECTORS AND PRESURE
!
      IF(DEBUG) WRITE(*,*) 'DEALLOCATING VELOCITIES VECTORS'
      ! INITIAL CONDITION AND PREVIOUS STEP RESULT
      DEALLOCATE(UO)
      DEALLOCATE(VO)
      ! VELOCITIES AFTER FIRST FRACTIONAL STEP
      DEALLOCATE(UP)
      DEALLOCATE(VP)
      ! VELOCITIES AFTER SECOND FRACTIONAL STEP
      DEALLOCATE(UPP)
      DEALLOCATE(VPP)
      ! VELOCITIES AFTER THIRD FRACTIONAL STEP
      DEALLOCATE(U)
      DEALLOCATE(V)
      ! PRESURE
      DEALLOCATE(P)
      ! PRESURE FOR PRINTING
      DEALLOCATE(PP)
      ! SINGULAR VECTORS OF L
      DEALLOCATE(UL)
      DEALLOCATE(VL)
      ! REGULARIZATION MATRIX
      DEALLOCATE(RM)
!
! DEALLOCATE DERIVATIVES VECTORS
!
      !FOR NON-LINEAR ADVECTION STEP
      DEALLOCATE(DUDX)
      DEALLOCATE(DUDY)
      DEALLOCATE(DVDX)
      DEALLOCATE(DVDY)
      !FOR RHS OF POISSON EQUATION
      DEALLOCATE(DUPDX)
      DEALLOCATE(DVPDY)
      !FOR PRESURE
      DEALLOCATE(DPDX)
      DEALLOCATE(DPDY)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE KILLEMALL



