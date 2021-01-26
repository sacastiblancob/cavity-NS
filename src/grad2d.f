!                    *****************
                     SUBROUTINE GRAD2D
!                    *****************
     & (DUDX,DVDY,UO,VO,UPBOUND,DOBOUND,LEBOUND,RIBOUND)

!
!***********************************************************************
! 2D-NAVIER STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!BRIEF    1) COMPUTE DIFFUSION COEFFICIENTS AND MATRIX.
!
!HISTORY  SERGIO CASTIBLANCO
!+        25/01/2021
!+        TRANSLATION FOR ORIGINAL MATLAB IMPLEMENTATION
!
!!
!!
!! POISSON EQUATION
!!
!!   DP2/D2X + DP2/D2Y = -(RHO/DT)*(DU*/DX + DV*/DY)
!!
!! GRID AND WHAT BOUND, UPBOUNDIN, DOBOUNDIN, LEBOUNDIN, AND RIBOUNDIN
!! MEANS
!!
!! ENUMARTION USED IN THE GRID
!!
!!    25   26   27   28   29   30
!!
!!    19   20   21   22   23   24
!!
!!    13   14   15   16   17   18
!!
!!    7    8    9    10   11   12
!!
!!    1    2    3    4    5    6
!!
!! UPBOUND(UB), DOBOUND(DB), LEBOUND(LB), RIBOUND(RB), BOUNDINT(BI)
!! INTERNAL NODES(I)
!!
!!    ULB  UB   UB   UB   UB   URB
!!
!!    LB   I    I    I    I    RB
!!
!!    LB   I    I    I    I    RB
!!
!!    LB   I    I    I    I    RB
!!
!!    DLB  DB   DB   DB   DB   DRB
!!
!! DUL --> TOP-LEFT NODE
!! DDL --> BOTTOM-LEFT NODE
!! DDR --> BOTTOM-RIGHT NODE
!! DUR --> TOP-RIGHT NODE
!!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| DUDX      |<--| D(UO)/DX VECTOR                                     |
!| DVDY      |<--| D(VO)/DY VECTOR                                     |
!| UO        |-->| UO                                                  |
!| VO        |-->| VO                                                  |
!| UPBOUND   |-->| TOP BOUNDARY INDICES                                |
!| DOBOUND   |-->| BOTTOM BOUNDARY INDICES                             |
!| LEBOUND   |-->| LEFT BOUNDARY INDICES                               |
!| RIBOUND   |-->| RIGHT BOUNDARY INDICES                              |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_NUMERICAL, ONLY: NX,NY,DX,DY
!
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION,DIMENSION(NX*NY),INTENT(IN) :: UO,VO
      DOUBLE PRECISION,DIMENSION(NX*NY),INTENT(OUT) :: DUDX,DVDY
      INTEGER, DIMENSION(NX-2), INTENT(IN) :: UPBOUND, DOBOUND
      INTEGER, DIMENSION(NY-2), INTENT(IN) :: LEBOUND, RIBOUND
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: I
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DO I = 1,NX*NY
!
!  DERIVATIVES IN THE CORNERS
!
        IF(I.EQ.1) THEN
          DUDX(I) = (1/DX)*(UO(I+1) - UO(I))
          DVDY(I) = (1/DY)*(VO(I+NX) - VO(I))
        ELSEIF(I.EQ.NX) THEN
          DUDX(I) = (1/DX)*(UO(I) - UO(I-1))
          DVDY(I) = (1/DY)*(VO(2*I) - VO(I))
        ELSEIF(I.EQ.(NX*NY-NX+1)) THEN
          DUDX(I) = (1/DX)*(UO(NX*NY-NX+2) - UO(NX*NY-NX+1))
          DVDY(I) = (1/DY)*(VO(NX*NY-NX+1) - VO(NX*NY-2*NX+1))
        ELSEIF(I.EQ.(NX*NY)) THEN
          DUDX(I) = (1/DX)*(UO(NX*NY) - UO(NX*NY-1))
          DVDY(I) = (1/DY)*(VO(NX*NY) - VO(NX*NY-NX))
!
!  DERIVATIVES IN THE BOUNDARIES
!
        ELSEIF(ANY(I.EQ.UPBOUND)) THEN
          DUDX(I) = (1/(2*DX))*(UO(I+1) - UO(I-1));
          DVDY(I) = (1/DY)*(VO(I) - VO(I-NX));
        ELSEIF(ANY(I.EQ.DOBOUND)) THEN
          DUDX(I) = (1/(2*DX))*(UO(I+1) - UO(I-1));
          DVDY(I) = (1/DY)*(VO(I+NX) - VO(I));
        ELSEIF(ANY(I.EQ.LEBOUND)) THEN
          DUDX(I) = (1/DX)*(UO(I+1) - UO(I));
          DVDY(I) = (1/(2*DY))*(VO(I+NX) - VO(I-NX));
        ELSEIF(ANY(I.EQ.RIBOUND)) THEN
          DUDX(I) = (1/DX)*(UO(I) - UO(I-1));
          DVDY(I) = (1/(2*DY))*(VO(I+NX) - VO(I-NX));
!
!  DERIVATIVES IN THE INTERNAL NODES
!
        ELSE
          DUDX(I) = (1/(2*DX))*(UO(I+1) - UO(I-1));
          DVDY(I) = (1/(2*DY))*(VO(I+NX) - VO(I-NX));
        ENDIF
      ENDDO



!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE GRAD2D



