!                    ******************
                     SUBROUTINE DIVER2D
!                    ******************
     & (DUDX,DUDY,DVDX,DVDY,UO,VO,BOUND,UPBOUNDI,DOBOUNDI,LEBOUNDI,
     &  RIBOUNDI)

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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!
!! THIS FUNCTION COMPUTES THE TERMS DU/DX, DU/DY, DV/DX, DV/DY WITH 
!! UPWIND SECOND ORDER SCHEME (IF POSIBLE)
!!
!! NONLINEAL ADVECTION EQUATION
!!
!!   DU/DT + U DU/DX + V DU/DY = 0
!!
!!   DV/DT + U DV/DX + V DV/DY = 0
!!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! ENTRIES
!!   UO, VECTOR WITH VELOCITIES IN X DIRECTION
!!   UV, VECTOR WITH VELOCITIES IN Y DIRECTION
!!   DX, DIFFERENTIAL IN X
!!   DY, DIFFERENTIAL IN Y
!!   N, NUMBER OF NODES IN X DIRECTION
!!   M, NUMBER OF NODES IN Y DIRECTION
!!   BOUND, INDICES OF THE BOUNDARY NODES (B)
!!   UPBOUNDIN, INDICES WITH THE TOP INTERNAL BOUNDARY (U)
!!   DOBOUNDIN, INDICES WITH THE BOTTOM INTERNAL BOUNDARY (D)
!!   LEBOUNDIN, INDICES WITH THE LEFT INTERNAL BOUNDARY (L)
!!   RIBOUNDIN, INDICES WITH THE RIGHT INTERNAL BOUNDARY (R)
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
!! BOUND(B), UPBOUNDIN(U), DOBOUNDIN(D), LEBOUNDIN(L), RIBOUNDIN(R),
!! INTERNAL NODES(I)
!!
!!    B    B    B    B    B    B
!!
!!    B    UL   U    U    UR   B
!!
!!    B    L    I    I    R    B
!!
!!    B    DL   D    D    DR   B
!!
!!    B    B    B    B    B    B
!!
!! UL --> TOP-LEFT INTERNAL
!! DL --> BOTTOM-LEFT INTERNAL
!! DR --> BOTTOM-RIGHT INTERNAL
!! UR --> TOP-RIGHT INTERNAL
!!
!! THE RESULTS IN THE BOUNDARY ARE SET TO ZERO, BECAUSE THEY ARE NO 
!! NEEDED TO SOLVE NON-LINEAL ADVECTION EQUATION (DUE TO THE BOUNDARY 
!! CONDITIONS) IN OTHER WORDS, THE DIVERGENCE IS ONLY COMPUTED IN THE 
!! INTERNAL NODES.!
!!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| DUDX      |<--| DU/DX VECTOR                                        |
!| DUDY      |<--| DU/DY VECTOR                                        |
!| DVDX      |<--| DV/DX VECTOR                                        |
!| DVDY      |<--| DV/DY VECTOR                                        |
!| UO        |-->| INITIAL/PREVIOUS VELOCITY IN X                      |
!| VO        |-->| INITIAL/PREVIOUS VELOCITY IN Y                      |
!| BOUND     |-->| BOUNDARY INDICES                                    |
!| UPBOUNDI  |-->| TOP INTERNAL BOUNDARY INDICES                       |
!| DOBOUNDI  |-->| BOTTOM INTERNAL BOUNDARY INDICES                    |
!| LEBOUNDI  |-->| LEFT INTERNAL BOUNDARY INDICES                      |
!| RIBOUNDI  |-->| RIGHT INTERNAL BOUNDARY INDICES                     |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_NUMERICAL, ONLY: NX,NY,DX,DY
!
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION,DIMENSION(NX*NY),INTENT(IN) :: UO,VO
      DOUBLE PRECISION,DIMENSION(NX*NY),INTENT(OUT) :: DUDX,DUDY,
     &    DVDX,DVDY
      INTEGER, DIMENSION(2*NX+2*(NY-2)), INTENT(IN) :: BOUND
      INTEGER, DIMENSION(NX-4), INTENT(IN) :: UPBOUNDI, DOBOUNDI
      INTEGER, DIMENSION(NY-4), INTENT(IN) :: LEBOUNDI, RIBOUNDI
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: I
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  DERIVATIVES IN THE INTERNAL NODES
!
      DO I=1,NX*NY
        IF((UO(I).GT.0).AND.(VO(I).GT.0)) THEN
          IF(ANY(I.EQ.BOUND)) THEN
            CONTINUE
          ELSEIF(ANY(I.EQ.DOBOUNDI).OR.(I.EQ.(2*NX-1)))THEN
            DUDX(I) = (1/(2*DX))*(3*UO(I) - 4*UO(I-1) + UO(I-2));
            DUDY(I) = (1/DY)*(UO(I) - UO(I-NX));
            DVDX(I) = (1/(2*DX))*(3*VO(I) - 4*VO(I-1) + VO(I-2));
            DVDY(I) = (1/DY)*(VO(I) - VO(I-NX));
          ELSEIF(ANY(I.EQ.LEBOUNDI).OR.(I.EQ.(NX*NY-2*NX+2))) THEN
            DUDX(I) = (1/DX)*(UO(I) - UO(I-1));
            DUDY(I) = (1/(2*DY))*(3*UO(I) - 4*UO(I-NX) + UO(I-2*NX));
            DVDX(I) = (1/DX)*(VO(I) - VO(I-1));
            DVDY(I) = (1/(2*DY))*(3*VO(I) - 4*VO(I-NX) + VO(I-2*NX));
          ELSEIF(I.EQ.(NX+2)) THEN
            DUDX(I) = (1/DX)*(UO(I) - UO(I-1));
            DUDY(I) = (1/DY)*(UO(I) - UO(I-NX));
            DVDX(I) = (1/DX)*(VO(I) - VO(I-1));
            DVDY(I) = (1/DY)*(VO(I) - VO(I-NX));
          ELSE
            DUDX(I) = (1/(2*DX))*(3*UO(I) - 4*UO(I-1) + UO(I-2));
            DUDY(I) = (1/(2*DY))*(3*UO(I) - 4*UO(I-NX) + UO(I-2*NX));
            DVDX(I) = (1/(2*DX))*(3*VO(I) - 4*VO(I-1) + VO(I-2));
            DVDY(I) = (1/(2*DY))*(3*VO(I) - 4*VO(I-NX) + VO(I-2*NX));
          ENDIF
        ELSEIF((UO(I).LE.0).AND.(VO(I).GT.0)) THEN
          IF(ANY(I.EQ.BOUND)) THEN
            CONTINUE
          ELSEIF(ANY(I.EQ.DOBOUNDI).OR.(I.EQ.(NX+2))) THEN
            DUDX(I) = (1/(2*DX))*(-3*UO(I) + 4*UO(I+1) - UO(I+2));
            DUDY(I) = (1/DY)*(UO(I) - UO(I-NX));
            DVDX(I) = (1/(2*DX))*(-3*VO(I) + 4*VO(I+1) - VO(I+2));
            DVDY(I) = (1/DY)*(VO(I) - VO(I-NX));
          ELSEIF(ANY(I.EQ.RIBOUNDI).OR.(I.EQ.(NX*NY-NX-1))) THEN
            DUDX(I) = (1/DX)*(-UO(I) + UO(I+1));
            DUDY(I) = (1/(2*DY))*(3*UO(I) - 4*UO(I-NX) + UO(I-2*NX));
            DVDX(I) = (1/DX)*(-VO(I) + VO(I+1));
            DVDY(I) = (1/(2*DY))*(3*VO(I) - 4*VO(I-NX) + VO(I-2*NX));
          ELSEIF(I.EQ.(2*NX-1)) THEN
            DUDX(I) = (1/DX)*(-UO(I) + UO(I+1));
            DUDY(I) = (1/DY)*(UO(I) - UO(I-NX));
            DVDX(I) = (1/DX)*(-VO(I) + VO(I+1));
            DVDY(I) = (1/DY)*(VO(I) - VO(I-NX));
          ELSE
            DUDX(I) = (1/(2*DX))*(-3*UO(I) + 4*UO(I+1) - UO(I+2));
            DUDY(I) = (1/(2*DY))*(3*UO(I) - 4*UO(I-NX) + UO(I-2*NX));
            DVDX(I) = (1/(2*DX))*(-3*VO(I) + 4*VO(I+1) - VO(I+2));
            DVDY(I) = (1/(2*DY))*(3*VO(I) - 4*VO(I-NX) + VO(I-2*NX));
          ENDIF
        ELSEIF((UO(I).GT.0).AND.(VO(I).LE.0)) THEN
          IF(ANY(I.EQ.BOUND)) THEN
            CONTINUE
          ELSEIF(ANY(I.EQ.UPBOUNDI).OR.(I.EQ.(NX*NY-NX-1))) THEN
            DUDX(I) = (1/(2*DX))*(3*UO(I) - 4*UO(I-1) + UO(I-2));
            DUDY(I) = (1/DY)*(-UO(I) + UO(I+NX));
            DVDX(I) = (1/(2*DX))*(3*VO(I) - 4*VO(I-1) + VO(I-2));
            DVDY(I) = (1/DY)*(-VO(I) + VO(I+NX));
          ELSEIF(ANY(I.EQ.LEBOUNDI).OR.(I.EQ.(NX+2))) THEN
            DUDX(I) = (1/DX)*(UO(I) - UO(I-1));
            DUDY(I) = (1/(2*DY))*(-3*UO(I) + 4*UO(I+NX) - UO(I+2*NX));
            DVDX(I) = (1/DX)*(VO(I) - VO(I-1));
            DVDY(I) = (1/(2*DY))*(-3*VO(I) + 4*VO(I+NX) - VO(I+2*NX));
          ELSEIF(I.EQ.(NX*NY-2*NX+2)) THEN
            DUDX(I) = (1/DX)*(UO(I) - UO(I-1));
            DUDY(I) = (1/DY)*(-UO(I) + UO(I+NX));
            DVDX(I) = (1/DX)*(VO(I) - VO(I-1));
            DVDY(I) = (1/DY)*(-VO(I) + VO(I+NX));
          ELSE
            DUDX(I) = (1/(2*DX))*(3*UO(I) - 4*UO(I-1) + UO(I-2));
            DUDY(I) = (1/(2*DY))*(-3*UO(I) + 4*UO(I+NX) - UO(I+2*NX));
            DVDX(I) = (1/(2*DX))*(3*VO(I) - 4*VO(I-1) + VO(I-2));
            DVDY(I) = (1/(2*DY))*(-3*VO(I) + 4*VO(I+NX) - VO(I+2*NX));
          ENDIF
        ELSE
          IF(ANY(I.EQ.BOUND)) THEN
            CONTINUE
          ELSEIF(ANY(I.EQ.UPBOUNDI).OR.(I.EQ.(NX*NY-2*NX+2))) THEN
            DUDX(I) = (1/(2*DX))*(-3*UO(I) + 4*UO(I+1) - UO(I+2));
            DUDY(I) = (1/DY)*(-UO(I) + UO(I+NX));
            DVDX(I) = (1/(2*DX))*(-3*VO(I) + 4*VO(I+1) - VO(I+2));
            DVDY(I) = (1/DY)*(-VO(I) + VO(I+NX));
          ELSEIF(ANY(I.EQ.RIBOUNDI).OR.(I.EQ.(2*NX-1))) THEN
            DUDX(I) = (1/DX)*(-UO(I) + UO(I+1));
            DUDY(I) = (1/(2*DY))*(-3*UO(I) + 4*UO(I+NX) - UO(I+2*NX));
            DVDX(I) = (1/DX)*(-VO(I) + VO(I+1));
            DVDY(I) = (1/(2*DY))*(-3*VO(I) + 4*VO(I+NX) - VO(I+2*NX));
          ELSEIF(I.EQ.(NX*NY-NX-1)) THEN
            DUDX(I) = (1/DX)*(-UO(I) + UO(I+1));
            DUDY(I) = (1/DY)*(-UO(I) + UO(I+NX));
            DVDX(I) = (1/DX)*(-VO(I) + VO(I+1));
            DVDY(I) = (1/DY)*(-VO(I) + VO(I+NX));
          ELSE
            DUDX(I) = (1/(2*DX))*(-3*UO(I) + 4*UO(I+1) - UO(I+2));
            DUDY(I) = (1/(2*DY))*(-3*UO(I) + 4*UO(I+NX) - UO(I+2*NX));
            DVDX(I) = (1/(2*DX))*(-3*VO(I) + 4*VO(I+1) - VO(I+2));
            DVDY(I) = (1/(2*DY))*(-3*VO(I) + 4*VO(I+NX) - VO(I+2*NX));
          ENDIF
        ENDIF
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE DIVER2D


