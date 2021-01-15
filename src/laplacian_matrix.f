!                 ***************************
                  SUBROUTINE LAPLACIAN_MATRIX
!                 ***************************
     & (DX,DY,BOUND,UPBOUND,DOBOUND,RIBOUND,LEBOUND,BOUNDINT,L)
!
!***********************************************************************
! NAVIER STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) COMPUTE LAPLACIAN MATRIX FOR SOLVE POISSON EQUATION.
!
!history  Sergio Castiblanco
!+        12/01/2021
!+        Translation for original Matlab implementation
!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| DX,DY     |-->| SPACE DIFFERENTIALS                                 |
!| BOUND     |-->| BOUNDARY INDICES                                    |
!| UPBOUND   |-->| TOP BOUNDARY INDICES                                |
!| DOBOUND   |-->| BOTTOM BOUNDARY INDICES                             |
!| RIBOUND   |-->| RIGHT BOUNDARY INDICES                              |
!| LEBOUND   |-->| LEFT BOUNDARY INDICES                               |
!| BOUNINT   |-->| INTERNAL BOUNDARY INDICES                           |
!| L         |<--| LAPLACIAN MATRIX                                    |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_NUMERICAL, ONLY:NX,NY,DEBUG,LU 
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION, INTENT(IN)                   :: DX,DY
      INTEGER, DIMENSION((2*NX+2*NY)-4), INTENT(IN)  :: BOUND
      INTEGER, DIMENSION(NX), INTENT(IN)             :: UPBOUND, DOBOUND
      INTEGER, DIMENSION(NY-2), INTENT(IN)           :: RIBOUND,LEBOUND
      INTEGER, DIMENSION(2*(NX-2)+2*(NY-4)), INTENT(IN) :: BOUNDINT
      DOUBLE PRECISION, DIMENSION(NX*NY,NX*NY), INTENT(OUT) :: L
!
! IN SUBROUTINE VARIABLES
!
      INTEGER :: I,J
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! COMPUTING LAPLACIAN OPERATOR MATRIX
!
      J=1
      DO I=1,NX*NY
        L(I,J) = -30*(1/(12*DX**2) + 1/(12*DY**2))
        IF (J.LT.NX*NY) THEN
          L(I,J+1) = (16/(12*DY**2))
        ENDIF
        IF (J.GT.1)THEN
          L(I,J-1) = (16/(12*DY**2))
        ENDIF
        IF (J.GT.2) THEN
          L(I,J-2) = -(1/(12*DY**2))
        ENDIF
        IF (J.LT.(NX*NY-1)) THEN
          L(I,J+2) = -(1/(12*DY**2))
        ENDIF
        IF (J.GT.NY) THEN
          L(I,J-NY) = (16/(12*DX**2))
        ENDIF
        IF (J.LT.(NX*NY-NY+1)) THEN
          L(I,J+NY) = (16/(12*DX**2))
        ENDIF
        IF (J.GT.2*NY) THEN
          L(I,J-2*NY) = -(1/(12*DX**2))
        ENDIF
        IF (J.LT.(NX*NY-2*NY+1)) THEN
          L(I,J+2*NY) = -(1/(12*DX**2))
        ENDIF
        J=J+1
      ENDDO
!
! ZEROING BOUNDARIES
!
      IF(DEBUG) WRITE(LU,*) 'ZEROING BOUNDARIES'
      L(BOUND,:) = 0.0D0
      L(BOUNDINT,:) = 0.0D0
!
! CHANGING VALUES FOR BOUNDINT
!
      IF(DEBUG) WRITE(LU,*) 'CHANGING VALUES FOR BOUNDINT'
      DO I = 1,SIZE(BOUNDINT)
        J = BOUNDINT(I)
        L(J,J) = -(2/(DX**2) + 2/(DY**2))
        L(J,J-NY) = 1/(DX**2)
        L(J,J+NY) = 1/(DX**2)
        L(J,J-1) = 1/(DY**2)
        L(J,J+1) = 1/(DY**2)
      ENDDO
!
! CHANGING VALUES FOR THE BOUNDARIES
!
      IF(DEBUG) WRITE(LU,*) 'CHANGING VALUES FOR BOUNDARIES'
      DO I = 1,SIZE(UPBOUND)
        J = UPBOUND(I)
        L(J,J) = 3D0/(2*DY)
        L(J,J+1) = -4/(2*DY)
        L(J,J+2) = 1/(2*DY)
      ENDDO
      DO I = 1,SIZE(DOBOUND)
        J = DOBOUND(I)
        L(J,J) = -3/(2*DY)
        L(J,J-1) = 4/(2*DY)
        L(J,J-2) = -1/(2*DY)
      ENDDO
      DO I = 1,SIZE(LEBOUND)
        J = LEBOUND(I)
        L(J,J) = -3/(2*DX)
        L(J,J+NY) = 4/(2*DX)
        L(J,J+2*NY) = -1/(2*DX)
      ENDDO
      DO I = 1,SIZE(RIBOUND)
        J = RIBOUND(I)
        L(J,J) = 3/(2*DX)
        L(J,J-NY) = -4/(2*DX)
        L(J,J-2*NY) = 1/(2*DX)
      ENDDO
!
      END SUBROUTINE LAPLACIAN_MATRIX
 
