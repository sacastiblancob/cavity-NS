!                 ***************************
                  SUBROUTINE DIFFUSION_MATRIX
!                 ***************************
     & (DX,DY,DT,NU,BOUND,UPBOUND,DOBOUND,RIBOUND,LEBOUND,BOUNDINT,K,
     &    SX,SY)
!
!***********************************************************************
! NAVIER STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) COMPUTE DIFFUSION COEFFICIENTS AND MATRIX.
!
!history  Sergio Castiblanco
!+        12/01/2021
!+        Translation for original Matlab implementation
!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| DX,DY     |-->| SPACE DIFFERENTIALS                                 |
!| DT        |-->| TIME DIFFERENTIAL                                   |
!| NU        |-->| KINEMATIC VISCOSITY                                 |
!| BOUND     |-->| BOUNDARY INDICES                                    |
!| UPBOUND   |-->| TOP BOUNDARY INDICES                                |
!| DOBOUND   |-->| BOTTOM BOUNDARY INDICES                             |
!| RIBOUND   |-->| RIGHT BOUNDARY INDICES                              |
!| LEBOUND   |-->| LEFT BOUNDARY INDICES                               |
!| BOUNINT   |-->| INTERNAL BOUNDARY INDICES                           |
!| K         |<--| DIFFUSION MATRIX                                    |
!| SX,SY     |<->| DIFF COEFFICIENTS                                   |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_NUMERICAL, ONLY:NX,NY,DEBUG,LU 
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION, INTENT(IN)                   :: DX,DY,DT,NU
      INTEGER, DIMENSION((2*NX+2*NY)-4), INTENT(IN)  :: BOUND
      INTEGER, DIMENSION(NX), INTENT(IN)             :: UPBOUND, DOBOUND
      INTEGER, DIMENSION, INTENT(IN)                 :: RIBOUND,LEBOUND
      DOUBLE PRECISION, INTENT(INOUT)                :: SX,SY
      INTEGER, DIMENSION(2*(NX-2)+2*(NY-4)), INTENT(IN) :: BOUNDINT
      DOUBLE PRECISION, DIMENSION(NX*NY,NX*NY), INTENT(OUT) :: K
!
! IN SUBROUTINE VARIABLES
!
      INTEGER :: I,J
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! SETTING DIFFUSSION STABILITY PARAMETERS
!
      IF(DEBUG) WRITE(LU,*) 'COMPUTING SX AND SY'
      SX = ((DT*NU)/(DX**2))
      SY = ((DT*NU)/(DY**2))
      IF(DEBUG) WRITE(LU,*) 'END COMPUTING SX AND SY'
!
! COMPUTING STTIFFNESS DIFFUSION MATRIX
!
      J=1
      DO I=1,NX,NY
        K(I,J) = -(1 + 30*SX/12 + 30*SY/12)
        IF (J.LT.NX*NY) THEN
          K(I,J+1) = (16*SY/12)
        ENDIF
        IF (J.GT.1)THEN
          K(I,J-1) = (16*SY/12)
        ENDIF
        IF (J.GT.2) THEN
          K(I,J-2) = -(SY/12)
        ENDIF
        IF (J.LT.(NX*NY-1) THEN
          K(I,J+2) = -(SY/12)
        ENDIF
        IF (J.GT.NY) THEN
          K(I,J-NY) = (16*SX/12)
        ENDIF
        IF (J.LT.(NX*NY-NY+1)) THEN
          K(I,J+NY) = (16*SX/12)
        ENDIF
        IF (J.GT.2*NY) THEN
          K(I,J-2*NY) = -(SX/12)
        ENDIF
        IF (J.LT.(NX*NY-2*NY+1) THEN
          K(I,J+2*NY) = -(SX/12)
        ENDIF
        J=J+1
      ENDDO
!
! ZEROING BOUNDARIES
!
      IF(DEBUG) WRITE(LU,*) 'ZEROING BOUNDARIES'
      K(BOUND,:) = 0.0D0
      K(BOUNDINT,:) = 0.0D0
!
! CHANGING VALUES FOR BOUNDINT
!
      IF(DEBUG) WRITE(LU,*) 'CHANGING VALUES FOR BOUNDINT'
      DO I = 1,SIZE(BOUNDINT)
        J = BOUNDINT(I)
        K(J,J) = -(1 + 2*SX + 2*SY)
        K(J,J-NY) = SX
        K(J,J+NY) = SX
        K(J,J-1) = SY
        K(J,J+1) = SY
      ENDDO
!
! CHANGING VALUES FOR THE BOUNDARIES
!
      IF(DEBUG) WRITE(LU,*) 'CHANGING VALUES FOR BOUNDARIES'
      DO I = 1,SIZE(UPBOUND)
        J = UPBOUND(I)
        K(J,J) = 1.0D0
      ENDDO
      DO I = 1,SIZE(DOBOUND)
        J = DOBOUND(I)
        K(J,J) = 1.0D0
      ENDDO
      DO I = 1,SIZE(LEBOUND)
        J = LEBOUND(I)
        K(J,J) = 1.0D0
      ENDDO
      DO I = 1,SIZE(RIBOUND)
        J = RIBOUND(I)
        K(J,J) = 1.0D0
      ENDDO
!
! DISPLAYING COMPUTED SX AND SY
!
      WRITE(LU,*) REPEAT('~',72)
      WRITE(LU,*) 'SX VALUE: ', SX
      WRITE(LU,*) 'SY VALUE: ', SY
      WRITE(LU,*) REPEAT('~',72)
!
      END SUBROUTINE DIFFUSION_MATRIX
 
