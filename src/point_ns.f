!                       *******************
                        SUBROUTINE POINT_NS
!                       *******************
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE DECLARATIONS_PHYSIC
      USE DECLARATIONS_NUMERICAL
!
      IMPLICIT NONE
!    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IN SUBROUTINE VARIABLES
      INTEGER :: I,J,M
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ALLOCATING GRID
!
      IF(DEBUG) WRITE(LU,*) 'ALLOCATING GRID, X AND Y'
      ALLOCATE(X(NX))
      X = 0D0
      ALLOCATE(Y(NY))
      Y = 0D0
!
! ALLOCATE BOUNDARY POSITION VECTORS
!
      IF(DEBUG) WRITE(LU,*) 'ALLOCATING BOUNDARY POSITIONAL VECTORS'
      ALLOCATE(POS(NX*NY))
      POS = 1
      ALLOCATE(MPOS(NY,NX))
      MPOS = 1
      ALLOCATE(BOUND((2*NX + 2*NY)-4))
      BOUND = 1
      ALLOCATE(UPBOUND(NX))
      UPBOUND = 1
      ALLOCATE(DOBOUND(NX))
      DOBOUND = 1
      ALLOCATE(RIBOUND(NY-2))
      RIBOUND = 1
      ALLOCATE(LEBOUND(NY-2))
      LEBOUND = 1
      ALLOCATE(BOUNDINT(2*(NX-2) + 2*(NY-4)))
      BOUNDINT = 1
!
! ALLOCATE STIFFNESS AND LAPLACE MATRICES
!
      IF(DEBUG) WRITE(LU,*) 'ALLOCATING STIFFNESS AND LAPLACE MATRICES'
      ALLOCATE(K(NX*NY,NX*NY))
      K = 0D0
      ALLOCATE(L(NX*NY,NX*NY))
      L = 0D0
!
! ALLOCATE VELOCITIES VECTORS
!
      IF(DEBUG) WRITE(LU,*) 'ALLOCATING VELOCITY VECTORS'
      ALLOCATE(UO(NX*NY))
      UO = 0D0
      ALLOCATE(VO(NX*NY))
      VO = 0D0
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  GRID - GRID - GRID - GRID
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! COMPUTING DIFFERENTIALS
!
      DX = (XMAX - XMIN)/(NX-1)
      DY = (YMAX - YMIN)/(NY-1)
!
! COMPUTING X AND Y
!
      IF(DEBUG) WRITE(LU,*) 'COMPUTING X'
      X(1) = XMIN
      DO I=2,NX
        X(I) = X(I-1) + DX
      ENDDO
!
      IF(DEBUG) WRITE(LU,*) 'COMPUTING Y'
      Y(1) = YMIN
      DO I=2,NY
        Y(I) = Y(I-1) + DY
      ENDDO
!
! ENUMERATION AND INDEXING
!
      IF(DEBUG) WRITE(LU,*) 'COMPUTING POS (ENUMERATION)'
      DO I=1,NX*NY
        POS(I) = I
      ENDDO
      M=1
      DO J=1,NX
        DO I=1,NY
          MPOS(I,J) = POS(M)
          M = M+1
        ENDDO
      ENDDO
!
! BOUNDARIES INDEXING
!
      IF(DEBUG) WRITE(LU,*) 'COMPUTING BOUNDARIES INDEXING'
!     UP BOUNDARY
      UPBOUND = MPOS(1,:)
!
!     BOTTOM BOUNDARY
      DOBOUND = MPOS(NY,:)
!
!     RIGHT BOUNDARY
      RIBOUND = MPOS(2:NY-1,NX)
!
!     LEFT BOUNDARY
      LEBOUND = MPOS(2:NY-1,1)
!
!     ALL BOUNDARIES TOGETHER
      BOUND = [UPBOUND,DOBOUND,RIBOUND,LEBOUND]
!
!     INTERNAL BOUNDARY
      BOUNDINT = [MPOS(2,2:NX-1),MPOS(NY-1,2:NX-1),
     &    MPOS(3:NY-2,NX-1), MPOS(3:NY-2,2)]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! END GRID - END GRID - END GRID
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TIME - TIME - TIME - TIME
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! SETTING INITIAL CONDITION
!
      IF(DEBUG) WRITE(LU,*) 'SETTING INITIAL CONDITION'
      UO(1:NX) = 1.0D0
!
! SETTING TIME PARAMETERS
!
      IF(DEBUG) WRITE(LU,*) 'COMPUTING DT AND NT'
      DT = MIN(DX,DY)*CFL/ABS(UO(1))
      NT = INT(FLOOR((TF-TO)/DT))
!
! ALLOCATE TIME VECTOR
!
      IF(DEBUG) WRITE(LU,*) 'ALLOCATING AND COMPUTING TIME VECTOR'
      ALLOCATE(T(NT))
!
! FILLING TIME VECTOR
!
      T(1) = TO
      DO I=2,NT
        T(I) = T(I-1) + DT
      ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! END TIME - END TIME - END TIME
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE POINT_NS





