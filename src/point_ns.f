!                       *******************
                        SUBROUTINE POINT_NS
!                       *******************
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!***********************************************************************
! NAVIER-STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) ALLOCATING AND COMPUTING INITIAL CONSTANTS AND SETUP.
!
!history  Sergio Castiblanco
!+        16/01/2021
!+        Translation for original Matlab implementation
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
! IN SUBROUTINE VARIABLES
      INTEGER :: I,J,K,M
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ALLOCATING GRID
!
      IF(DEBUG) WRITE(*,*) 'ALLOCATING GRID, X AND Y'
      ALLOCATE(X(NX))
      X = 0D0
      ALLOCATE(Y(NY))
      Y = 0D0
!
! ALLOCATE BOUNDARY POSITION VECTORS
!
      IF(DEBUG) WRITE(*,*) 'ALLOCATING BOUNDARY POSITIONAL VECTORS'
      ALLOCATE(POS(NX*NY))
      POS = 1
      ALLOCATE(MPOS(NY,NX))
      MPOS = 1
!
!  EXTERNAL BOUNDARIES
      ALLOCATE(BOUND(2*NX + 2*(NY-2)))
      BOUND = 1
      ALLOCATE(UPBOUND(NX-2))
      UPBOUND = 1
      ALLOCATE(DOBOUND(NX-2))
      DOBOUND = 1
      ALLOCATE(LEBOUND(NY-2))
      LEBOUND = 1
      ALLOCATE(RIBOUND(NY-2))
      RIBOUND = 1
!  INTERNAL BOUNDARIES
!
      ALLOCATE(BOUNDINT(2*(NX-2) + 2*(NY-4)))
      BOUNDINT = 1
      ALLOCATE(UPBOUNDI(NX-4))
      UPBOUNDI = 1
      ALLOCATE(DOBOUNDI(NX-4))
      DOBOUNDI = 1
      ALLOCATE(RIBOUNDI(NY-4))
      RIBOUNDI = 1
      ALLOCATE(LEBOUNDI(NY-4))
      LEBOUNDI = 1
!
! ALLOCATE VELOCITIES VECTORS AND PRESURE
!
      IF(DEBUG) WRITE(*,*) 'ALLOCATING VELOCITIES VECTORS'
      ! INITIAL CONDITION AND PREVIOUS STEP RESULT
      ALLOCATE(UO(NX*NY))
      UO = 0D0
      ALLOCATE(VO(NX*NY))
      VO = 0D0
      ! VELOCITIES AFTER FIRST FRACTIONAL STEP
      ALLOCATE(UP(NX*NY))
      UP = 0D0
      ALLOCATE(VP(NX*NY))
      VP = 0D0
      ! VELOCITIES AFTER SECOND FRACTIONAL STEP
      ALLOCATE(UPP(NX*NY))
      UPP = 0D0
      ALLOCATE(VPP(NX*NY))
      VPP = 0D0
      ! VELOCITIES AFTER THIRD FRACTIONAL STEP
      ALLOCATE(U(NX*NY))
      U = 0D0
      ALLOCATE(V(NX*NY))
      V = 0D0
      ! PRESURE
      ALLOCATE(P(NX*NY))
      P = 0D0
      ! PRESURE FOR PRINTING
      ALLOCATE(PP(NX*NY))
      PP = 0D0
      ! SINGULAR VECTORS OF L
      ALLOCATE(UL(NX*NY))
      UL = 0D0
      ALLOCATE(VL(NX*NY))
      VL = 0D0
      ! REGULARIZATION MATRIX
      ALLOCATE(RM(NX*NY,NX*NY))
      RM = 0D0
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
      IF(DEBUG) WRITE(*,*) 'COMPUTING X'
      X(1) = XMIN
      DO I=2,NX
        X(I) = X(I-1) + DX
      ENDDO
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING Y'
      Y(1) = YMIN
      DO I=2,NY
        Y(I) = Y(I-1) + DY
      ENDDO
!
! ENUMERATION AND INDEXING
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING POS (ENUMERATION)'
      DO I=1,NX*NY
        POS(I) = I
      ENDDO
      M=1
      DO I = NY,1,-1
        DO J = 1,NX
          MPOS(I,J) = POS(M)
          M = M+1
        ENDDO
      ENDDO
      ! ! !DO I = 1,NY
      ! ! !  WRITE(*,*) 'MPOS',I,':',MPOS(I,:)
      ! ! !ENDDO
!
! BOUNDARIES INDEXING
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING BOUNDARIES INDEXING'
!     INTERNAL BOUNDARIES
      K=1 
      DO I = (NX*NY-2*NX+3),(NX*NY-NX-2)
        UPBOUNDI(K) = I
        K=K+1
      ENDDO
      K=1
      DO I =(NX+3),(2*NX-2)
        DOBOUNDI(K) = I
        K = K+1
      ENDDO
      K=1
      DO I=(2*NX+2),(NX*NY-3*NX+2),NX
        LEBOUNDI(K) = I
        K=K+1
      ENDDO
      K=1
      DO I=(3*NX-1),(NX*NY-2*NX-1),NX
        RIBOUNDI(K) = I
        K=K+1
      ENDDO
      BOUNDINT = [MPOS(2,2:NX-1),MPOS(NY-1,2:NX-1),MPOS(3:NY-2,2),
     &    MPOS(3:NY-2,NX-1)]

!     UP BOUNDARY
      K=1
      DO I=(NX*NY-NX+2),(NX*NY-1)
        UPBOUND(K) = I
        K = K+1
      ENDDO
!
!     BOTTOM BOUNDARY
      K=1
      DO I=(2),(NX-1)
        DOBOUND(K) = I
        K=K+1
      ENDDO
!
!     RIGHT BOUNDARY
      K=1
      DO I=(2*NX),(NX*NY-NX),NX
        RIBOUND(K) = I
        K=K+1
      ENDDO
!
!     LEFT BOUNDARY
      K=1
      DO I=(NX+1),(NX*NY-2*NX+1),NX
        LEBOUND(K) = I
        K=K+1
      ENDDO
!
!     ALL BOUNDARIES TOGETHER
      BOUND = [MPOS(NY,:),MPOS(1,:),MPOS(2:NY-1,1),MPOS(2:NY-1,NX)]
!
!  SORTING BOUNDARY INDICES VECTORS
      CALL QUICKSORT(BOUND,1,SIZE(BOUND))
      CALL QUICKSORT(LEBOUND,1,SIZE(LEBOUND))
      CALL QUICKSORT(RIBOUND,1,SIZE(RIBOUND))
      CALL QUICKSORT(UPBOUND,1,SIZE(UPBOUND))
      CALL QUICKSORT(DOBOUND,1,SIZE(DOBOUND))
      CALL QUICKSORT(LEBOUNDI,1,SIZE(LEBOUNDI))
      CALL QUICKSORT(RIBOUNDI,1,SIZE(RIBOUNDI))
      CALL QUICKSORT(UPBOUNDI,1,SIZE(UPBOUNDI))
      CALL QUICKSORT(DOBOUNDI,1,SIZE(DOBOUNDI))
      CALL QUICKSORT(BOUNDINT,1,SIZE(BOUNDINT))
      ! ! !WRITE(*,*) 'BOUND ',BOUND
      ! ! !WRITE(*,*) 'LEBOUND ',LEBOUND
      ! ! !WRITE(*,*) 'RIBOUND ',RIBOUND
      ! ! !WRITE(*,*) 'UPBOUND ',UPBOUND
      ! ! !WRITE(*,*) 'DOBOUND ',DOBOUND
      ! ! !WRITE(*,*) 'LEBOUNDI ',LEBOUNDI
      ! ! !WRITE(*,*) 'RIBOUNDI ',RIBOUNDI
      ! ! !WRITE(*,*) 'UPBOUNDI ',UPBOUNDI
      ! ! !WRITE(*,*) 'DOBOUNDI ',DOBOUNDI
      ! ! !WRITE(*,*) 'BOUNDINT ',BOUNDINT
!      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! END GRID - END GRID - END GRID
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! INITIAL CONDITION
!
      DO I=1,SIZE(UPBOUND)
        UO(UPBOUND(I)) = 1.0D0
      ENDDO
      UO(NX*NY) = 1.0D0
      UO(NX*NY-NX+1) = 1.0D0

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TIME - TIME - TIME - TIME
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! SETTING TIME PARAMETERS
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING DT AND NT'
!  WRITING CHANGE IN DT IF ANY
      IF(CFL.GT.1.2D0) THEN
        CFL = 1.2D0
      ENDIF
      DT = (MIN(DX,DY)*CFL)/MAXVAL(ABS(UO))
      WRITE(*,*) REPEAT('~',72)
      WRITE(*,*) 'COURANT NUMBER: ',CFL
      WRITE(*,*) 'DIFFERENTIALS IN SPACE, DX, DY: ',DX, DY
      WRITE(*,*) 'TIME STEP: ',DT
      WRITE(*,*) REPEAT('~',72)
      NT = INT(FLOOR((TF-TO)/DT))
!
! ALLOCATE TIME VECTOR
!
      IF(DEBUG) WRITE(*,*) 'ALLOCATING AND COMPUTING TIME VECTOR'
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





