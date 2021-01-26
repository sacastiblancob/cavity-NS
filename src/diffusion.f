!                    ********************
                     SUBROUTINE DIFFUSION
!                    ********************
     & (U,V,UPP,VPP,SX,SY,KM,MNITERD,NTIDX,NTIDY,TOLCG,BOUND,UPBOUND,
     &  DOBOUND,LEBOUND,RIBOUND,UPBOUNDI,DOBOUNDI,LEBOUNDI,RIBOUNDI)
!
!***********************************************************************
! 2D-NAVIER STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!BRIEF    1) COMPUTE DIFFUSION COEFFICIENTS AND MATRIX.
!
!HISTORY  SERGIO CASTIBLANCO
!+        16/01/2021
!+        TRANSLATION FOR ORIGINAL MATLAB IMPLEMENTATION
!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| U,V       |<->| U AND V SOLUTION OF THIRD FRACTIONAL STEP           |
!| UPP,VPP   |-->| U AND V SOLUTION OF SECOND FRACTIONAL STEP          |
!| SX,SY     |-->| DIFFUSIVITY DISCRETE COEFFICIENTS                   |
!| KM        |-->| DIFFUSION MATRIX                                    |
!| MNITERD   |-->| MAXIMUM NUMBER OF ITERATIONS FOR CG SOLVER          |
!| NTIDX,Y   |<->| NUMBER OF ITERATIONS TAKEN BY CG SOLVER             |
!| TOLCG     |-->| TOLERANCE FOR THE RESIDUAL OF CG SOLVER             |
!| BOUND     |-->| BOUNDARY INDICES                                    |
!| UPBOUND   |-->| TOP BOUNDARY INDICES                                |
!| DOBOUND   |-->| BOTTOM BOUNDARY INDICES                             |
!| RIBOUND   |-->| RIGHT BOUNDARY INDICES                              |
!| LEBOUND   |-->| LEFT BOUNDARY INDICES                               |
!| UPBOUNDI  |-->| TOP INTERNAL BOUNDARY INDICES                       |
!| DOBOUNDI  |-->| BOTTOM INTERNAL BOUNDARY INDICES                    |
!| RIBOUNDI  |-->| RIGHT INTERNAL BOUNDARY INDICES                     |
!| LEBOUNDI  |-->| LEFT INTERNAL BOUNDARY INDICES                      |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE
      USE DECLARATIONS_NUMERICAL, ONLY:NX,NY,DEBUG 
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(INOUT) :: U,V
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(IN) :: UPP,VPP
      DOUBLE PRECISION, INTENT(IN) :: SX,SY,TOLCG
      TYPE(CSC_OBJ), INTENT(IN) :: KM
      INTEGER, INTENT(IN) :: MNITERD
      INTEGER, INTENT(OUT) :: NTIDX, NTIDY
      INTEGER, DIMENSION(2*NX + 2*(NY-2)), INTENT(IN) :: BOUND
      INTEGER, DIMENSION(NX-2), INTENT(IN) :: UPBOUND, DOBOUND
      INTEGER, DIMENSION(NY-2), INTENT(IN) :: LEBOUND, RIBOUND
      INTEGER, DIMENSION(NX-4), INTENT(IN) :: UPBOUNDI, DOBOUNDI
      INTEGER, DIMENSION(NY-4), INTENT(IN) :: LEBOUNDI, RIBOUNDI
!
! IN SUBROUTINE VARIABLES
!
      INTEGER :: J,L
      DOUBLE PRECISION :: AP, AX, AY
      DOUBLE PRECISION, DIMENSION((NX-2)*(NY-2)) :: DRHSX, DRHSY
      DOUBLE PRECISION, DIMENSION((NX-2)*(NY-2)) :: UD,VD
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      AP = (1 + 2*SX + 2*SY)
      AX = -SX
      AY = -SY
!
!  COMPUTING RHS'S
!
      L = 1
      DO J = 1,NX*NY
        IF(ANY(J.EQ.BOUND)) THEN
            CONTINUE
        ELSEIF(ANY(J.EQ.UPBOUNDI)) THEN
            DRHSX(L) = UPP(J) - AY*1.D0
            DRHSY(L) = VPP(J)
            L=L+1;
        ELSEIF(ANY(J.EQ.RIBOUNDI)) THEN
            DRHSX(L) = UPP(J)
            DRHSY(L) = VPP(J)
            L=L+1;
        ELSEIF(ANY(J.EQ.DOBOUNDI)) THEN
            DRHSX(L) = UPP(J)
            DRHSY(L) = VPP(J)
            L=L+1;
        ELSEIF(ANY(J.EQ.LEBOUNDI)) THEN
            DRHSX(L) = UPP(J)
            DRHSY(L) = VPP(J)
            L=L+1;
        ELSEIF(J.EQ.(NX+2)) THEN
            DRHSX(L) = UPP(J)
            DRHSY(L) = VPP(J)
            L=L+1;
        ELSEIF(J.EQ.(2*NX-1)) THEN
            DRHSX(L) = UPP(J)
            DRHSY(L) = VPP(J)
            L=L+1;
        ELSEIF(J.EQ.(NX*NY-2*NX+2)) THEN
            DRHSX(L) = UPP(J) - AY*1.D0
            DRHSY(L) = VPP(J)
            L=L+1;
        ELSEIF(J.EQ.(NX*NY-NX-1)) THEN
            DRHSX(L) = UPP(J) - AY*1.D0
            DRHSY(L) = VPP(J)
            L=L+1;
        ELSE
            DRHSX(L) = UPP(J)
            DRHSY(L) = VPP(J)
            L=L+1;
        ENDIF
      ENDDO
!
!  SOLVING DIFFUSION EQUATION
!
      UD = 1.D0
      VD = 1.D0
      IF(DEBUG) WRITE(*,*) 'CALLING CG SOLVER FOR DIFFUSION'
      CALL CSC_CG(KM,DRHSX,UD,MNITERD,NTIDX,TOLCG,(NX-2)*(NY-2))
      CALL CSC_CG(KM,DRHSY,VD,MNITERD,NTIDY,TOLCG,(NX-2)*(NY-2))
      IF(DEBUG) WRITE(*,*) 'EXIT CG SOLVER FOR DIFFUSION'
      IF(DEBUG) WRITE(*,*) 'ITERATIONS, X, Y: ', NTIDX, NTIDY
!
!  COMPUTING U AND V (UD (OR VD) + BOUNDARY CONDITIONS)
!
      L = 1
      DO J = 1,NX*NY
        IF(ANY(J.EQ.UPBOUND)) THEN
            U(J) = 1.D0
            V(J) = 0.D0
        ELSEIF(ANY(J.EQ.RIBOUND)) THEN
            U(J) = 0.D0
            V(J) = 0.D0
        ELSEIF(ANY(J.EQ.DOBOUND)) THEN
            U(J) = 0.D0
            V(J) = 0.D0
        ELSEIF(ANY(J.EQ.LEBOUND)) THEN
            U(J) = 0.D0
            V(J) = 0.D0
        ELSEIF(J.EQ.1) THEN
            U(J) = 0.D0
            V(J) = 0.D0
        ELSEIF(J.EQ.(NX)) THEN
            U(J) = 0.D0
            V(J) = 0.D0
        ELSEIF(J.EQ.(NX*NY-NX+1)) THEN
            U(J) = 1.D0
            V(J) = 0.D0
        ELSEIF(J.EQ.(NX*NY)) THEN
            U(J) = 1.D0
            V(J) = 0.D0
        ELSE
            U(J) = UD(L)
            V(J) = VD(L)
            L=L+1;
        ENDIF
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE DIFFUSION




