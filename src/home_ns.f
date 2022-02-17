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
      INTEGER :: TI,I,K
      INTEGER :: WTI = 0
      DOUBLE PRECISION :: CPUTINIT
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  READING USER INPUT FILE
!
      NAMELIST /NSCONF/ NU,RHO,NX,NY,XMIN,XMAX,YMIN,YMAX,CFL,TO,TF,
     &    PCD,WSSORD,TOLCG,MNITERD,SOLSING,PCS,KSSS,WSSORS,TOLSING,
     &    SING,MNITERM,SOLP,PCP,KSSP,WSSORP,TOLSOR,MNITERS,
     &    TIEP,WTIME,ISBOUND,BOUNDFILE,ISSTART,STARTFILE,DEBUG
!
      OPEN(1010, FILE = "nsconf.nml", STATUS = 'OLD')
      READ(1010, NML = NSCONF)
      CLOSE(1010)
      IF(DEBUG) WRITE(*,*) 'EXIT READING USER ENTRIES'
!
! COMPUTING INITIAL CPU-TIME
      CALL CPU_TIME(CPUTINIT)
!
!  WRITING HEADERS
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO WRITE_HEADERS'
      CALL WRITE_HEADERS
      IF(DEBUG) WRITE(*,*) 'EXIT WRITE_HEADERS'
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  ALLOCATING MEMORY INITIAL CONDITION, GRID, ETC.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  POINT_NS (GRID)
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO POINT_NS'
      CALL POINT_NS
      IF(DEBUG) WRITE(*,*) 'EXIT POINT_NS'
!
!  LECBOUND (INITIAL CONDITION, CFL AND TIME
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO LECBOUND'
      CALL LECBOUND
      IF(DEBUG) WRITE(*,*) 'EXIT LECBOUND'
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
      CALL POINT_REGULARIZATION(SING,TOLSING,MNITERM,W,LM,LPM,
     &    LQM,UL,VL,RM)
      IF(DEBUG) WRITE(*,*) 'EXIT POINT_REGULARIZATION'
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  TIME LOOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO TI = 1,SIZE(T)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF(DEBUG) WRITE(*,*) 'TIME LOOP INIT ON TIME = ',T(TI),':',TI
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF(TI.NE.1) THEN
!!!!!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! UPDATING BOUNDARIES TO NEXT TIME STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IF(DEBUG) WRITE(*,*) 'UPDATING PREVIOUS SOLUTION BOUNDARIES'
      IF(ISBOUND) THEN
      K = 1
      UO(NX*NY-NX+1) = BCONDU(TI,K)
      VO(NX*NY-NX+1) = BCONDV(TI,K)
      DO I = 1,SIZE(UPBOUND)
        K = K+1
        UO(UPBOUND(I)) = BCONDU(TI,K)
        VO(UPBOUND(I)) = BCONDV(TI,K)
      ENDDO
      UO(NX*NY) = BCONDU(TI,K)
      VO(NX*NY) = BCONDV(TI,K)
      ENDIF
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  FIRST FRACTIONAL STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  COMPUTING DERVIATIVES
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO FIRST FRACTIONAL STEP'
      CALL DIVER2D(DUDX,DUDY,DVDX,DVDY,UO,VO,BOUND,UPBOUNDI,
     &     DOBOUNDI,LEBOUNDI,RIBOUNDI)
      IF(DEBUG) WRITE(*,*) 'EXIT FIRST FRACTIONAL STEP'
!
!  COMPUTING UP AND VP, RESULT OF NON-LINEAL ADVECTION
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING UP AND VP'
      DO I = 1,NX*NY
        UP(I) = UO(I) - DT*(UO(I)*DUDX(I) + VO(I)*DUDY(I))
        VP(I) = VO(I) - DT*(UO(I)*DVDX(I) + VO(I)*DVDY(I))
      ENDDO
      IF(DEBUG) WRITE(*,*) 'END COMPUTING UP AND VP'
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  SECOND FRACTIONAL STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  COMPUTING DUPDX AND DVPDY
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING DUPDX AND DVPDY'
      CALL GRAD2D(DUPDX,DVPDY,UP,VP,UPBOUND,DOBOUND,LEBOUND,RIBOUND)
      IF(DEBUG) WRITE(*,*) 'END COMPUTING DUPDX AND DVPDY'
!
!  COMPUTING REGULARIZATION
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO POISSON'
      CALL POISSON(DUPDX,DVPDY,RM,LM,LPM,LQM,MNITERS,NITS,TOLSOR,
     &    TIEP,P,PP)
      IF(DEBUG) WRITE(*,*) 'EXIT FROM POISSON'
!
!  COMPUTING PRESURE GRADIENT
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING PRESURE GRADIENT'
      CALL GRAD2D(DPDX,DPDY,P,P,UPBOUND,DOBOUND,LEBOUND,RIBOUND)
      IF(DEBUG) WRITE(*,*) 'END COMPUTING PRESURE GRADIENT'
!
!  COMPUTING UPP AND VPP
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING UPP AND VPP'
      DO I=1,NX*NY
        UPP(I) = UP(I) - (DT/RHO)*DPDX(I)
        VPP(I) = VP(I) - (DT/RHO)*DPDY(I)
      ENDDO
      IF(DEBUG) WRITE(*,*) 'END COMPUTING UPP AND VPP'
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  THIRD FRACTIONAL STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  COMPUTING DIFFUSION
!
      IF(DEBUG) WRITE(*,*) 'GOING INTO DIFFUSION'
      CALL DIFFUSION(U,V,UO,VO,UPP,VPP,SX,SY,KM,MNITERD,NTIDX,NTIDY,
     &    TOLCG,BOUND,UPBOUND,DOBOUND,LEBOUND,RIBOUND,UPBOUNDI,DOBOUNDI,
     &    LEBOUNDI,RIBOUNDI)
      IF(DEBUG) WRITE(*,*) 'EXIT FROM DIFFUSION'
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! UPDATING REYNOLDS NUMBER
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      RE = MAX(MAXVAL(ABS(U)),MAXVAL(ABS(V)))*(MIN(XMAX-XMIN,YMAX-YMIN))
      RE = RE/NU
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!!!!
      ENDIF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  UPDATING VARIABLES AND WRITING OUTPUTS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IF(DEBUG) WRITE(*,*) 'CALLING UPDATE AND WRITING'
      CALL UPDATE_AND_WRITE(X,Y,U,V,PP,UO,VO,T,TI,NTIME,
     &     WTIME,NITS,NTIDX,NTIDY,TOLSOR,TOLCG,WTI,RE,CPUTINIT)
      IF(DEBUG) WRITE(*,*) 'EXIT UPDATE AND WRITING'

!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  END TIME LOOP
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  DEALLOCATING
!
      IF(DEBUG) WRITE(*,*) 'DEALLOCATING'
      CALL KILLEMALL
      IF(DEBUG) WRITE(*,*) 'HAPPY END :)"'
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      STOP 0
      END
