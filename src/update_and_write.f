!                 ***************************
                  SUBROUTINE UPDATE_AND_WRITE
!                 ***************************
     & (X,Y,U,V,PP,UO,VO,T,TI,NTIME,WTIME,NITS,NTIDX,NTIDY,
     &    TOLSOR,TOLCG,WTI,RE,CPUTINIT)
!
!***********************************************************************
! 2D-NAVIER-STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) UPDATE VARIABLES, AND WRITE IF WTIME.
!
!history  Sergio Castiblanco
!+        26/01/2021
!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| X,Y       |-->| COORDINATES                                         |
!| U,V       |-->| SOLUTION OF VELOCITY                                |
!| PP        |-->| PRESURE FOR PLOTTING                                |
!| UO,VO     |<->| SOLUTION TO BE UPDATED                              |
!| BOUND     |-->| BOUNDARY INDICES                                    |
!| T         |-->| VECTOR WITH TIMES                                   |
!| TI        |-->| TIME STEP AT WE ARE                                 |
!| NTIME     |<->| NEXT WRITING TIME                                   |
!| WTIME     |-->| TIME INTERVAL FOR WRITING                           |
!| NITS      |-->| NUMBER OF ITERATIONS OF SOR SOLVER (PRESURE)        |
!| NTIDX     |-->| NUMBER OF ITERATIONS OF CG SOLVER (DIFFUSION OF U)  |
!| NTIDY     |-->| NUMBER OF ITERATIONS OF CG SOLVER (DIFFUSION OF V)  |
!| TOLSOR    |-->| GIVEN SOR SOLVER TOLERANCE                          |
!| TOLCG     |-->| GIVEN CG SOLVER TOLERANCE                           |
!| WTI       |-->| HANDLING WRITE TIME                                 |
!| RE        |-->| REYNOLDS NUMBER                                     |
!| CPUTINIT  |-->| CPU INITAL TIME                                     |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE DECLARATIONS_NUMERICAL, ONLY:NX,NY,NT,DT
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NX) :: X
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NY) :: Y
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NX*NY) :: U,V,PP 
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NX*NY) :: UO, VO
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NT) :: T
      INTEGER, INTENT(IN) :: TI
      DOUBLE PRECISION, INTENT(INOUT) :: NTIME
      DOUBLE PRECISION, INTENT(IN) :: WTIME, CPUTINIT
      INTEGER, INTENT(IN) :: NITS,NTIDX,NTIDY
      DOUBLE PRECISION, INTENT(IN) :: TOLSOR,TOLCG,RE
      INTEGER, INTENT(INOUT) :: WTI
!
! IN SUBROUTINE VARIABLES
!
      INTEGER :: I,J,K
      DOUBLE PRECISION :: C1,C2,C3,C4,C5,C6,CPUT
      CHARACTER(LEN=80) :: FOUT
      CHARACTER(LEN=9)  :: FMT2
      CHARACTER(LEN=34) :: FSPEC
      DOUBLE PRECISION :: EPSI = 1E-6
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  WRITING INITIAL CONDITION IF TI == 1
!
      IF(TI.EQ.1) THEN
        WTI = 10
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! WRITING IN TERMINAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        WRITE(*,*) REPEAT('~',72)
        WRITE(*,*) 'INIT DONE'
        WRITE(*,*) 'WRITING INITIAL CONDITION'
        WRITE(*,*) 'REYNOLDS NUMBER: ',RE
        WRITE(*,*) 'TIME STEP VALUE: ',DT
        WRITE(*,*) 'INITIAL TIME: ',T(1)
        WRITE(*,*) 'FINAL TIME: ',T(SIZE(T))
        WRITE(*,*) REPEAT('~',72)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        FMT2 = '(i4.4,a4)'
        FSPEC = './res/2D_Navier-Stokes_'
        FOUT = FSPEC
!
        WRITE(UNIT=FOUT(24:34),FMT=FMT2) WTI,'.dat'
!
        OPEN(WTI,FILE=FOUT)
!
!  HEADER
        WRITE(WTI,*) 'TITLE = "NAVIER STOKES 2D - FD"'
        WRITE(WTI,*) 'VARIABLES = "X" "Y" "U" "V" "P" "U_MAG"'
        WRITE(WTI,*) ' ZONE F=POINT, I=',NX,', J= ',NY
!
!  WRITING RESULTS
        K = 1
        DO J = 1,NY
          DO I = 1,NX
            C1 = X(I)
            C2 = Y(J)
            C3 = UO(K)
            C4 = VO(K)
            C5 = PP(K)
            C6 = SQRT(C3*C3 + C4*C4)
            WRITE(WTI,'(*(F18.8))') C1, C2, C3, C4, C5, C6
            K = K + 1
          ENDDO
        ENDDO
        CLOSE(WTI)
        WTI = WTI + 1
        NTIME = T(1) + WTIME
      ELSE !UPDATE VARIABLES AND WRITING IN REGULAR SIMULATION
!
!  UPDATING UO,VO
!
        DO I = 1,SIZE(UO)
          UO(I) = U(I)
          VO(I) = V(I)
        ENDDO
!
!  WRITING IF IS TIME TO WRITE
!
        IF((T(TI).GT.(NTIME-EPSI)).OR.(T(TI).EQ.T(NT))) THEN
!
! COMPUTING CPU TIME
          CALL CPU_TIME(CPUT)        
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! WRITING IN TERMINAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          WRITE(*,*) REPEAT('~',72)
          WRITE(*,*) 'TIME: ',T(TI),'(s)'
          WRITE(*,*) 'SIMULATION TIME: ',CPUT-CPUTINIT,'(s)'
          WRITE(*,*) 'TIME STEP: ',TI,' OF ',SIZE(T)
          WRITE(*,*) 'TIME STEP VALUE: ',DT
          WRITE(*,*) 'REYNOLDS NUMBER: ',RE
          WRITE(*,*) 'LAST POISSON SOLVER TOL FOR THE RESIDUAL: ' 
          WRITE(*,*) TOLSOR
          WRITE(*,*) 'LAST POISSON SOLVER NUMBER OF ITERATIONS: '
          WRITE(*,*) NITS
          WRITE(*,*) 'LAST DIFFUSION SOLVER TOL FOR THE RESIDUAL: ' 
          WRITE(*,*) TOLCG
          WRITE(*,*) 'LAST DIFFUSION SOLVER NUMBER OF ITERATIONS: '
          WRITE(*,*) 'U: ',NTIDX,', V: ',NTIDY
          WRITE(*,*) 'WRITING FILE NUMBER: ', (WTI-10)
          WRITE(*,*) REPEAT('~',72)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  WRITING FILES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          FMT2 = '(i4.4,a4)'
          FSPEC = './res/2D_Navier-Stokes_'
          FOUT = FSPEC
!
          WRITE(UNIT=FOUT(24:34),FMT=FMT2) WTI,'.dat'
!
! HERE IS CALLING STOP 1 WITHOUT REASON
          OPEN(WTI,FILE=FOUT)
!
!  HEADER
          WRITE(WTI,*) 'TITLE = "NAVIER STOKES 2D - FD"'
          WRITE(WTI,*) 'VARIABLES = "X" "Y" "U" "V" "P" "U_MAG"'
          WRITE(WTI,*) ' ZONE F=POINT, I=',NX,', J= ',NY
!
!  WRITING RESULTS
          K = 1
          DO J = 1,NY
            DO I = 1,NX
              C1 = X(I)
              C2 = Y(J)
              C3 = UO(K)
              C4 = VO(K)
              C5 = PP(K)
              C6 = SQRT(C3*C3 + C4*C4)
              WRITE(WTI,'(*(F18.8))') C1, C2, C3, C4, C5, C6
              K = K + 1
            ENDDO
          ENDDO
          CLOSE(WTI)
          WTI = WTI + 1
          NTIME = NTIME + WTIME
        ENDIF
      ENDIF
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE UPDATE_AND_WRITE

