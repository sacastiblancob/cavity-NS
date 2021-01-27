!                       *******************
                        SUBROUTINE LECBOUND
!                       *******************
     & (UO,VO,P,ISBOUND,BOUNDFILE,ISSTART,STARTFILE,T,TO,TF,XMAX,XMIN,X)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!***********************************************************************
! NAVIER-STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) READING BOUNDARY AND INITIAL CONDITION FILES,
!            ALSO INTERPOLATES
!
!history  Sergio Castiblanco
!+        26/01/2021
!+        Translation for original Matlab implementation
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| UO,VO,P   |<--| VECTORS WITH INITIAL CONDITIONS                     |
!| ISBOUND   |-->| LOGICAL IF USER BOUNDARY CONDITIONS                 |
!| BOUNDFILE |-->| NAME OF BOUNDARY FILE                               |
!| ISSTART   |-->| LOGICAL IF HOT START                                |
!| STARTFILE |-->| NAME OF HOT START FILE                              |
!| T         |-->| VECTOR WITH TIMES                                   |
!| TO        |-->| INITIAL TIME                                        |
!| TF        |-->| FINAL TIME                                          |
!| XMAX,XMIN |-->| MAX AND MIN X COORDINATES                           |
!| X         |-->| VECTOR WITH X COORDINATES                           |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE CSC_STORAGE
      USE DECLARATIONS_NUMERICAL, ONLY:NX,NY,DT,NT,BCOND
!
      IMPLICIT NONE
!    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! INTENT VARIABLES
!
      DOUBLE PRECISION, DIMENSION(NX*NY), INTENT(INOUT) :: UO,VO,P
      LOGICAL, INTENT(IN) :: ISBOUND, ISSTART
      CHARACTER(LEN=250), INTENT(INOUT) :: BOUNDFILE, STARTFILE
      DOUBLE PRECISION, DIMENSION(NT), INTENT(IN) :: T
      DOUBLE PRECISION, INTENT(IN) :: TO, TF, XMAX, XMIN
      DOUBLE PRECISION, DIMENSION(NX), INTENT(IN) :: X
!
! IN SUBROUTINE VARIABLES
!
      INTEGER :: I,K,NLINES
      DOUBLE PRECISION :: R
      !DUMMY TO READ LINE IN INPUT FILE
      DOUBLE PRECISION, DIMENSION(23) :: BVAL
      !DUMMY MATRIX TO READ INPUTS
      DOUBLE PRECISION, DIMENSION(1000,23) :: BREAD
      !DUMMY MATRIX TO STORE ACTUAL INPUTS
      DOUBLE PRECISION, ALLOCATABLE :: BINP(:,:)
      !MATRIX WITH BOUNDARY CONDITIONS INTERPOLATED IN (X)
      DOUBLE PRECISION, ALLOCATABLE :: BCUX(:,:), BCVX(:,:)
      !MATRIX TO MAKE SPLINE INTERPOLATION
      TYPE(CSC_OBJ) :: SIM
      !VECTOR WITH X COORDINATES OF THE 11 POINTS IN INPUT FILE
      DOUBLE PRECISION, DIMENSION(11) :: XCOORD
      !DUMMY DOUBLE VARIABLES
      DOUBLE PRECISION :: DDUM, DXC
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  READING INITIAL CONDITION
!
      IF(ISBOUND) THEN
!
!  ALLOCATING BCOND, FIXED UP TO 1000,23
!
        BREAD = -999.D0
!
! READING BOUNDARY FILE
!
        OPEN (UNIT=1, FILE = TRIM(BOUNDFILE), STATUS='OLD')
        K = 1
        DO I = 1,1000
          IF((I.EQ.1).OR.(I.EQ.2)) THEN
            READ(1,*, END=10)
          ELSE
            READ(1,*, END=10) BVAL
            BREAD(K,:) = BVAL
            K = K+1
            ! ! !WRITE(*,*) BCOND(I,:)
          ENDIF
        ENDDO
10      CLOSE(1)    
!
! TAKING VALUES FROM BREAD TO BINP
!
        NLINES = 0
        I = 1
        R = 0
        DO WHILE((R.LT.-999.00001).OR.(R.GT.-998.9999))
          R = BREAD(I,1)
          NLINES = NLINES + 1
          I = I+1
        ENDDO
        NLINES = NLINES-1
        ! ! !WRITE(*,*) 'NLINES...', NLINES
        ALLOCATE(BINP(NLINES,23))
        DO I = 1,NLINES
          BINP(I,:) = BREAD(I,:)
          WRITE(*,*) 'BINP', BINP(I,:)
        ENDDO
!
! INTEPOLATION OVER GRID POINTS (FILLING BCUX AND BCVX)
!
        ALLOCATE(BCUX(NLINES,NX))
        ALLOCATE(BCVX(NLINES,NX))
        CALL ALL_CSC(SIM, 40, 40, (18*9)+14,'SIM   ')
!
! FILLING MATRIX TO SOLVE SPLINE INTERPOLATION
!
        XCOORD = 0.D0
        DXC = (XMAX-XMIN)/10
        IF((XMIN.GT.(-0.001)).AND.(XMIN.LT.0.001)) THEN
          XCOORD(1) = XMIN + 0.001
        ELSE
          XCOORD(1) = XMIN
        ENDIF
        IF((XMAX.GT.(-0.001)).AND.(XMAX.LT.0.001)) THEN
          XCOORD(11) = XMAX - 0.001
        ELSE
          XCOORD(11) = XMAX
        ENDIF
        DO I=2,10
          DDUM = XCOORD(I-1) + DXC
          IF((DDUM.GT.(-0.001)).AND.(DDUM.LT.0.001)) THEN
            XCOORD(I) = DDUM + 0.001
          ELSE
            XCOORD(I) = DDUM
          ENDIF
        ENDDO
        WRITE(*,*) 'XCOORD', XCOORD
        !SIM%C(1) = 1
        !SIM%C(
        !DO J = 1,SIM%NC
        !  DO I = SIM

      ENDIF

!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE LECBOUND
