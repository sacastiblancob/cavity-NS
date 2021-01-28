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
      USE DECLARATIONS_NUMERICAL, ONLY:DEBUG,NX,NY,DT,NT,BCOND
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
      INTEGER :: I,J,K,NLINES
      ! NUMBER OF ITERATIONS IN SOR SOLVER
      INTEGER :: NIU,NIV
      !AUXILIAR DOUBLE
      DOUBLE PRECISION :: R
      !DUMMY TO READ LINE IN INPUT FILE
      DOUBLE PRECISION, DIMENSION(23) :: BVAL
      !DUMMY MATRIX TO READ INPUTS
      DOUBLE PRECISION, DIMENSION(1000,23) :: BREAD
      !DUMMY MATRIX TO STORE ACTUAL INPUTS
      DOUBLE PRECISION, ALLOCATABLE :: BINP(:,:)
      !MATRIX WITH BOUNDARY CONDITIONS INTERPOLATED IN (X)
      DOUBLE PRECISION, ALLOCATABLE :: BCUX(:,:), BCVX(:,:)
      !MATRICES TO MAKE SPLINE INTERPOLATION
      TYPE(CSC_OBJ) :: MS
      !RIGHT HAND SIDES TO SOLVE INTERPOLATION SYSTEM
      DOUBLE PRECISION, DIMENSION(11) :: RHSU,RHSV
      !POLYNOMIAL COEFFICIENTS RESULT
      DOUBLE PRECISION, DIMENSION(11) :: AU,AV,CU,CV
      DOUBLE PRECISION, DIMENSION(10) :: BU,BV,DU,DV
      !VECTOR WITH TRANSFORMED X COORDINATES
      DOUBLE PRECISION, DIMENSION(NX) :: XT
      !DUMMY DOUBLE VARIABLES
      DOUBLE PRECISION :: DXC
      !TOLERANCE FOR SOR SOLVER
      DOUBLE PRECISION :: TOLU = 1E-14
      DOUBLE PRECISION :: TOLV = 1E-14
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
          ! ! !WRITE(*,*) 'BINP', BINP(I,:)
        ENDDO
!
! INTEPOLATION OVER GRID POINTS (FILLING BCUX AND BCVX)
!
        XT = X - MINVAL(X)
        XT = 10*(XT/MAXVAL(XT))
        ! ! !WRITE(*,*) 'XT ',XT
        ALLOCATE(BCUX(NLINES,NX))
        ALLOCATE(BCVX(NLINES,NX))
!
! FILLING MATRIX TO SOLVE SPLINE INTERPOLATION
!
        IF(DEBUG) WRITE(*,*) 'COMPUTING SPLINE INTERPOLATION'
        !DXC = (XMAX-XMIN)/10
        DXC = 1.D0
        CALL ALL_CSC(MS, 11, 11, (3*9)+4,'MS    ')
        MS%C(1)=1
        DO I = 1,11
          IF(I.EQ.1) THEN
            MS%C(I+1) = MS%C(I)+2
            MS%R(1:2) = (/1,2/)
            MS%V(1:2) = (/2.D0,2.D0/)
          ELSEIF(I.EQ.11) THEN
            MS%C(I+1) = MS%C(I)+2
            MS%V(30:31) = (/2.D0,2D0/)
            MS%R(30:31) = (/10,11/)
          ELSE
            MS%C(I+1) = MS%C(I)+3
            MS%R((I-1)*3:(I-1)*3+2) = (/I-1,I,I+1/)
            MS%V((I-1)*3:(I-1)*3+2) = (/1.D0,4.D0,1.D0/)
          ENDIF
        ENDDO
        ! ! !WRITE(*,*) 'MS%C ',MS%C
        ! ! !WRITE(*,*) 'MS%R ',MS%R
        ! ! !WRITE(*,*) 'MS%V ',MS%V
        DO I=1,NLINES
        !I=1
!!!!
        !
        ! COMPUTING A'S COEFFICIENTS
        !
        AU = BINP(I,2:12)
        AV = BINP(I,13:23)
        !
        ! COMPUTING C'S COEFFICIENTS
        !
        CU = 0.D0
        CU(1) = 1.D0
        CV = 0.D0
        CV(1) = 1.D0
        RHSU = 0.D0
        RHSV = 0.D0
        DO J = 1,11
          IF(J.EQ.1) THEN
            RHSU(J) = 3*(AU(2)-AU(1))
            RHSV(J) = 3*(AV(2)-AV(1))
          ELSEIF(J.EQ.11) THEN
            RHSU(J) = -3*(AU(11)-AU(10))
            RHSV(J) = -3*(AV(11)-AV(10))
          ELSE
            RHSU(J) = 3*(AU(J+1)-AU(J)) - 3*(AU(J)-AU(J-1))
            RHSV(J) = 3*(AV(J+1)-AV(J)) - 3*(AV(J)-AV(J-1))
          ENDIF
        ENDDO
        CALL CSC_CG(MS,RHSU,CU,1000,NIU,TOLU,11)
        CALL CSC_CG(MS,RHSV,CV,1000,NIV,TOLV,11)
        DO J=1,11
          IF(ABS(CU(J)).LT.TOLU) THEN
            CU(J) = 0.D0
          ENDIF
          IF(ABS(CV(J)).LT.TOLV) THEN
            CV(J) = 0.D0
          ENDIF
        ENDDO
        ! ! !WRITE(*,*) 'NIU ',NIU
        ! ! !WRITE(*,*) 'NIV ',NIV
        ! ! !WRITE(*,*) 'CU ',CU
        ! ! !WRITE(*,*) 'CV ',CV
        ! ! !WRITE(*,*) 'RHSU ', RHSU
        ! ! !WRITE(*,*) 'RHSV ', RHSV
        !
        ! COMPUTING B'S COEFFICIENTS
        !
        BU = 0.D0
        BV = 0.D0
        DO J = 1,10
          BU(J) = (AU(J+1)-AU(J)) - (1D0/3D0)*(2*CU(J)+CU(J+1))
          BV(J) = (AV(J+1)-AV(J)) - (1D0/3D0)*(2*CV(J)+CV(J+1))
        ENDDO
        !
        ! COMPUTING D'S COEFFICIENTS
        !
        DU = 0.D0
        DV = 0.D0
        DO J = 1,10
          DU(J) = (1D0/3D0)*(CU(J+1)-CU(J))
          DV(J) = (1D0/3D0)*(CV(J+1)-CV(J))
        ENDDO
        ! ! !WRITE(*,*) 'AU ',AU
        ! ! !WRITE(*,*) 'BU ',BU
        ! ! !WRITE(*,*) 'CU ',CU
        ! ! !WRITE(*,*) 'DU ',DU
        ! ! !WRITE(*,*) 'AV ',AV
        ! ! !WRITE(*,*) 'BV ',BV
        ! ! !WRITE(*,*) 'CV ',CV
        ! ! !WRITE(*,*) 'DV ',DV
        !
        ! INTERPOLATING ON X
        !
        R = 1.00001
        K = 1
        DO J = 1,NX
          IF(XT(J).LT.R) THEN
            BCUX(I,J) = AU(K) + BU(K)*(XT(J)-(K-1)) +
     &      CU(K)*(XT(J)-(K-1))**2 + DU(K)*(XT(J)-(K-1))**3
            BCVX(I,J) = AV(K) + BV(K)*(XT(J)-(K-1)) +
     &      CV(K)*(XT(J)-(K-1))**2 + DV(K)*(XT(J)-(K-1))**3
          ELSE
            K = K + 1
            R = R + 1.D0
            BCUX(I,J) = AU(K) + BU(K)*(XT(J)-(K-1)) +
     &      CU(K)*(XT(J)-(K-1))**2 + DU(K)*(XT(J)-(K-1))**3
            BCVX(I,J) = AV(K) + BV(K)*(XT(J)-(K-1)) +
     &      CV(K)*(XT(J)-(K-1))**2 + DV(K)*(XT(J)-(K-1))**3
          ENDIF
        ENDDO
        ! ! !WRITE(*,*) 'BCUX ', BCUX(I,:)
        ! ! !WRITE(*,*) 'BCVX ', BCVX(I,:)
!
!!!!
        ENDDO


      ENDIF
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE LECBOUND
