!                       *******************
                        SUBROUTINE LECBOUND
!                       *******************
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!***********************************************************************
! NAVIER-STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) READING BOUNDARY AND INITIAL CONDITION FILES,
!            ALSO INTERPOLATES, AND SET CFL AND COMPUTATIONAL TIME
!
!history  Sergio Castiblanco
!+        26/01/2021
!+        Translation for original Matlab implementation
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE CSC_STORAGE
      USE DECLARATIONS_PHYSIC
      USE DECLARATIONS_NUMERICAL
!
      IMPLICIT NONE
!    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! IN SUBROUTINE VARIABLES
!
      INTEGER :: I,J,K,INX,INY
      INTEGER :: NLINES = 1
      ! NUMBER OF ITERATIONS IN SOR SOLVER
      INTEGER :: NIU,NIV
      !AUXILIAR DOUBLE
      DOUBLE PRECISION :: R
      !DUMMY TO READ LINE IN INPUT FILE
      DOUBLE PRECISION, DIMENSION(23) :: BVAL
      !DUMMY MATRIX TO READ INPUTS
      DOUBLE PRECISION, DIMENSION(1000,23) :: BREAD
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
      !DUMMY CHARACTER
      CHARACTER(LEN=40) :: DCHA
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  READING BOUNDARY CONDITION FILE
!
      IF(ISBOUND) THEN
!
!  ALLOCATING BCOND, FIXED UP TO 1000,23
!
        BREAD = -999.D0
!
! READING BOUNDARY FILE
!
        OPEN (UNIT=10001, FILE = TRIM(BOUNDFILE), STATUS='OLD')
        K = 1
        DO I = 1,1000
          IF((I.EQ.1).OR.(I.EQ.2)) THEN
            READ(10001,*, END=10)
          ELSE
            READ(10001,*, END=10) BVAL
            BREAD(K,:) = BVAL
            K = K+1
            ! ! !WRITE(*,*) BCOND(I,:)
          ENDIF
        ENDDO
10      CLOSE(10001)    
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
        IF(NLINES.LT.2) THEN
          WRITE(*,*) 'ERROR IN BOUNDARY FILE, NOT ENOUGH INFORMATION, NO
     &T ENOUGH LINES :('
          WRITE(*,*) 'ERROR IN: ',TRIM(BOUNDFILE)
          STOP 1
        ENDIF
        ALLOCATE(BINP(NLINES,23))
        ALLOCATE(TINP(NLINES))
        DO I = 1,NLINES
          BINP(I,:) = BREAD(I,:)
          ! ! !WRITE(*,*) 'BINP', BINP(I,:)
        ENDDO
        TINP = BINP(:,1)
        !WRITE(*,*) 'TINP ', TINP, (TINP(1).EQ.TO)
        !
        ! PROOF ABOUT QUALITY OF INPUTS
        !
        IF(ABS(TINP(1) - TO).GT.0.0001) THEN
          WRITE(*,*) '!ERROR!, BOUNDARY FILE FIRST TIME ENTRY IS NOT THE
     & SAME AS INITIAL TIME GIVEN IN NSCONF.NML FILE :('
          WRITE(*,*) '!EXIT!, ERROR IN: ',TRIM(BOUNDFILE)
          STOP 1
        ENDIF
        ! ! !WRITE(*,*) 'TINP ', TINP(NLINES), TF, TINP(NLINES).LT.TF
        IF(TINP(NLINES).LE.TF) THEN
          WRITE(*,*) '!ERROR!, LAST TIME ENTRY IN BOUNDARY FILE IS LESS
     &OR EQUAL THAN FINAL TIME GIVEN IN NSCONF.NML FILE :('
          WRITE(*,*) 'PLEASE ADD RESPECTIVE ADDITIONAL LINES TO ENSURE T
     &HAT FINAL TIME IN BOUNDARY FILE IS GREATER THAN CONF FINAL TIME'
          WRITE(*,*) '!EXIT!, ERROR IN: ',TRIM(BOUNDFILE)
          STOP 1
        ENDIF
      ENDIF
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! READING INITIAL CONDITION FILE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF(ISSTART) THEN
        OPEN(UNIT=10002,FILE=TRIM(STARTFILE),STATUS='OLD')
        READ(10002,*)
        READ(10002,*)
        READ(10002,*) DCHA,DCHA,DCHA,INX,DCHA,INY
        ! ! !WRITE(*,*) 'DCHA ',INX,INY
        IF((INX.NE.NX).AND.(INY.NE.NY)) THEN
          CLOSE(UNIT=10002)
          WRITE(*,*) 'BAD INITIAL CONDITION, NUMBER OF NODES IN INITIAL 
     &CONDITION IS NOT THE SAME AS THE CONFIGURATION GIVEN ONES :('
          WRITE(*,*) 'CONF: ',NX,NY,'.NE. HOT START: ',INX,INY
          WRITE(*,*) 'ERROR IN: ',TRIM(STARTFILE)
          STOP 1
        ENDIF
        !
        ! FILLING INTIAL CONDITIONS
        !
        DO I = 1,NX*NY
          READ(10002,*) R,R,UO(I),VO(I),P(I),R
        ENDDO
        ! ! !WRITE(*,*) 'UO',UO
        ! ! !WRITE(*,*) 'VO',VO
        ! ! !WRITE(*,*) 'P',P
        CLOSE(UNIT=10002)
      ENDIF
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SETTLING CFL AND COMPUTATIONAL TIME
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! INITIAL CONDITION IF NO HOT START FILE
!
      IF(.NOT.ISSTART) THEN
        UO = 0D0
        VO = 0D0
        P = 1D0
        DO I=1,SIZE(UPBOUND)
          UO(UPBOUND(I)) = 1.0D0
        ENDDO
        UO(NX*NY) = 1.0D0
        UO(NX*NY-NX+1) = 1.0D0
      ENDIF
!
! REYNOLDS NUMBER
!
      RE = MAX(MAXVAL(ABS(UO)),MAXVAL(ABS(VO)))
      RE = RE*(MIN(XMAX-XMIN,YMAX-YMIN))
      RE = RE/NU
!
! TIME - TIME - TIME - TIME
!
! SETTING TIME PARAMETERS
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING DT AND NT'
      !
      ! CHANGING CFL IF IS GREATER THAN 1.2
      !
      IF(CFL.GT.1.2D0) THEN
        CFL = 1.2D0
      ENDIF
      !
      !COMPUTING DT WITH CFL, THE MINIMAL DIFFERENTIAL, AND MAXIMUM
      !INPUT VELOCITY
      !
      IF(ISBOUND) THEN
        DT = (MIN(DX,DY)*CFL)
        DT = DT/MAX(MAXVAL(ABS(UO)),MAXVAL(ABS(BINP(:,2:23))))
      ELSE
        DT = (MIN(DX,DY)*CFL)/MAXVAL(ABS(UO))
      ENDIF 
      !
      ! NUMBER OF TIME STEPS
      !
      NT = INT(FLOOR((TF-TO)/DT))
      !
      ! WRITING IN TERMINAL CFL, DX, DY, DT, NT
      !
      WRITE(*,*) REPEAT('~',72)
      WRITE(*,*) 'COURANT NUMBER: ',CFL
      WRITE(*,*) 'DIFFERENTIALS IN SPACE, DX, DY: ',DX, DY
      WRITE(*,*) 'TIME STEP: ',DT
      WRITE(*,*) 'TOTAL NUMBER OF TIME STEPS: ',NT
      WRITE(*,*) REPEAT('~',72)
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
!
! ALLOCATING BOUNDARY CONDITIONS MATRIX
!
      ALLOCATE(BCONDU(NT,NX))
      BCONDU = 0.D0
      ALLOCATE(BCONDV(NT,NX))
      BCONDV = 0.D0
!
! END TIME - END TIME - END TIME
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! INTERPOLATION IN SPACE FOR TOP BOUNDARY
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF(ISBOUND) THEN
!
! INTEPOLATION OVER GRID POINTS (FILLING BCUX AND BCVX)
!
        XT = X - MINVAL(X)
        XT = 10*(XT/MAXVAL(XT))
        ! ! !WRITE(*,*) 'XT ',XT
        ALLOCATE(BCUX(NLINES,NX))
        ALLOCATE(BCVX(NLINES,NX))
!
! FILLING MATRIX TO SOLVE C'S SPLINE INTERPOLATION
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
!
! INTERPOLATING IN DIMENSION X (FILLING BCUX AND BCVX WITH SPLINE)
!
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
        CALL CSC_CG(MS,RHSU,CU,2000,NIU,TOLU,11)
        CALL CSC_CG(MS,RHSV,CV,2000,NIV,TOLV,11)
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
        DO J=1,11
          ! CHECKING A ZEROS
          IF(ABS(AU(J)).LT.TOLU) THEN
            AU(J) = 0.D0
          ENDIF
          IF(ABS(AV(J)).LT.TOLV) THEN
            AV(J) = 0.D0
          ENDIF
          ! CHECKING C ZEROS
          IF(ABS(CU(J)).LT.TOLU) THEN
            CU(J) = 0.D0
          ENDIF
          IF(ABS(CV(J)).LT.TOLV) THEN
            CV(J) = 0.D0
          ENDIF
        ENDDO
        DO J=1,10
          ! CHECKING B ZEROS
          IF(ABS(BU(J)).LT.TOLU) THEN
            BU(J) = 0.D0
          ENDIF
          IF(ABS(BV(J)).LT.TOLV) THEN
            BV(J) = 0.D0
          ENDIF
          ! CHECKING D ZEROS
          IF(ABS(DU(J)).LT.TOLU) THEN
            DU(J) = 0.D0
          ENDIF
          IF(ABS(DV(J)).LT.TOLV) THEN
            DV(J) = 0.D0
          ENDIF
        ENDDO
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
!
! INTERPOLATING IN TIME WITH LINEAR INTERPOLATION
!
        K = 1
        R = TINP(K+1)
        DO I = 1,NT
        !I=1
!!!
        IF(T(I).LT.R) THEN
          DO J = 1,NX
            BCONDU(I,J) = BCUX(K,J) + (T(I)-TINP(K))*((BCUX(K+1,J) -
     &    BCUX(K,J))/(TINP(K+1) - TINP(K)))
            BCONDV(I,J) = BCVX(K,J) + (T(I)-TINP(K))*((BCVX(K+1,J) -
     &    BCVX(K,J))/(TINP(K+1) - TINP(K)))
          ENDDO
          ! ! !WRITE(*,*) 'ITIMU ',I,T(I),BCONDU(I,:)
          ! ! !WRITE(*,*) 'ITIMV ',I,T(I),BCONDV(I,:)
        ELSE
          K = K + 1
          R = TINP(K+1)
          DO J = 1,NX
            BCONDU(I,J) = BCUX(K,J) + (T(I)-TINP(K))*((BCUX(K+1,J) -
     &    BCUX(K,J))/(TINP(K+1) - TINP(K)))
            BCONDV(I,J) = BCVX(K,J) + (T(I)-TINP(K))*((BCVX(K+1,J) -
     &    BCVX(K,J))/(TINP(K+1) - TINP(K)))
          ENDDO
          ! ! !WRITE(*,*) 'ITIMU ',I,T(I),BCONDU(I,:)
          ! ! !WRITE(*,*) 'ITIMV ',I,T(I),BCONDV(I,:)
        ENDIF


!!!
        ENDDO
      ENDIF
!
! MODIFICATION OVER INITIAL CONDITION IF ISBOUND BUT NOT ISSTART
!
      IF(ISBOUND.AND.(.NOT.ISSTART)) THEN
        K = 1
        UO(NX*NY-NX+1) = BCONDU(1,K)
        VO(NX*NY-NX+1) = BCONDV(1,K)
        DO I = 1,SIZE(UPBOUND)
          K = K+1
          UO(UPBOUND(I)) = BCONDU(1,K)
          VO(UPBOUND(I)) = BCONDV(1,K)
        ENDDO
        UO(NX*NY) = BCONDU(1,K)
        VO(NX*NY) = BCONDV(1,K)
      ENDIF

!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE LECBOUND
