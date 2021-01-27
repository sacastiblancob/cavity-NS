!                       *******************
                        SUBROUTINE LECBOUND
!                       *******************
     & (UO,VO,P,ISBOUND,BOUNDFILE,ISSTART,STARTFILE,T,TO,TF,XMAX,XMIN)
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
!
! IN SUBROUTINE VARIABLES
!
      INTEGER :: I,NLINES
      DOUBLE PRECISION, DIMENSION(23) :: BVAL
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  READING INITIAL CONDITION
!
      IF(ISBOUND) THEN
!
!  ALLOCATING BCOND, FIXED UP TO 1000,23
!
        ALLOCATE(BCOND(1000,23))
        BCOND = 0.D0
!
! READING BOUNDARY FILE
!
        OPEN (UNIT=1, FILE = TRIM(BOUNDFILE), STATUS='OLD')
        DO I = 1,1000
          IF((I.EQ.1).OR.(I.EQ.2)) THEN
            READ(1,*, END=10)
          ELSE
            READ(1,*, END=10) BVAL
            BCOND(I,:) = BVAL
            ! ! !WRITE(*,*) BCOND(I,:)
          ENDIF
        ENDDO
10      CLOSE(1)    
      ENDIF
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE LECBOUND
