!                        ******************
                         SUBROUTINE ALL_CSC
!                        ******************
!
     &(OBJ, N, M, NZ, NAM)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! CSC STORAGE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!brief     ALLOCATES A CSC STRUCTURE
!
!history   Sergio Castiblanco
!+         16/01/2021
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| OBJ      |<--| CSC_OBJECT                                           |
!| N        |-->| NUMBER OF ROWS                                       |
!| M        |-->| NUMBER OF COLUMNS                                    |
!| NZ       |-->| NUMBER OF NONZEROS                                   |
!| NAM      |-->| NAME OF OBJECT                                       |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_ALL_CSC => ALL_CSC
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! INTENT VARIABLES
!
      TYPE(CSC_OBJ), INTENT(INOUT) :: OBJ
      INTEGER, INTENT(IN)          :: N, M, NZ
      CHARACTER(LEN=6), INTENT(IN) :: NAM
!
! IN SUBROUTINE VARIABLES
!
      INTEGER :: ERR
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      ERR = 0
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!     NUMBER OF NONZEROS
      OBJ%NZ = NZ
!
!     NUMBER OF COLUMNS
      OBJ%NC = M
!
!     NUMBER OF ROWS
      OBJ%NR = N
!
!     NAME
      OBJ%NAM = NAM
!
!     VALUES VECTOR
      IF(.NOT.ASSOCIATED(OBJ%V)) THEN
        ALLOCATE(OBJ%V(OBJ%NZ),STAT=ERR)
      ELSE
        WRITE(*,*) 'DEALLOCATING AND REALLOCATING %V OF OBJ ', OBJ%NAM
        DEALLOCATE(OBJ%V)
        ALLOCATE(OBJ%V(OBJ%NZ),STAT=ERR)
      ENDIF
!
!     CHECKING ALLOCATE
      IF(ERR.NE.0) THEN
        WRITE(*,*) 'ALLOCATE ERROR IN ', OBJ%NAM
        STOP 1
      ENDIF
      ERR = 0
!
!     FILLING WITH ZEROS
      OBJ%V = 0.0D0
!
!     ROW INDICES VECTOR
      IF(.NOT.ASSOCIATED(OBJ%R)) THEN
        ALLOCATE(OBJ%R(OBJ%NZ),STAT=ERR)
      ELSE
        WRITE(*,*) 'DEALLOCATING AND REALLOCATING %R OF OBJ ', OBJ%NAM
        DEALLOCATE(OBJ%R)
        ALLOCATE(OBJ%R(OBJ%NZ),STAT=ERR)
      ENDIF
!
!     CHECKING ALLOCATE
      IF(ERR.NE.0) THEN
        WRITE(*,*) 'ALLOCATE ERROR IN ', OBJ%NAM
        STOP 1
      ENDIF
      ERR = 0
!
!     FILLING WITH ZEROS
      OBJ%R = 0 
!
!     COLUMN RELATED INDICES VECTOR
      IF(.NOT.ASSOCIATED(OBJ%C)) THEN
        ALLOCATE(OBJ%C(OBJ%NC+1),STAT=ERR)
      ELSE
        WRITE(*,*) 'DEALLOCATING AND REALLOCATING %C OF OBJ ', OBJ%NAM
        DEALLOCATE(OBJ%C)
        ALLOCATE(OBJ%C(OBJ%NC+1),STAT=ERR)
      ENDIF
!
!     CHECKING ALLOCATE
      IF(ERR.NE.0) THEN
        WRITE(*,*) 'ALLOCATE ERROR IN ', OBJ%NAM
        STOP 1
      ENDIF
!
!     FILLING WITH NC+1
      OBJ%C = 0
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      RETURN
      END








