!                     ******************
                      SUBROUTINE CSC_SUM
!                     ******************
     & (MA,MB,MC,NAMC)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) SUM MATRIX A AND B IN CSC STORAGE
!
!history  Sergio Castiblanco
!+        17/01/2021
!+        Translation for original Matlab implementation
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| MA       |-->| MATRIX A IN CSC                                      |
!| MB       |-->| MATRIX B IN CSC                                      |
!| MC       |<->| MATRIX C IN CSC (C = A + B)                          |
!| NAMC     |-->| NAME OF MATRIX C                                     |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_SUM => CSC_SUM
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE(CSC_OBJ), INTENT(IN) :: MA,MB
      TYPE(CSC_OBJ), INTENT(INOUT) :: MC
      CHARACTER(LEN=6), INTENT(IN) :: NAMC
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: ERR,I,J,K,IVA,IVB
      !ROWS OF A, B AND C PER COLUMN
      INTEGER, ALLOCATABLE :: RA(:), RB(:), RAB(:)
      !VALUES OF A AND B PER COLUMN
      DOUBLE PRECISION, ALLOCATABLE :: VA(:), VB(:)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      ERR = 0
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  CHECKING IF DIMENSION ARE THE SAME
!
      IF((MA%NC.NE.MB%NC).AND.(MA%NR.NE.MB%NR)) THEN
        WRITE(*,*) '!!!ERROR!!! DIMENSIONS FOR SUM ARE NOT THE SAME!!!'
        WRITE(*,*) 'BETWEEN ', MA%NAM,' AND ', MB%NAM
        STOP 1
      ENDIF
!
!  COMPUTING NEEDED VALUES FOR ALLOCATE AND ALLOCATING MC
!
      MC%NC = MA%NC
      MC%NR = MA%NR
      MC%NAM = NAMC
!
!     COLUMN RELATED INDICES VECTOR
      IF(.NOT.ASSOCIATED(MC%C)) THEN
        ALLOCATE(MC%C(MC%NC+1),STAT=ERR)
      ELSE
        WRITE(*,*) 'DEALLOCATING AND REALLOCATING %C OF OBJ ', MC%NAM
        DEALLOCATE(MC%C)
        ALLOCATE(MC%C(MC%NC+1),STAT=ERR)
      ENDIF
!
!     CHECKING ALLOCATE
      IF(ERR.NE.0) THEN
        WRITE(*,*) 'ALLOCATE ERROR IN ', MC%NAM
        STOP 1
      ENDIF
      ERR = 0
!
!     FILLING WITH ONES
      MC%C = 1
!
!     INITALIZING MC%NZ
      MC%NZ = 0
!
      DO J = 1,MA%NC
        !J = 1
        ALLOCATE(RA(SIZE( MA%R(MA%C(J):(MA%C(J+1)-MA%C(1))) )))
        RA = MA%R(MA%C(J):(MA%C(J+1)-MA%C(1)))
        ALLOCATE(RB(SIZE( MB%R(MB%C(J):(MB%C(J+1)-MB%C(1))) )))
        RB = MB%R(MB%C(J):(MB%C(J+1)-MB%C(1))) 
        CALL UNION(RA,RB,RAB)
        ! ! !WRITE(*,*) 'RA ',RA
        ! ! !WRITE(*,*) 'RB ',RB
        ! ! !WRITE(*,*) 'RAB ', RAB
        ALLOCATE(VA(SIZE( MA%R(MA%C(J):(MA%C(J+1)-MA%C(1))) )))
        VA = MA%V(MA%C(J):(MA%C(J+1)-MA%C(1)))
        ALLOCATE(VB(SIZE( MB%R(MB%C(J):(MB%C(J+1)-MB%C(1))) )))
        VB = MB%V(MB%C(J):(MB%C(J+1)-MB%C(1)))
        IVA = 1
        IVB = 1
        DO I = 1,SIZE(RAB)
          IF((ANY(RAB(I).EQ.RA)).AND.(ANY(RAB(I).EQ.RB))) THEN
            IF((VA(IVA)+VB(IVB)).EQ.0) THEN
              IVA = IVA + 1
              IVB = IVB + 1
            ELSE
              MC%NZ = MC%NZ + 1;
              MC%C(J+1:MC%NC+1) = MC%C(J+1:MC%NC+1) + 1
              IVA = IVA + 1
              IVB = IVB + 1
            ENDIF
          ELSE
            MC%NZ = MC%NZ + 1;
            MC%C(J+1:MC%NC+1) = MC%C(J+1:MC%NC+1) + 1
            IF(ANY(RAB(I).EQ.RA)) THEN
              IVA = IVA + 1
            ENDIF
            IF(ANY(RAB(I).EQ.RB)) THEN
              IVB = IVB + 1
            ENDIF
          ENDIF
          ! ! !WRITE(*,*) 'RAB(I)', RAB(I)
        ENDDO
        !
        ! DEALLOCATING DYNAMIC VARIABLES
        DEALLOCATE(RA)
        DEALLOCATE(RB)
        DEALLOCATE(RAB)
        DEALLOCATE(VA)
        DEALLOCATE(VB)
      ENDDO
!
!  ALLOCATING VECTOR OF ROWS AND VALUES OF MC
!
!
!     VALUES VECTOR
      IF(.NOT.ASSOCIATED(MC%V)) THEN
        ALLOCATE(MC%V(MC%NZ),STAT=ERR)
      ELSE
        WRITE(*,*) 'DEALLOCATING AND REALLOCATING %V OF OBJ ', MC%NAM
        DEALLOCATE(MC%V)
        ALLOCATE(MC%V(MC%NZ),STAT=ERR)
      ENDIF
!
!     CHECKING ALLOCATE
      IF(ERR.NE.0) THEN
        WRITE(*,*) 'ALLOCATE ERROR IN ', MC%NAM
        STOP 1
      ENDIF
      ERR = 0
!
!     FILLING WITH ZEROS
      MC%V = 0.0D0
!
!     ROW INDICES VECTOR
      IF(.NOT.ASSOCIATED(MC%R)) THEN
        ALLOCATE(MC%R(MC%NZ),STAT=ERR)
      ELSE
        WRITE(*,*) 'DEALLOCATING AND REALLOCATING %R OF OBJ ', MC%NAM
        DEALLOCATE(MC%R)
        ALLOCATE(MC%R(MC%NZ),STAT=ERR)
      ENDIF
!
!     CHECKING ALLOCATE
      IF(ERR.NE.0) THEN
        WRITE(*,*) 'ALLOCATE ERROR IN ', MC%NAM
        STOP 1
      ENDIF
      ERR = 0
!
!     FILLING WITH ZEROS
      MC%R = 0 
!
!  FILLING MC%R AND MC%B BY BASICALLY REPEATING LAST DO LOOP
!
      K = 1
      DO J = 1,MA%NC
        !J = 1
        ALLOCATE(RA(SIZE( MA%R(MA%C(J):(MA%C(J+1)-MA%C(1))) )))
        RA = MA%R(MA%C(J):(MA%C(J+1)-MA%C(1)))
        ALLOCATE(RB(SIZE( MB%R(MB%C(J):(MB%C(J+1)-MB%C(1))) )))
        RB = MB%R(MB%C(J):(MB%C(J+1)-MB%C(1))) 
        CALL UNION(RA,RB,RAB)
        ! ! !WRITE(*,*) 'RA ',RA
        ! ! !WRITE(*,*) 'RB ',RB
        ! ! !WRITE(*,*) 'RAB ', RAB
        ALLOCATE(VA(SIZE( MA%R(MA%C(J):(MA%C(J+1)-MA%C(1))) )))
        VA = MA%V(MA%C(J):(MA%C(J+1)-MA%C(1)))
        ALLOCATE(VB(SIZE( MB%R(MB%C(J):(MB%C(J+1)-MB%C(1))) )))
        VB = MB%V(MB%C(J):(MB%C(J+1)-MB%C(1)))
        IVA = 1
        IVB = 1
        DO I = 1,SIZE(RAB)
          IF((ANY(RAB(I).EQ.RA)).AND.(ANY(RAB(I).EQ.RB))) THEN
            IF((VA(IVA)+VB(IVB)).EQ.0) THEN
              IVA = IVA + 1
              IVB = IVB + 1
            ELSE
              MC%V(K) = VA(IVA) + VB(IVB);
              MC%R(K) = RAB(I)
              K = K + 1
              IVA = IVA + 1
              IVB = IVB + 1
            ENDIF
          ELSE
            IF(ANY(RAB(I).EQ.RA)) THEN
              MC%V(K) = VA(IVA)
              MC%R(K) = RAB(I)
              K = K + 1
              IVA = IVA + 1
            ENDIF
            IF(ANY(RAB(I).EQ.RB)) THEN
              MC%V(K) = VB(IVB)
              MC%R(K) = RAB(I)
              K = K + 1
              IVB = IVB + 1
            ENDIF
          ENDIF
        ENDDO
        !
        ! DEALLOCATING DYNAMIC VARIABLES
        DEALLOCATE(RA)
        DEALLOCATE(RB)
        DEALLOCATE(RAB)
        DEALLOCATE(VA)
        DEALLOCATE(VB)
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_SUM
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!                     ****************
                      SUBROUTINE UNION 
!                     ****************
     & (VA,VB,VC)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) COMPUTES THE UNION BETWEEN VECTORS A AND B
!         !!! WATNING: THE ELEMENTS IN A AND B MUST NOT BE REPITED IN
!               ITS OWN VECTOR
!
!history  Sergio Castiblanco
!+        17/01/2021
!+        Translation for original Matlab implementation
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| VA       |-->| VECTOR A                                             |
!| VB       |-->| VECTOR B                                             |
!| VC       |<->| VECTOR C (C = A U B)                                 |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_UNION => UNION
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      INTEGER, INTENT(IN) :: VA(:), VB(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: VC(:)
!
! IN SUBROUTINE VARIBALES
!
      INTEGER :: DA,DB,DC,I,J
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!  ALLOCATING VC
      DA = SIZE(VA)
      DB = SIZE(VB)
      DC = MAX(SIZE(VA),SIZE(VB))
      IF(SIZE(VA).EQ.(MIN(SIZE(VA),SIZE(VB)))) THEN
        DO I = 1,SIZE(VA)
          IF(.NOT.(ANY(VA(I).EQ.VB))) THEN
            DC = DC + 1
          ENDIF
        ENDDO
      ELSE
        DO I = 1,SIZE(VB)
          IF(.NOT.(ANY(VB(I).EQ.VA))) THEN
            DC = DC + 1
          ENDIF
        ENDDO
      ENDIF
      ALLOCATE(VC(DC))
!
!  FILLING VC
!
      IF(SIZE(VA).EQ.(MIN(SIZE(VA),SIZE(VB)))) THEN
        J = SIZE(VB)
        VC(1:J) = VB
        DO I = 1,SIZE(VA)
          IF(.NOT.(ANY(VA(I).EQ.VB))) THEN
            J = J + 1
            VC(J) = VA(I)
          ENDIF
        ENDDO
      ELSE
        J = SIZE(VA)
        VC(1:J) = VA
        DO I = 1,SIZE(VB)
          IF(.NOT.(ANY(VB(I).EQ.VA))) THEN
            J = J + 1
            VC(J) = VB(I)
          ENDIF
        ENDDO
      ENDIF
!
!  SORTING VC
!
      CALL QUICKSORT(VC,1,DC)
!      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE UNION




