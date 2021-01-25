!                     *******************
                      SUBROUTINE CSC_DIAG
!                     *******************
     & (NCH,NRH,H2,D,MA,NAMMA)
!
!***********************************************************************
! 2D-DIFFUSION SOLVER SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) CREATES A MATRIX WITH CSC STORAGE BY DIAGONALS
!
!history  Sergio Castiblanco
!+        16/01/2021
!+        Translation for original Matlab implementation
!
!+
!!
!!This function make a matrix in CSC storage by diagonals, the entries
!!  are:
!!      h2 --> matrix whose columns are gonna be the diagonals of A
!!      d ---> row-vector with the numbers of the diagonals in h2 cols
!!             0 for the main diagonal, >0 for lower part and >0 for
!!             upper part
!!             none of the diagonals are mandatory
!!      NRH,NCH ---> Number of rows and columns of h2 respectively
!!
!!      warning: if you have zeros in h2 you should remove them with 
!!          another function.
!!      NOTE: This function computes a square matrix with the given
!!      diagonals.
!!
!!      Example:
!!                    h2              d
!!               | 6  1  6 |      [-1 0 1]
!!               | 7  2  7 |
!!               | 8  3  8 |
!!               | 9  4  9 |
!!               | 10 5 10 |
!!      Returns:
!!                       A
!!              | 1  6  0  0  0 |
!!              | 6  2  7  0  0 |
!!              | 0  7  3  8  0 |
!!              | 0  0  8  4  9 |
!!              | 0  0  0  9  5 |
!!
!!      Note that 10's of h2(5,1) and h2(5,3) have dissapeared, those could
!!      be zero anyway
!!
!!      Of course A is stored in CSC, so the function returns
!!      Av, Ar, Ac (values, rows, columns)
!!
!!      Sergio A. Castiblanco B. - Avanced Numerical Methods
!!      Pontificia Universidad Javeriana - Bogota
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| NCH       |-->| NUMBER OF COLUMNS OF THE MATRIX WITH DIAGONALS      |
!| NRH       |-->| NUMBER OF ROWS OF THE MATRIX WITH DIAGONALS         |
!| H2        |-->| MATRIX WITH DIAGONALS                               |
!| D         |-->| VECTOR WITH DIAGONALS POSITIONS                     |
!| MA        |<->| RESULT MATRIX IN CSC                                |
!| NAMMA     |-->| NAME OF RESULT MATRIX IN CSC                        |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_DIAG => CSC_DIAG
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      INTEGER, INTENT(IN)          :: NCH
      INTEGER, INTENT(IN)          :: NRH
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NRH,NCH) :: H2
      INTEGER, INTENT(IN), DIMENSION(NCH)              :: D
      TYPE(CSC_OBJ), INTENT(INOUT)                     :: MA
      CHARACTER(LEN=6), INTENT(IN)                     :: NAMMA
!
!  IN SUBROUTINE VARIABLES
!
      !DEBUGGER INTERNAL OPTION
      LOGICAL :: DEBUG = .FALSE.
      !LOOP INTEGER
      INTEGER :: I,J,K,C,P,ROW,COL
      !PRINCIPAL DIAGONAL VALUE
      INTEGER :: DP
      !NUMBER OF COLUMNS OF MA
      INTEGER :: NCA
      !NUMBER OF NONZEROS OF MA
      INTEGER :: NZA
      !DUMMY VECTORS OF THE SAME SIZE OF D
      INTEGER, DIMENSION(NCH) :: CV
      !DUMMY LOGICAL INTEGER VECTOR OF THE SAME SIZE OF D
      INTEGER, DIMENSION(NCH) :: DLD
      !DUMMY INDEXING VECTOR
      INTEGER, ALLOCATABLE :: DLI(:), HROWS(:)
      !HCOLS
      DOUBLE PRECISION, ALLOCATABLE :: HCOLS(:,:)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!  DETERMINE PRINCIPAL DIAGONAL
      IF(DEBUG) WRITE(*,*) 'DETERMINING PRINCIPAL DIAGONAL'
      DP = MINVAL(ABS(D))
      DO I=1,NCH
        J = D(I)
        IF(ABS(J).EQ.DP) THEN
          DP=J
          EXIT
        ENDIF
      ENDDO
!
!  DETERMINE NUMBER OF NONZEROS OF RESULT MA MATRIX
      IF(DEBUG) WRITE(*,*) 'DETERMINING NUMBER OF NONZEROS OF A'
      NCA = NRH + ABS(DP)
      NZA = 0
      DO I=1,NCH
        NZA = NZA - ABS(D(I)) + NCA
      ENDDO
!
!  ALLOCATING MA
      IF(DEBUG) WRITE(*,*) 'GOING INTO ALL_CSC OF A'
      CALL ALL_CSC(MA, NCA, NCA, NZA, NAMMA)
      IF(DEBUG) WRITE(*,*) 'EXIT ALL_CSC OF A'
!
!  FILLING MATRIX A
      IF(DEBUG) WRITE(*,*) 'FILLING MATRIX A'
      !INIT OF MA%C(1) AND INTEGERS FOR LOOPING
      MA%C(1) = 1
      K = 0
      C = -NCA
      CV = D*0
      P = 1
      !LOOP OVER THE COLUMNS
      DO I=1,NCA
        !
        !COMPUTING DLD ((d<=k).*(d>c))
        CALL LOG2INT(NCH,(D.LE.K).AND.(D.GT.C),DLD)
        !
        !FILLING MA%C ON I+1
        MA%C(I+1) = MA%C(I) + SUM(DLD)
        !
        !COMPUTING Ar(Ac(i):Ac(i+1)-1) = sort(abs(d(logical((d<=k).*(d>c)))-1-k));
        CALL LOG2IND(NCH,((D.LE.K).AND.(D.GT.C)),DLI)
        DLI = ABS(D(DLI)-1-K)
        CALL QUICKSORT(DLI,1,SIZE(DLI))
        ! ! !WRITE(*,*) 'DLI', DLI
        MA%R(MA%C(I):MA%C(I+1)-1) = DLI
        !
        !COMPUTING CV = CV + ((d<=k).*(d>c))
        CV = CV + DLD
        !
        !COMPUTING HCOLS, WHICH IS DYNAMICALLY ALLOCATED
        CALL LOG2IND(NCH,((D.LE.K).AND.(D.GT.C)),DLI)
        ALLOCATE(HCOLS(NRH,SIZE(DLI)))
        HCOLS = H2(:,DLI)
        ! ! !WRITE(*,*) 'CV', CV
        ! ! !WRITE(*,*) 'H2(:,logical((d<=k).*(d>c)))', DLI, HCOLS
        ! ! !WRITE(*,*) 'cv(logical((d<=k).*(d>c)))', CV(DLI)
        !
        !COMPUTING HROWS
        CALL LOG2IND(NCH,((D.LE.K).AND.(D.GT.C)),HROWS)
        HROWS = CV(HROWS)
        HROWS = HROWS(SIZE(HROWS):1:-1)
        ! ! !WRITE(*,*) 'HROWS',HROWS
        !
        !UPDATING ROW AND COL
        ROW = 1
        COL = SIZE(HROWS)
        DO J=MA%C(I),(MA%C(I+1)-1)
          !
          !COMPUTING MA%V(P)
          MA%V(P) = HCOLS(HROWS(ROW),COL)
          !
          !UPDATING ROW, COL, P
          ROW = ROW + 1
          COL = COL - 1
          P = P + 1
        ENDDO
        !
        !UPDATING K, C
        K = K + 1
        C = C + 1
        DEALLOCATE(HCOLS)
      ENDDO
      ! ! !WRITE(*,*) 'AV ', MA%V
      ! ! !WRITE(*,*) 'AR ', MA%R
      ! ! !WRITE(*,*) 'AC ', MA%C     
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_DIAG
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!                     ******************
                      SUBROUTINE LOG2INT
!                     ******************
     & (DIMI,LOGI,INTE)
!
!***********************************************************************
! 2D-DIFFUSION SOLVER SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) FROM LOGICALS TO INTEGERS (.TRUE. > 1; .FALSE. > 0)
!
!history  Sergio Castiblanco
!+        17/01/2021
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_LOG2INT => LOG2INT
      IMPLICIT NONE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      INTEGER, INTENT(IN) :: DIMI
      INTEGER, DIMENSION(DIMI), INTENT(OUT) :: INTE
      LOGICAL, DIMENSION(DIMI), INTENT(IN) :: LOGI
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      INTEGER :: I
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DO I=1,DIMI
        IF(LOGI(I)) THEN
          INTE(I) = 1
        ELSE
          INTE(I) = 0
        ENDIF
      ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE LOG2INT
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!                     ******************
                      SUBROUTINE LOG2IND
!                     ******************
     & (DIMI,LOGI,INTE)
!
!***********************************************************************
! 2D-DIFFUSION SOLVER SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) FROM LOGICALS TO INTEGERS INDICES
!
!history  Sergio Castiblanco
!+        17/01/2021
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_LOG2IND => LOG2IND
      IMPLICIT NONE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      INTEGER, INTENT(IN) :: DIMI
      INTEGER, ALLOCATABLE, INTENT(OUT) :: INTE(:)
      LOGICAL, DIMENSION(DIMI), INTENT(IN) :: LOGI
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      INTEGER :: I,J,DIMI2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!  DEALLOCATING INTE IF IS ALLOCATED
      IF(ALLOCATED(INTE)) THEN
        DEALLOCATE(INTE)
      ENDIF
!
!  COMPUTING INDEXING
      DIMI2 = 0
      DO I=1,DIMI
        IF(LOGI(I)) THEN
          DIMI2 = DIMI2 + 1
        ENDIF
      ENDDO
      ! ALLOCATING INTE
      ALLOCATE(INTE(DIMI2))
      J = 1
      DO I=1,DIMI
        IF(LOGI(I)) THEN
          INTE(J) = I
          J = J+1
        ENDIF
      ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE LOG2IND
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
! changed to work with integers entries only (Sergio_Castiblanco)
!
      recursive subroutine quicksort(a, first, last)
      implicit none
      ! original ! double precision  x, t, a(*)
      ! modified
      double precision x
      integer t, a(*)
      !
      integer first, last
      integer i, j
!
      x = a( (first+last) / 2 )
      i = first
      j = last
      do
        do while (a(i) < x)
          i=i+1
        end do
        do while (x < a(j))
          j=j-1
        end do
        if (i >= j) exit
        t = a(i);  a(i) = a(j);  a(j) = t
        i=i+1
        j=j-1
      end do
      if (first < i-1) call quicksort(a, first, i-1)
      if (j+1 < last)  call quicksort(a, j+1, last)
      end subroutine quicksort


