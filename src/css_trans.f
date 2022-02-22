!                     ********************
                      SUBROUTINE CSS_TRANS
!                     ********************
     & (AV,AC,AR,NZ,NR,NC,TV,TC,TR)
!
!***********************************************************************
! CSS PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) COMPUTES THE TRANSPOSE OF MATRIX A IN CSS (CSR or CSC)
!
!history  Sergio Castiblanco
!+        21/02/2022
!+        Translation for original Matlab implementation
!
!+
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| AV       |-->| VALUES OF MATRIX A IN CSS                            |
!| AC       |-->| ROW INDICES IF CSC, COLUMN INDICES IF CSR            |
!| AR       |-->| COLUMN INDICES IF CSC, ROW INDICES IF CSR            |
!| NZ       |-->| NUMBER OF NON-ZEROS ENTRIES OF A                     |
!| NR       |-->| NUMBER OF ROWS OF A (CSR), NUMBER OF COLUMNS (CSC)   |
!| NC       |-->| NUMBER OF COLUMNS OF A(CSR), NUMBER OF ROWS(CSC)     |
!| TV       |<--| VALUES OF MATRIX A^T IN CSS                          |
!| TC       |<--| RESPECTIVE ROW OR COLUMN INDICES                     |
!| TR       |<--| RESPECTIVE COLUMN OR ROW INDICES                     |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSS_TRANS => CSS_TRANS
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      INTEGER, INTENT(IN) :: NZ,NR,NC
      DOUBLE PRECISION, DIMENSION(NZ), INTENT(IN) :: AV
      INTEGER, DIMENSION(NZ), INTENT(IN) :: AC
      INTEGER, DIMENSION(NR+1), INTENT(IN) :: AR
      DOUBLE PRECISION, DIMENSION(NZ), INTENT(OUT) :: TV
      INTEGER, DIMENSION(NZ), INTENT(OUT) :: TC
      INTEGER, DIMENSION(NC+1), INTENT(OUT) :: TR
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: MH,I,J,K,JP,IAA,IAB
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      MH = NC + 1
!
      DO I=1,NZ
        J = AC(I)+2
        IF(J.LE.MH)THEN
          TR(J) = TR(J)+1
        ENDIF
      ENDDO
!
      TR(1) = AR(1)
      TR(2) = AR(1)
!
      DO I=3,MH
        TR(I) = TR(I) + TR(I-1)
      ENDDO
!
      DO I=1,NR
        IAA = AR(I)
        IAB = AR(I+1)-AR(1)
        IF(IAB.LT.IAA)THEN
          EXIT
        ENDIF
        DO JP=IAA,IAB
          J = AC(JP) + 1
          K = TR(J)
          TC(K) = I
          TV(K) = AV(JP)
          TR(J) = K+1
        ENDDO
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSS_TRANS

