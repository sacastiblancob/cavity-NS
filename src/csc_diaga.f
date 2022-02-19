!                    ********************
                     SUBROUTINE CSC_DIAGA
!                    ********************
     & (MA,DA)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) SUBSTRACT DIAGONAL OF MATRIX A STORED IN CSC
!
!history  Sergio Castiblanco
!+        16/02/2022
!+        Translation for original Matlab implementation
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| MA       |-->| MATRIX A IN CSC                                      |
!| DA       |<->| VECTOR WITH DIAGONAL ELEMENTS OF MATRIX A            |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_DIAGA => CSC_DIAGA
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE(CSC_OBJ), INTENT(IN) :: MA
      DOUBLE PRECISION, DIMENSION(MA%NC), INTENT(INOUT) :: DA
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: J, I, M
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
! SEARCHING AND STORING DIAGONAL ELEMENTS
!
      M = MA%NC
      DO J=1,M
        DO I=MA%C(J),(MA%C(J+1)-MA%C(1))
          IF(MA%R(I).EQ.J)THEN
            DA(J)=MA%V(I)
          ENDIF
        ENDDO
      ENDDO
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_DIAGA
