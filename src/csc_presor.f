!                    *********************
                     SUBROUTINE CSC_PRESOR
!                    *********************
     & (MA,MP,MQ,W,NAMP,NAMQ)
!
!***********************************************************************
! CSC PACKAGE - SERGIO CASTIBLANCO
!***********************************************************************
!
!brief    1) PREPARES THE MATRICES FOR SUCCESIVE OVER-RELAXATION SOLVER
!
!history  Sergio Castiblanco
!+        25/01/2021
!+        Translation for original Matlab implementation
!
!+
!!
!! P and Q:
!! Classic iterative methods use the next recurrence equation to solve 
!! the linear system of equations Ax=b:
!! 
!!               P * x(t+1) = Q * x(t) + b
!!
!!       For SOR:
!!            P = ((1/w)*D + L)
!!            Q = ((1/w - 1)*D - U); where: L = lower part of A
!!                                          U = upper part of A
!!                                          D = diagonal of A
!!
!! For the particular case when w = 1, we are in Gauss-Seidel method
!!
!!
!!      Sergio A. Castiblanco B. - Advanced Numerical Methods 
!!      Pontificia Universidad Javeriana - Bogota
!!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| MA       |-->| MATRIX A IN CSC                                      |
!| MP       |<->| MATRIX P IN CSC                                      |
!| MQ       |<->| MATRIX Q IN CSC                                      |
!| W        |<->| OVER-RELAXATION COEFFICIENT                          |
!| NAMP     |-->| NAME OF MATRIX P                                     |
!| NAMQ     |-->| NAME OF MATRIX Q                                     |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE, EX_CSC_PRESOR => CSC_PRESOR
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      TYPE(CSC_OBJ), INTENT(IN) :: MA
      TYPE(CSC_OBJ), INTENT(INOUT) :: MP,MQ
      DOUBLE PRECISION, INTENT(IN) :: W
      CHARACTER(LEN=6), INTENT(IN) :: NAMP,NAMQ
!
!  IN SUBROUTINE VARIABLES
!
      INTEGER :: J, I, QP, PP
      INTEGER :: NZP, NZQ
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! COMPUTING NZP AND NZQ
!
      NZP = 0
      NZQ = 0
!
      DO J = 1,MA%NC
        DO I = MA%C(J),MA%C(J+1)-1
          IF(MA%R(I).LT.J) THEN
            NZQ = NZQ + 1
          ELSEIF(MA%R(I).GT.J) THEN
            NZP = NZP + 1
          ELSEIF(MA%R(I).EQ.J) THEN
            IF ((W.GT.(1.D0-1E-8)).AND.(W.LT.1.D0+1E-8)) THEN
              NZP = NZP + 1
            ELSE
              NZP = NZP + 1
              NZQ = NZQ + 1
            ENDIF
          ENDIF
        ENDDO
      ENDDO
! 
! ALLOCATING MP AND MQ
!
      CALL ALL_CSC(MP, MA%NR, MA%NC, NZP, NAMP)
      CALL ALL_CSC(MQ, MA%NR, MA%NC, NZQ, NAMQ)
!
! FILLING MP AND MQ
!
      MQ%C(1) = 1
      MP%C(1) = 1
!
      QP = 1
      PP = 1
      IF((W.GT.(1.D0-1E-8)).AND.(W.LT.1.D0+1E-8)) THEN
        DO J = 1,MA%NC
          DO I = MA%C(J),MA%C(J+1)-1
            IF(MA%R(I).LT.J) THEN
              MQ%V(QP) = -MA%V(I)
              MQ%R(QP) = MA%R(I)
              QP = QP + 1
            ELSE
              MP%V(PP) = MA%V(I)
              MP%R(PP) = MA%R(I) 
              PP = PP + 1
            ENDIF
          ENDDO
          MQ%C(J+1) = QP
          MP%C(J+1) = PP
        ENDDO
      ELSE
        DO J = 1,MA%NC
          DO I = MA%C(J),MA%C(J+1)-1
            IF(MA%R(I).LT.J) THEN
              MQ%V(QP) = -MA%V(I)
              MQ%R(QP) = MA%R(I)
              QP = QP + 1
            ELSEIF(MA%R(I).GT.J) THEN
              MP%V(PP) = MA%V(I)
              MP%R(PP) = MA%R(I) 
              PP = PP + 1
            ELSE
              MQ%V(QP) = ((1D0/W) - 1D0)*MA%V(I)
              MQ%R(QP) = MA%R(I)
              QP = QP + 1
              MP%V(PP) = (1D0/W)*MA%V(I)
              MP%R(PP) = MA%R(I)
              PP = PP + 1
            ENDIF
          ENDDO
          MQ%C(J+1) = QP
          MP%C(J+1) = PP
        ENDDO
      ENDIF
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_PRESOR


