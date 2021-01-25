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
      INTEGER :: AJ, BJ, P, AF, BF, IA, IB, J
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
! ALLOCATING MC
      CALL ALL_CSC(MC, MA%NR*MB%NR, MA%NC*MB%NC, MA%NZ*MB%NZ,NAMC)


!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE CSC_PRESOR
