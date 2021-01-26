!                 ***************************
                  SUBROUTINE LAPLACIAN_MATRIX
!                 ***************************
     & (DX,DY,LM)
!
!***********************************************************************
! 2D-NAVIER STOKES SOLVER - FINITE DIFFERENCES
!***********************************************************************
!
!brief    1) COMPUTE DIFFUSION COEFFICIENTS AND MATRIX.
!
!history  Sergio Castiblanco
!+        16/01/2021
!+        Translation for original Matlab implementation
!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| DX,DY     |-->| SPACE DIFFERENTIALS                                 |
!| LM        |<--| DIFFUSION MATRIX                                    |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE
      USE DECLARATIONS_NUMERICAL, ONLY:NX,NY,DEBUG 
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION, INTENT(IN)                  :: DX,DY
      TYPE(CSC_OBJ), INTENT(OUT)                    :: LM
!
! IN SUBROUTINE VARIABLES
!
      !MATRIX COEFFICIENTS
      DOUBLE PRECISION :: AP,AX,AY
!
! AUXILIAR MATRICES
!
      !MATRIX OF VALUES FOR L11
      DOUBLE PRECISION, DIMENSION(NY,1) :: HO
      !DIAGONAL INDICES FOR L11 AND L22 (=0, MAIN DIAGONAL)
      INTEGER, DIMENSION(1) :: DPO
      !MATRIX L11
      TYPE(CSC_OBJ) :: L11
      !MATRIX OF VALUES FOR L12
      DOUBLE PRECISION, DIMENSION(NX,3) ::  H1
      !DIAGONALS INDICES FOR L12
      INTEGER, DIMENSION(3) :: D1
      !MATRIX L12
      TYPE(CSC_OBJ) :: L12
      !MATRIX OF VALUES FOR L21
      DOUBLE PRECISION, DIMENSION(NY-1,2) :: H2
      !DIAGONAL INDICES FOR L21
      INTEGER, DIMENSION(2) :: D2
      !MATRIX L21
      TYPE(CSC_OBJ) :: L21
      !MATRIX OF VALUES FOR L22
      DOUBLE PRECISION, DIMENSION(NX,1) :: H3
      !MATRIX L22
      TYPE(CSC_OBJ) :: L22
      !MATRIX L1
      TYPE(CSC_OBJ) :: L1
      !MATRIX L2
      TYPE(CSC_OBJ) :: L2
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! COMPUTING LAPLACIAN MATRIX MATRIX
!
!
! COMPUTING MATRIX COEFFICIENTS
      IF(DY.GE.DX) THEN
        AP = -(2 + 2*((DX**2)/(DY**2)))
        AX = 1.D0
        AY = (DX**2)/(DY**2)
      ELSE
        AP = -(2*((DY**2)/(DX**2)) + 2)
        AX = (DY**2)/(DX**2)
        AY = 1.D0
      ENDIF

!
! COMPUTING MATRIX BY DIAGONALS AND KRONECKER PRODUCTS
      IF(DEBUG) WRITE(*,*) 'L11 MATRIX'
      !
      ! COMPUTING L11 MATRIX
      HO = 1.D0
      DPO = 0
      CALL CSC_DIAG(1,NY,HO,DPO,L11,'L11   ')
      ! ! !WRITE(*,*) 'L11V ', L11%V
      ! ! !WRITE(*,*) 'L11R ', L11%R
      ! ! !WRITE(*,*) 'L11C ', L11%C     
      !
      ! COMPUTING L12 MATRIX
      IF(DEBUG) WRITE(*,*) 'COMPUTING L12'
      ! H1
      H1 = 1.0D0
      H1(:,1) = AX*H1(:,1)
      H1(:,2) = AP*H1(:,2)
      H1(:,3) = AX*H1(:,3)
      H1(1,3) = 2.D0*AX
      H1(NX-1,1) = 2.D0*AX
      ! D1
      D1 = (/-1, 0, 1/)
      ! CSC_DIAG (H1,D1)
      CALL CSC_DIAG(3,NX,H1,D1,L12,'L12   ')
      ! ! !WRITE(*,*) 'L12V ', L12%V
      ! ! !WRITE(*,*) 'L12R ', L12%R
      ! ! !WRITE(*,*) 'L12C ', L12%C     
      !
      ! COMPUTING L21 MATRIX
      IF(DEBUG) WRITE(*,*) 'COMPUTING L21'
      H2 = 1.0D0
      H2(1,2) = 2.D0
      H2(NY-1,1) = 2.D0
      D2 = (/-1, 1/)
      ! CSC_DIAG (H2,D2)
      CALL CSC_DIAG(2,NY-1,H2,D2,L21,'L21   ')
      ! ! !WRITE(*,*) 'L21V ', L21%V
      ! ! !WRITE(*,*) 'L21R ', L21%R
      ! ! !WRITE(*,*) 'L21C ', L21%C     
      !
      ! COMPUTING L22 MATRIX
      IF(DEBUG) WRITE(*,*) 'COMPUTING L22'
      H3 = 1.0D0
      H3 = H3*AY
      ! CSC_DIAG(H3,0)
      CALL CSC_DIAG(1,NX,H3,DPO,L22,'L22   ')
      ! ! !WRITE(*,*) 'L22V ', L22%V
      ! ! !WRITE(*,*) 'L22R ', L22%R
      ! ! !WRITE(*,*) 'L22C ', L22%C     
!
!  COMPUTING KRONECKER PRODUCTS
!
      ! FOR MATRIX L1
      CALL CSC_KRON(L11,L12,L1,'L1    ')
      ! ! !WRITE(*,*) 'L1V ', L1%V
      ! ! !WRITE(*,*) 'L1R ', L1%R
      ! ! !WRITE(*,*) 'L1C ', L1%C     
      ! ! !WRITE(*,*) 'L1NC,L1NR,L1NZ ', L1%NC,L1%NR,L1%NZ
      ! FOR MATRIX L2
      CALL CSC_KRON(L21,L22,L2,'L2    ')
      ! ! !WRITE(*,*) 'L2V ', L2%V
      ! ! !WRITE(*,*) 'L2R ', L2%R
      ! ! !WRITE(*,*) 'L2C ', L2%C
      ! ! !WRITE(*,*) 'L2NC,L2NR,L2NZ ', L2%NC,L2%NR,L2%NZ
!
!  COMPUTING THE LAPLACIAN MATRIX AS THE SUM OF L1 AND L2
!
      CALL CSC_SUM(L1,L2,LM,'LM    ') 
      ! ! !WRITE(*,*) 'LV ', LM%V
      ! ! !WRITE(*,*) 'LR ', LM%R
      ! ! !WRITE(*,*) 'LC ', LM%C
      ! ! !WRITE(*,*) 'LNC,LNR,LNZ ', LM%NC,LM%NR,LM%NZ
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      END SUBROUTINE LAPLACIAN_MATRIX
 
