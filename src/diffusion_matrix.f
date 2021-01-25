!                 ***************************
                  SUBROUTINE DIFFUSION_MATRIX
!                 ***************************
     & (DX,DY,DT,VX,VY,KM,SX,SY)
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
!| DT        |-->| TIME DIFFERENTIAL                                   |
!| VX,VY     |-->| DIFFUSION COEFFICIENTS                              |
!| BOUND     |-->| BOUNDARY INDICES                                    |
!| UPBOUND   |-->| TOP BOUNDARY INDICES                                |
!| DOBOUND   |-->| BOTTOM BOUNDARY INDICES                             |
!| RIBOUND   |-->| RIGHT BOUNDARY INDICES                              |
!| LEBOUND   |-->| LEFT BOUNDARY INDICES                               |
!| KM        |<--| DIFFUSION MATRIX                                    |
!| SX,SY     |<->| DIFF COEFFICIENTS                                   |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE CSC_STORAGE
      USE DECLARATIONS_NUMERICAL, ONLY:NX,NY,DEBUG 
      IMPLICIT NONE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      DOUBLE PRECISION, INTENT(IN)                  :: DX,DY,DT,VX,VY
      DOUBLE PRECISION, INTENT(INOUT)               :: SX,SY
      TYPE(CSC_OBJ), INTENT(OUT)                    :: KM
!
! IN SUBROUTINE VARIABLES
!
      !MATRIX COEFFICIENTS
      DOUBLE PRECISION :: AP,AX,AY
!
! AUXILIAR MATRICES
!
      !MATRIX OF VALUES FOR K11
      DOUBLE PRECISION, DIMENSION(NY-2,1) :: HO
      !DIAGONAL INDICES FOR K11 AND K22 (=0, MAIN DIAGONAL)
      INTEGER, DIMENSION(1) :: DPO
      !MATRIX K11
      TYPE(CSC_OBJ) :: K11
      !MATRIX OF VALUES FOR K12
      DOUBLE PRECISION, DIMENSION(NX-2,3) ::  H1
      !DIAGONALS INDICES FOR K12
      INTEGER, DIMENSION(3) :: D1
      !MATRIX K12
      TYPE(CSC_OBJ) :: K12
      !MATRIX OF VALUES FOR K21
      DOUBLE PRECISION, DIMENSION(NY-3,2) :: H2
      !DIAGONAL INDICES FOR K21
      INTEGER, DIMENSION(2) :: D2
      !MATRIX K21
      TYPE(CSC_OBJ) :: K21
      !MATRIX OF VALUES FOR K22
      DOUBLE PRECISION, DIMENSION(NX-2,1) :: H3
      !MATRIX K22
      TYPE(CSC_OBJ) :: K22
      !MATRIX K1
      TYPE(CSC_OBJ) :: K1
      !MATRIX K2
      TYPE(CSC_OBJ) :: K2
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! SETTING DIFFUSSION STABILITY PARAMETERS
!
      IF(DEBUG) WRITE(*,*) 'COMPUTING SX AND SY'
      SX = ((DT*VX)/(DX**2))
      SY = ((DT*VY)/(DY**2))
      IF(DEBUG) WRITE(*,*) 'END COMPUTING SX AND SY'
!
! COMPUTING STTIFFNESS DIFFUSION MATRIX
!
!
! COMPUTING MATRIX COEFFICIENTS
      AP = (1 + 2*SX + 2*SY)
      AX = -SX
      AY = -SY
!
! COMPUTING MATRIX BY DIAGONALS AND KRONECKER PRODUCTS
      IF(DEBUG) WRITE(*,*) 'K11 MATRIX'
      !
      ! COMPUTING K11 MATRIX
      HO = 1D0
      DPO = 0
      CALL CSC_DIAG(1,NY-2,HO,DPO,K11,'K11   ')
      ! ! !WRITE(*,*) 'K11V ', K11%V
      ! ! !WRITE(*,*) 'K11R ', K11%R
      ! ! !WRITE(*,*) 'K11C ', K11%C     
      !
      ! COMPUTING K12 MATRIX
      IF(DEBUG) WRITE(*,*) 'COMPUTING K12'
      ! H1
      H1 = 1.0D0
      H1(:,1) = AX*H1(:,1)
      H1(:,2) = AP*H1(:,2)
      H1(:,3) = AX*H1(:,3)
      ! D1
      D1 = (/-1, 0, 1/)
      ! CSC_DIAG (H1,D1)
      CALL CSC_DIAG(3,NX-2,H1,D1,K12,'K12   ')
      ! ! !WRITE(*,*) 'K12V ', K12%V
      ! ! !WRITE(*,*) 'K12R ', K12%R
      ! ! !WRITE(*,*) 'K12C ', K12%C     
      !
      ! COMPUTING K21 MATRIX
      IF(DEBUG) WRITE(*,*) 'COMPUTING K21'
      H2 = 1.0D0
      D2 = (/-1, 1/)
      ! CSC_DIAG (H2,D2)
      CALL CSC_DIAG(2,NY-3,H2,D2,K21,'K21   ')
      ! ! !WRITE(*,*) 'K21V ', K21%V
      ! ! !WRITE(*,*) 'K21R ', K21%R
      ! ! !WRITE(*,*) 'K21C ', K21%C     
      !
      ! COMPUTING K22 MATRIX
      IF(DEBUG) WRITE(*,*) 'COMPUTING K22'
      H3 = 1.0D0
      H3 = H3*AY
      ! CSC_DIAG(H3,0)
      CALL CSC_DIAG(1,NX-2,H3,DPO,K22,'K22   ')
      ! ! !WRITE(*,*) 'K22V ', K22%V
      ! ! !WRITE(*,*) 'K22R ', K22%R
      ! ! !WRITE(*,*) 'K22C ', K22%C     
!
!  COMPUTING KRONECKER PRODUCTS
!
      ! FOR MATRIX K1
      CALL CSC_KRON(K11,K12,K1,'K1    ')
      ! ! !WRITE(*,*) 'K1V ', K1%V
      ! ! !WRITE(*,*) 'K1R ', K1%R
      ! ! !WRITE(*,*) 'K1C ', K1%C     
      ! ! !WRITE(*,*) 'K1NC,K1NR,K1NZ ', K1%NC,K1%NR,K1%NZ
      ! FOR MATRIX K2
      CALL CSC_KRON(K21,K22,K2,'K2    ')
      ! ! !WRITE(*,*) 'K2V ', K2%V
      ! ! !WRITE(*,*) 'K2R ', K2%R
      ! ! !WRITE(*,*) 'K2C ', K2%C
      ! ! !WRITE(*,*) 'K2NC,K2NR,K2NZ ', K2%NC,K2%NR,K2%NZ
!
!  COMPUTING THE DIFFUSION MATRIX AS THE SUM OF K1 AND K2
!
      CALL CSC_SUM(K1,K2,KM,'KM    ') 
      ! ! !WRITE(*,*) 'KV ', KM%V
      ! ! !WRITE(*,*) 'KR ', KM%R
      ! ! !WRITE(*,*) 'KC ', KM%C
      ! ! !WRITE(*,*) 'KNC,KNR,KNZ ', KM%NC,KM%NR,KM%NZ
! DISPLAYING COMPUTED SX AND SY
!
      WRITE(*,*) REPEAT('~',72)
      WRITE(*,*) 'SX VALUE: ', SX
      WRITE(*,*) 'SY VALUE: ', SY
      WRITE(*,*) REPEAT('~',72)
!
      END SUBROUTINE DIFFUSION_MATRIX
 
