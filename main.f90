program cavity_NS

! ################################################################################################################################ !
! Navier Stokes - Cavity case
! ################################################################################################################################ !
! 
! Made by Sergio Castiblanco
! Universidad Nacional de Colombia - Pontificia Universidad Javeriana
! 
! ################################################################################################################################ !

USE sam				!Physical variables
USE gandalf			!Time and space discretization variables

!Grid
real, dimension(1:nx) :: x
real, dimension(1:ny) :: y
!Boundary nodes (Internal and external{up,down,left,right})
integer, dimension((2*nx+2*ny)-4) :: bound
integer, dimension(nx) :: upbound, dobound
integer, dimension(ny-2) :: ribound, lebound
integer, dimension(2*(nx-2) + 2*(ny-4)) :: boundint
!Stiffness and Laplace matrices
real, dimension(nx*ny,nx*ny) :: K=0.0, L=0.0
!Velocities vectors
real, dimension(nx*ny) :: uo = 0.0
real, dimension(nx*ny) :: vo = 0.0
!Time discretization variables
real :: dt
real, dimension(:), allocatable :: t
integer :: nt
!Diffusion stabilty parameters
real :: Sx
real :: Sy
!Dummy variables
integer :: j
integer, dimension(nx*ny) :: dum = 0

! GRID GENERATION
call eru(x,y)

! BOUNDARY NODES IDENTIFICATION
call istari(bound,upbound,dobound,ribound,lebound,boundint)

! INITIAL CONDITION FOR VELOCITY
uo(1:nx) = 1.0

! SETTING TIME PARAMETERS
dt = min(dx,dy)*CFL/abs(uo(1))
nt = ((tf-to)/dt)

allocate(t(nt))
t(1) = to
do i=2,nt
	t(i) = t(i-1) + dt
enddo
write(*,*) dt, nt

!COMPUTING STIFNESS MATRIX (DIFUSSION STEP) AND ITS STABILITY TERMS
call peregrin(dt,bound,upbound,dobound,ribound,lebound,boundint,K,Sx,Sy)

!BIEN hasta aca, siguiente definir matriz de laplaciano de la presion
call strider(bound,upbound,dobound,ribound,lebound,boundint,L)

open(9,file='num.txt')
open(10,file='num2.txt')
do i=1,nx*ny
	write(9,*) K(i,:)
	write(10,*) L(i,:)
enddo 
close(10)
close(9)




deallocate(t)

endprogram cavity_NS
