subroutine istari(bound,upbound,dobound,ribound,lebound,boundint)

! ################################################################################################################################ !
! BOUNDARY POSITIONS ARRAYS
! ################################################################################################################################ !
!As the blue istari, that went to the deep boundaries of the Middle Earth

USE gandalf

!integer, parameter :: nx,ny
integer, dimension((2*nx+2*ny)-4) :: bound
integer, dimension(nx) :: upbound, dobound
integer, dimension(ny-2) :: ribound, lebound
integer, dimension(2*(nx-2) + 2*(ny-4)) :: boundint

pos = (/(i,i=1,nx*ny)/)

upbound = pos(1:nx)
dobound = pos(nx*ny-nx+1:nx*ny)
ribound(1) = 2*nx
do i=2,(ny-2)
	ribound(i) = ribound(i-1) + nx
enddo
!ribound = pos(2*nx:nx:nx*ny-nx)
lebound(1) = nx+1
do i=2,(ny-2)
	lebound(i) = lebound(i-1) + nx
enddo
!lebound = pos(nx+1:nx:nx*ny-2*nx+1)
bound = [upbound,dobound,ribound,lebound]
boundint = [pos(nx+2:2*nx-1),pos(nx*ny-2*nx+2:nx*ny-nx-1),ribound(2:ny-3)-1,lebound(2:ny-3)+1]

!write(*,*) pos
!write(*,*) bound
!write(*,*) boundint
!write(*,*) posmod 

endsubroutine istari
