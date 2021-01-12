subroutine strider(bound,upbound,dobound,ribound,lebound,boundint,L)

USE gandalf
USE sam

integer, dimension((2*nx+2*ny)-4) :: bound
integer, dimension(nx) :: upbound, dobound
integer, dimension(ny-2) :: ribound, lebound
integer, dimension(2*(nx-2) + 2*(ny-4)) :: boundint
real, dimension(nx*ny,nx*ny) :: L
integer :: j

!Setting laplace operator matrix for presure treatment
j=1
do i=1,nx*ny
	L(i,j) = -30*(1/(12*dx**2) + 1/(12*dy**2))
	if (j<nx*ny) then	
		L(i,j+1) = (16/(12*dy**2))	
	end if
	if (j>1) then	
		L(i,j-1) = (16/(12*dy**2))
	end if	
	if (j>2) then
		L(i,j-2) = -(1/(12*dy**2))
	end if	
	if (j<(nx*ny-1)) then
		L(i,j+2) = -(1/(12*dy**2))
	end if	
	if (j>ny) then
		L(i,j-ny) = (16/(12*dx**2))
	end if
	if (j<(nx*ny-ny+1)) then
		L(i,j+ny) = (16/(12*dx**2))
	end if
	if (j>2*ny) then
		L(i,j-2*ny) = -(1/(12*dx**2))
	end if
	if (j<(nx*ny-2*ny+1)) then
		L(i,j+2*ny) = -(1/(12*dx**2))
	end if
	j=j+1
enddo

L(bound,:) = 0.0
L(boundint,:) = 0.0

do i=1,size(boundint)
	j = boundint(i)
	L(j,j) = (-2/(dx**2)) + (-2/(dy**2))
    L(j,j-ny) = 1/(dx**2)
    L(j,j+ny) = 1/(dx**2)
    L(j,j-1) = 1/(dy**2)
    L(j,j+1) = 1/(dy**2)
enddo

j=0
do i=1,size(upbound)
	j = upbound(i)
!    L(j,j) = 1/(dy);
!    L(j,j+1) = -1/(dy);
	L(j,j) = 3/(2*dy);
    L(j,j+1) = -4/(2*dy);
    L(j,j+2) = 1/(2*dy);
enddo
j=0
do i=1,size(dobound)
	j = dobound(i)
!    L(j,j) = -1/(dy)
!    L(j,j-1) = 1/(dy)
    L(j,j) = -3/(2*dy)
    L(j,j-1) = 4/(2*dy)
    L(j,j-2) = -1/(2*dy)
enddo
j=0
do i=1,size(lebound)
	j = lebound(i)
    L(j,j+2*ny) = -1/(2*dx)
    L(j,j+ny) = 4/(2*dx)
    L(j,j) = -3/(2*dx)
!    L(j,j) = -1/dx
!    L(j,j+ny) = 1/dx
!    L(j,j+1) = 1/(2*dy)
!    L(j,j-1) = -1/(2*dy)
enddo
j=0
do i=1,size(ribound)
	j = ribound(i)    
    L(j,j-2*ny) = 1/(2*dx)
    L(j,j-ny) = -4/(2*dx)
    L(j,j) = 3/(2*dx)
!    L(j,j) = 1/dx
!    L(j,j-ny) = -1/dx
!    L(j,j+1) = 1/(2*dy)
!    L(j,j-1) = -1/(2*dy)
enddo

endsubroutine strider


