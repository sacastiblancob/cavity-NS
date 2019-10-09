subroutine peregrin(dt,bound,upbound,dobound,ribound,lebound,boundint,K,Sx,Sy)

USE gandalf
USE sam

real :: Sx
real :: Sy
real :: dt
integer, dimension((2*nx+2*ny)-4) :: bound
integer, dimension(nx) :: upbound, dobound
integer, dimension(ny-2) :: ribound, lebound
integer, dimension(2*(nx-2) + 2*(ny-4)) :: boundint
real, dimension(nx*ny,nx*ny) :: K
integer :: j

!Setting difussion stability parameters 
Sx = ((dt*nu)/(dx**2))		!X direction
Sy = ((dt*nu)/(dy**2))		!Y direction

!Setting sttifness difussion matrix
j=1
do i=1,nx*ny
	K(i,j) = -(1 + 30*Sx/12 + 30*Sy/12)
	if (j<nx*ny) then	
		K(i,j+1) = (16*Sy/12)	
	end if
	if (j>1) then	
		K(i,j-1) = (16*Sy/12)
	end if	
	if (j>2) then
		K(i,j-2) = -(Sy/12)
	end if	
	if (j<(nx*ny-1)) then
		K(i,j+2) = -(Sy/12)
	end if	
	if (j>ny) then
		K(i,j-ny) = (16*Sx/12)
	end if
	if (j<(nx*ny-ny+1)) then
		K(i,j+ny) = (16*Sx/12)
	end if
	if (j>2*ny) then
		K(i,j-2*ny) = -(Sx/12)
	end if
	if (j<(nx*ny-2*ny+1)) then
		K(i,j+2*ny) = -(Sx/12)
	end if
	j=j+1
enddo

K(bound,:) = 0.0
K(boundint,:) = 0.0

do i=1,size(boundint)
	j = boundint(i)
	K(j,j) = -(1 + 2*Sx + 2*Sy)
    K(j,j-ny) = Sx
    K(j,j+ny) = Sx
    K(j,j-1) = Sy
    K(j,j+1) = Sy
enddo

j=0
do i=1,size(upbound)
	j = upbound(i)
    K(j,j) = 1.0D0;
!     Ku(i,i) = 1;
!     Ku(i,i) = 1/dy;
!     Ku(i,i+1) = -1/dy;
enddo
j=0
do i=1,size(dobound)
	j = dobound(i)
    K(j,j) = 1.0D0;
!     Ku(i,i) = 1;
!     Ku(i,i) = -1/dy;
!     Ku(i,i-1) = 1/dy;
enddo
j=0
do i=1,size(lebound)
	j = lebound(i)    
	K(j,j) = 1.0D0;
!     Kv(i,i) = 1;
!     Kv(i,i) = -1/dx;
!     Kv(i,i+n) = 1/dx;
enddo
j=0
do i=1,size(ribound)
	j = ribound(i)    
	K(j,j) = 1.0D0;
!     Kv(i,i) = 1;
!     Kv(i,i) = 1/dy;
!     Kv(i,i-n) = -1/dy;
enddo

write(*,*) Sx, Sy

endsubroutine peregrin


