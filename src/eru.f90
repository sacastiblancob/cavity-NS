subroutine eru(x,y)

! ################################################################################################################################ !
! GRID GENERATOR
! ################################################################################################################################ !
!As eru, the creator of Middle Earth
USE gandalf
real, dimension(1:nx) :: x
real, dimension(1:ny) :: y

x(1) = xmin
do i=2,nx
	x(i) = x(i-1) + dx
end do

y(1) = ymin
do i=2,ny
	y(i) = y(i-1) + dy
end do

endsubroutine eru
