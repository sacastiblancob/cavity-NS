module gandalf

! ################################################################################################################################ !
! NUMERICAL PARAMETERS
! ################################################################################################################################ !
!As gandalf, the magician that move all behind threads

! Number of nodes in both directions X and Y
integer, parameter :: nx = 6
integer, parameter :: ny = 6

! Dimensions of the computational space (grid)
real, parameter :: xmin = 0.0
real, parameter :: xmax = 1.0
real, parameter :: ymin = 0.0
real, parameter :: ymax = 1.0

! Space differentials
real :: dx = (ymax-ymin)/(nx-1)
real :: dy = (ymax-ymin)/(ny-1)

! Grid
integer, dimension(1:nx*ny) :: pos

! Time discretization
real, parameter :: to = 0.0
real, parameter :: tf = 2.0

! Numerical constants
real :: CFL = 0.25

endmodule gandalf

