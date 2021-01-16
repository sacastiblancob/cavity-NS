!----------------------------------------------------------------------!
!This subroutine builds the vector product of ALHS*PHI for the SEM
!on Quadrilateral Elements for the Shallow Water Equations on the Sphere
!for the Geopotential Semi-Implicit.
!Written by Francis X. Giraldo on 8/99
!           Naval Research Laboratory
!           Global Modeling Section
!           Monterey, CA 93943-5502
!----------------------------------------------------------------------!
      subroutine lhs_gmres(alhs,x,ain,npoin)

      real, dimension(npoin) :: alhs,x
      real, dimension(npoin,npoin) :: ain
      do i=1,npoin
        aout = 0.
        do j=1,npoin
          aout = aout + ain(i,j)*x(j) 
        enddo
        alhs(i) = aout
      enddo

      end subroutine lhs_gmres



