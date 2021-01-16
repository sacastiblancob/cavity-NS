!-----------------------------------------------------------------------
!      GMRES Solver: 
!             Original Version by P.F. Fischer, Argonne National Lab
!             Modified Version by F.X. Giraldo, Naval Research Lab
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine solve_gmres(x,b,ain,npoin)
!     include 'param.h'

      integer, parameter :: mmax = 300
      !Solver arrays
      real ain(npoin,npoin)

      real adiag(npoin), x(npoin), b(npoin), lambda
      real w(npoin), v(npoin,mmax)

      real h(mmax,mmax), s(mmax), c(mmax), gamma(mmax)
      real h_loc, h_glob, twonorm, glsc2
      logical ifdone

      !set machine tolerances
      solver_tol = 1.e-8
      n=npoin
      one = 1.
      eps = 1.e-20
      if (one+eps == one) eps = 1.e-14
      if (one+eps == one) eps = 1.e-7
      if (one+eps == one) eps = 1.e-16
!      eps=1.0e-7
      eps=solver_tol

!-Intrinsic GMRES constants
      gamma1 = -1
      
!-Set up preconditioner matrix (PD: See Heath's book.
!-Use Diagonal preconditioner)
      do i=1,npoin
        adiag(i) = -ain(i,i)
      enddo

      !Initialize
      call rzero(x,n)
 
      !Solve M*w=b
      call precon(w,b,adiag,n)

      !Compute L2 Norm
      gamma(1) = twonorm(w,n)

      if (gamma(1) == 0) return
      t = 1./gamma(1)
      call cmult2(v,w,t,n)

      if (gamma1 < 0) then
         etol  = gamma(1)*eps
         etol2 = gamma(1)*eps*eps
      else    
         etol  = gamma1*eps     
         etol2 = gamma1*eps*eps
      end if

      !Begin Arnoldi w/ modified gram-schmidt
      maxiter = mmax
      do k=1,maxiter
         write(*,*) k
         !construct v(:,k+1)=A*v(:,k)
         call lhs_gmres(v(:,k+1),v(:,k),ain,npoin)


         !Solve M*w=v(:,k+1)
         call precon(w,v(:,k+1),adiag,n)

         do j=1,k
            h(j,k) = glsc2(w,v(:,j),n)
            t      = -h(j,k)
            call add2s2(w,v(:,j),t,n)
         end do !j

         !Compute L2 Norm
         h(k+1,k) = twonorm(w,n)

         if (abs(h(k+1,k)) <  etol2) then
            ifdone=.true.
         else
            ifdone=.false.
            t = 1./h(k+1,k)
            call cmult2(v(:,k+1),w,t,n)
         end if

         !apply Given's rotations to new column of H
         do i=1,k-1
            t = h(i,k)
            h(i  ,k) =  c(i)*t + s(i)*h(i+1,k)
            h(i+1,k) = -s(i)*t + c(i)*h(i+1,k)
         end do !i
         den        =  sqrt(  h(k,k)*h(k,k) + h(k+1,k)*h(k+1,k)  )
         c(k)       =  h(k  ,k) / den
         s(k)       =  h(k+1,k) / den
         h(k,k)     =  c(k)*h(k,k)+s(k)*h(k+1,k)
         gamma(k+1) = -s(k)*gamma(k)
         gamma(k  ) =  c(k)*gamma(k)
         if (ifdone .or. abs(gamma(k+1)) < etol) exit
      end do !k

      !if we're here, we exceeded the max iteration, reduce k by 1
      if (k > maxiter) then
         print*,' MAXITER exceeded'
         k = k-1
      end if

      !Compute solution via back substitution
!     write(*,'("     k Tol     = ",i5,e16.8)')k,gamma(k+1)/gamma1
      kiter1=kiter0
      kiter0=k
      kiter=kiter + k
      do i=k,1,-1
         t = gamma(i)
         do j=k,i+1,-1
            t = t - h(i,j)*c(j)
         end do !j
         c(i) = t/h(i,i)
      end do !i

      !Sum up Arnoldi vectors
      do i=1,k
         call add2s2(x,v(:,i),c(i),n)
      end do !i

      end subroutine solve_gmres
!-----------------------------------------------------------------------
      subroutine add2s2(a,b,c1,n)
!     include 'param.h'
      real a(n), b(n)

      do i=1,n
        a(i)=a(i)+c1*b(i)
      end do  

      end subroutine add2s2
!-----------------------------------------------------------------------
      function twonorm(x,n)
!     include 'param.h'
      real x(n), twonorm

      tscal = x(1)*x(1)
      do i=2,n
         tscal = tscal+x(i)*x(i)
      end do
      twonorm=tscal
      if (twonorm.gt.0) twonorm = sqrt(twonorm)

      end function twonorm
!-----------------------------------------------------------------------
      function glsc2(x,y,n)
!     include 'param.h'
      real x(n), y(n), glsc2

      tscal = x(1)*y(1)
      do i=2,n
         tscal = tscal+x(i)*y(i)
      end do
      glsc2=tscal

      end function glsc2
!-----------------------------------------------------------------------
      subroutine cmult2(a,b,c,n)
!     include 'param.h'
      real a(n), b(n), c

      do i = 1, n
         a(i) = b(i)*c
      end do

      end subroutine cmult2
!-----------------------------------------------------------------------
      subroutine rzero(a,n)
!     include 'param.h'
      real a(n)

      do i = 1, n
         a(i) = 0.0
      end do   

      end subroutine rzero
!-----------------------------------------------------------------------
      subroutine precon(y,x,a,n)
!     include 'param.h'
      real y(n), x(n), a(n)

      do i=1,n
         y(i)=x(i)/a(i)
      end do !i

      end subroutine precon
!-----------------------------------------------------------------------
      subroutine cmult(a,c1,n)
      real a(n), c1

      do i=1,n
        a(i)=a(i)*c1
      end do  

      end subroutine cmult
!-----------------------------------------------------------------------
      subroutine col2(a,b,n)
      real a(n), b(n)

      do i=1,n
        a(i)=a(i)*b(i)
      end do  

      end subroutine col2
!-----------------------------------------------------------------------
      subroutine invcol2(a,b,n)
      real a(n), b(n)

      do i=1,n
        a(i)=a(i)/b(i)
      end do  

      end subroutine invcol2
!-----------------------------------------------------------------------
      subroutine add2(a,b,n)
      real a(n), b(n)

      do i=1,n
        a(i)=a(i)+b(i)
      end do  

      end subroutine add2
!-----------------------------------------------------------------------
      subroutine chsign(a,n)
      real a(n)

      do i=1,n
        a(i)=-a(i)
      end do  

      end subroutine chsign
!-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      real a(n), b(n)

      do i=1,n
        a(i)=b(i)
      end do  

      end subroutine copy
!-----------------------------------------------------------------------
      subroutine sub2(a,b,n)
      real a(n), b(n)

      do i=1,n
        a(i)=a(i)-b(i)
      end do  

      end subroutine sub2
!-----------------------------------------------------------------------
      function vlsc2(x,y,n)
      real x(n), y(n), vlsc2

      tscal = x(1)*y(1)
      do i=2,n
         tscal = tscal+x(i)*y(i)
      end do
      vlsc2 = tscal

      end function vlsc2
!-----------------------------------------------------------------------
      subroutine proj(xk,n,nk,b,nvectors)  ! L2 projection of b onto xk
!     include 'param.h'
      real xk(n,0:nvectors,2), b(n), alpha(nvectors), work(nvectors)

      if (nk.eq.0) then
         call rzero(xk(:,0,1),n)
         call rzero(xk(:,0,2),n)
         return
      endif

      tb = twonorm(b,n)
      tk = twonorm(xk(:,1,2),n)
!      if (irank==0) write(6,*) '2-norm:',nk,tb,tk
      do k=1,nk
        work(k) = vlsc2(xk(:,k,2),b,n)
      end do  
      alpha=work
!      if (irank==0) write(6,1) nk,(alpha(k),k=1,nk)
!    1 format('proj:',i4,1p10e10.2)

      call cmult2 (xk(:,0,1),xk(:,1,1),alpha(1),n)
      call cmult2 (xk(:,0,2),xk(:,1,2),alpha(1),n)

      do k=2,nk
         call add2s2 (xk(:,0,1),xk(:,k,1),alpha(k),n)
         call add2s2 (xk(:,0,2),xk(:,k,2),alpha(k),n)
      enddo

c     do i=1,nk
c     do j=1,nk
c        alpha(j) = glsc2(xk(:,i,2),xk(:,j,2),ip,n)
c     enddo
c     if (irank==0) write(6,2) i,(alpha(j),j=1,nk)
c     enddo
c   2 format('Aij:',i4,1p10e10.2)

      end subroutine proj
!-----------------------------------------------------------------------
