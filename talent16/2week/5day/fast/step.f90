module step
   use random
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   logical, parameter :: wrt=.true. ! writes results out every step if true
   real(kind=r8), private, parameter :: delta=1.0
   integer(kind=i4), parameter :: neq=1000
   integer(kind=i4), parameter :: nav=10000000
   integer(kind=i4), parameter :: nsteps=1000 !steps/print
   integer(kind=i4), parameter :: L=10
   integer(kind=i4), parameter :: Lt=10
   real(kind=r8), parameter :: M2=1.0
   real(kind=r8), parameter :: lambda=0.0
   real(kind=r8), parameter :: pi=3.14159265358979323846
   integer(kind=i8) :: irn
contains
!   subroutine setstep(Lin,Ltin)
!   integer(kind=i4) :: Lin,Ltin
!   end subroutine

   subroutine metropolis(phi,SEold,accepted,niter)
   real(kind=r8), dimension(0:L-1,0:L-1,0:L-1,0:Lt-1) :: phi
   real(kind=r8) :: phinew,phiold
   real(kind=r8) :: ratio,acceptance,rand1
   integer(kind=i4) :: mx,my,mz,mt
   real(kind=r8) :: act
   integer(kind=i4) :: niter
   real(kind=r8) :: accepted,SEold,SEnew
   do mx=1,L-1
      do my=1,L-1
         do mz=1,L-1
            do mt=1,Lt-1
               phiold=phi(mx,my,mz,mt)
               phinew=phi(mx,my,mz,mt)+delta*(ran(irn)-0.5_r8)
               SEnew=SEold
               SEnew = SEnew + &
                  lambda/24.0_r8*(phinew**4.0_r8-phiold**4.0_r8)
               SEnew = SEnew + &
                  0.5_r8*(8.0_r8+M2)*(phinew**2.0_r8-phiold**2.0_r8)
               SEnew = SEnew - &
                 (phi(mod(mx+1,L),my,mz,mt)+phi(mod(mx-1+L,L),my,mz,mt))*(phinew-phiold)
               SEnew = SEnew - &
                 (phi(mx,mod(my+1,L),mz,mt)+phi(mx,mod(my-1+L,L),mz,mt))*(phinew-phiold)
               SEnew = SEnew - &
                 (phi(mx,my,mod(mz+1,L),mt)+phi(mx,my,mod(mz-1+L,L),mt))*(phinew-phiold)
               SEnew = SEnew - &
                 (phi(mx,my,mz,mod(mt+1,Lt))+phi(mx,my,mz,mod(mt-1+Lt,Lt)))*(phinew-phiold)
               niter=niter+1_i4
               if (ran(irn).lt.exp(-SEnew+SEold)) then
                  accepted = accepted + 1.0_r8
                  phi(mx,my,mz,mt) = phinew
                  SEold = SEnew
              endif
            enddo
         enddo
      enddo
   enddo
   end subroutine metropolis

   function SE(phi)
   real(kind=r8), dimension(0:L-1,0:L-1,0:L-1,0:Lt-1) :: phi
   real(kind=r8) :: SE
   integer(kind=i4) :: nx,ny,nz,nt
   SE=0.0_r8
   do nx=0,L-1
      do ny=0,L-1
         do nz=0,L-1
            do nt=0,Lt-1
               SE = SE + phi(nx,ny,nz,nt)*(-phi(mod(nx+1,L),ny,nz,nt) &
                  -phi(nx,mod(ny+1,L),nz,nt) &
                  -phi(nx,ny,mod(nz+1,L),nt) &
                  -phi(nx,ny,nz,mod(nt+1,Lt)) &
                  +0.5_r8*(8.0_r8+M2)*phi(nx,ny,nz,nt) &
                  +lambda/24.0_r8*phi(nx,ny,nz,nt)**3.0_r8)
            enddo
         enddo
      enddo
   enddo
   end function SE
end module step
