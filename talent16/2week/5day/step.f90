module step
   use random
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), parameter :: mtype=2 !1=move one, 2=move all, 3=move time slice
   logical, parameter :: wrt=.true. ! writes results out every step if true
   real(kind=r8), private, parameter :: delta1=1.5
   real(kind=r8), private, parameter :: deltaall=0.01
   real(kind=r8), private, parameter :: deltats=0.05
   integer(kind=i4), parameter :: neq=1000
   integer(kind=i4), parameter :: nav=10000000
   integer(kind=i4), parameter :: nsteps=1000 !steps/print
   integer(kind=i4), parameter :: L=10
   integer(kind=i4), parameter :: Lt=10
   real(kind=r8), parameter :: M2=-10.0
   real(kind=r8), parameter :: lambda=24.0
   real(kind=r8), parameter :: a=1.0_r8
   real(kind=r8), parameter :: pi=3.14159265358979323846
   integer(kind=i8) :: irn
contains
!   subroutine setstep(Lin,Ltin)
!   integer(kind=i4) :: Lin,Ltin
!   end subroutine

   subroutine metropolis1(phi,SEold,accepted,niter)
   real(kind=r8), dimension(0:L-1,0:L-1,0:L-1,0:Lt-1) :: phi
   real(kind=r8), dimension(0:L-1,0:L-1,0:L-1,0:Lt-1) :: phinew
   real(kind=r8) :: ratio,acceptance,rand1
   integer(kind=i4) :: nx,ny,nz,nt
   real(kind=r8) :: act
   integer(kind=i4) :: xs,ys,zs,ts,niter !slice that get changed
   real(kind=r8) :: accepted,SEold,SEnew
   xs=nint(ran(irn)*10.0_r8-0.5_r8) !integer from 0 to 9
   ys=nint(ran(irn)*10.0_r8-0.5_r8) !integer from 0 to 9
   zs=nint(ran(irn)*10.0_r8-0.5_r8) !integer from 0 to 9
   ts=nint(ran(irn)*10.0_r8-0.5_r8) !integer from 0 to 9
   phinew=phi
   phinew(xs,ys,zs,ts)=phi(xs,ys,zs,ts)+delta1*(ran(irn)-0.5_r8)
   SEnew=SE(phinew)
   ratio=exp(-SEnew+SEold)
   niter=niter+1_i4
   if (ran(irn).lt.ratio) then
      phi=phinew
      accepted=accepted+1.0_r8
      SEold=SEnew
   endif
   end subroutine metropolis1

   subroutine metropolisall(phi,SEold,accepted,niter)
   real(kind=r8), dimension(0:L-1,0:L-1,0:L-1,0:Lt-1) :: phi
   real(kind=r8), dimension(0:L-1,0:L-1,0:L-1,0:Lt-1) :: phinew
   real(kind=r8) :: ratio,acceptance,rand1
   integer(kind=i4) :: nx,ny,nz,nt
   real(kind=r8) :: act
   integer(kind=i4) :: niter
   real(kind=r8) :: accepted,SEold,SEnew
   phinew=phi
   do nx=0,L-1
      do ny=0,L-1
         do nz=0,L-1
            do nt=1,Lt-1
               phinew(nx,ny,nz,nt)=phi(nx,ny,nz,nt)+deltaall*(ran(irn)-0.5_r8)
            enddo
         enddo
      enddo
   enddo
   SEnew=SE(phinew)
   ratio=exp(-SEnew+SEold)
   niter=niter+1_i4
   if (ran(irn).lt.ratio) then
      phi=phinew
      accepted=accepted+1.0_r8
      SEold=SEnew
   endif
   end subroutine metropolisall

   subroutine metropolists(phi,SEold,accepted,niter)
   real(kind=r8), dimension(0:L-1,0:L-1,0:L-1,0:Lt-1) :: phi
   real(kind=r8), dimension(0:L-1,0:L-1,0:L-1,0:Lt-1) :: phinew
   real(kind=r8) :: ratio,acceptance,rand1
   integer(kind=i4) :: nx,ny,nz,ts
   real(kind=r8) :: act
   integer(kind=i4) :: niter
   real(kind=r8) :: accepted,SEold,SEnew
   ts=nint(ran(irn)*10.0_r8-0.5_r8) !integer from 0 to 9
   phinew=phi
   do nx=0,L-1
      do ny=0,L-1
         do nz=0,L-1
            phinew(nx,ny,nz,ts)=phi(nx,ny,nz,ts)+deltats*(ran(irn)-0.5_r8)
         enddo
      enddo
   enddo
   SEnew=SE(phinew)
   ratio=exp(-SEnew+SEold)
   niter=niter+1_i4
   if (ran(irn).lt.ratio) then
      phi=phinew
      accepted=accepted+1.0_r8
      SEold=SEnew
   endif
   end subroutine metropolists

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
