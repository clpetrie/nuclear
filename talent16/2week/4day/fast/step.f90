module step
   use random
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, parameter :: delta=1.5
   integer(kind=i4), parameter :: neq=1000
   integer(kind=i4), parameter :: nav=100000
   integer(kind=i4), parameter :: nsteps=1000 !steps/block
   integer(kind=i4), parameter :: L=10
   integer(kind=i4), parameter :: Lt=10
   real(kind=r8), parameter :: M=1.0
   real(kind=r8), parameter :: a=1.0_r8
   real(kind=r8), parameter :: pi=3.14159265358979323846
   integer(kind=i8) :: irn
contains
!   subroutine setstep(Lin,Ltin)
!   integer(kind=i4) :: Lin,Ltin
!   end subroutine

   subroutine metropolis(phi,SEold,accepted,niter)
   real(kind=r8), dimension(0:L-1,0:L-1,0:L-1,0:Lt-1) :: phi
   real(kind=r8) :: ratio,acceptance,rand1
   integer(kind=i4) :: nx,ny,nz,nt,niter
   real(kind=r8) :: act
   real(kind=r8) :: accepted,SEold,SEnew
   real(kind=r8) :: phinew,phiold
   do nx=0,L-1
      do ny=0,L-1
         do nz=0,L-1
            do nt=0,L-1
               phiold=phi(nx,ny,nz,nt)
               phinew=phi(nx,ny,nz,nt)+delta*(ran(irn)-0.5_r8)
               SEnew=SEold
               SEnew = SEnew + (phinew-phiold)*( &
                  -(phi(mod(nx+1,L),ny,nz,nt)) &
                  -(phi(nx,mod(ny+1,L),nz,nt)) &
                  -(phi(nx,ny,mod(nz+1,L),nt)) &
                  -(phi(nx,ny,nz,mod(nt+1,Lt)))) &
                  +0.5_r8*(8.0_r8+M**2)*(phinew**2.0_r8-phiold**2.0_r8)
               ratio=exp(-SEnew+SEold)
               niter=niter+1_i4
               if (ran(irn).lt.ratio) then
                  phi(nx,ny,nz,nt)=phinew
                  accepted=accepted+1.0_r8
                  SEold=SEnew
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
                  +0.5_r8*(8.0_r8+M**2)*phi(nx,ny,nz,nt))
            enddo
         enddo
      enddo
   enddo
   end function SE
end module step
