program day4
use random
use step
implicit none

integer, parameter :: i4=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(15)
integer, parameter :: r8=selected_real_kind(15,9)
integer(kind=i8), allocatable :: x(:)
real(kind=r8), allocatable :: anst(:)
integer(kind=i4) :: i,s,nx,ny,nz,nt,kx,ky,kz,kt
integer(kind=i4) :: mx,my,mz,mt
real(kind=r8), allocatable :: phi(:,:,:,:)
real(kind=r8), allocatable :: expbin(:)
real(kind=r8) :: avecount,accepted
integer(kind=i4) :: niter
irn=17
allocate(anst(0:L-1))
allocate(phi(0:L-1,0:L-1,0:L-1,0:Lt-1))
allocate(expbin(0:L-1))
!call getstep(L,Lt,M)

! Exact solution
anst=0.0_r8
do nx=0,L-1
   do kx=0,L-1
      do ky=0,L-1
         do kz=0,L-1
            do kt=0,Lt-1
               anst(nx)=anst(nx)+1.0_r8/(2.0_r8*(1-cos(2.0_r8*pi*kx/L) &
                  +1-cos(2.0_r8*pi*ky/L) &
                  +1-cos(2.0_r8*pi*kz/L) &
                  +1-cos(2.0_r8*pi*kt/Lt))+M**2)*cos(2.0_r8*pi*kx/L*nx)
            enddo
         enddo
      enddo
   enddo
   anst(nx)=anst(nx)*1.0_r8/(L**3*Lt)
enddo

! Monte Carlo solution
do nx=0,L-1
   do ny=0,L-1
      do nz=0,L-1
         do nt=0,Lt-1
            phi(nx,ny,nz,nt)=2.0_r8*(ran(irn)-0.5) !initialize from -1 to 1
         enddo
      enddo
   enddo
enddo

! equilibrium blocks
accepted=0.0_r8
do i=1,neq
   call metropolis(phi,accepted)
   if (i.eq.neq) then
      write(*,*)
      write(*,*) 'acceptance=',accepted/neq
   endif
enddo

write(*,*)
write(*,*) 'equilibrium done!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

expbin=0.0_r8
avecount=0.0_r8
niter=0_i4
do i=1,nav
   accepted=0.0_r8
   do s=1,nsteps
      call metropolis(phi,accepted)
      niter=niter+1_i4
   enddo
   write(*,*) 'iteration=',niter
   write(*,*) 'acceptance=',accepted/nsteps
   do nx = 0,L-1
      do mt = 0,Lt-1
         do mz = 0,L-1
            do my = 0,L-1
               do mx = 0,L-1
                  expbin(nx)=expbin(nx)+phi(mod(mx+nx,L),my,mz,mt)*phi(mx,my,mz,mt)/(L**3)
               enddo
            enddo
         enddo
      enddo
      avecount=avecount+1.0_r8
      write(*,*) nx,expbin(nx)/avecount,anst(nx)
   enddo
   write(*,*)
   write(*,*)
enddo

!expbin=expbin/nav
end program day4
