program day2
use random
use step
implicit none
integer, parameter :: i4=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(15)
integer, parameter :: r8=selected_real_kind(15,9)

integer(kind=i4), parameter :: L=6
integer(kind=i4), parameter :: Lt=50
real(kind=r8), parameter :: c=-0.2
real(kind=r8), parameter :: m=938.92 !MeV
real(kind=r8), parameter :: a=0.01 !MeV^(-1)
real(kind=r8), parameter :: at=0.01 !MeV^(-1)
real(kind=r8), parameter :: alphat=at/a
real(kind=r8), allocatable :: s(:,:)
integer(kind=i4) :: nx,nt
real(kind=r8), dimension(0:L-1,0:Lt) :: vrel
real(kind=r8), dimension(0:Lt) :: overlap
irn=17
allocate(s(0:L-1,0:Lt-1))
call setstep(L,Lt,m,a,at)

!calculate exactly
do nx = 0,L-1
   vrel(nx,0) = 1.0_r8/sqrt(1.0_r8*L)
enddo
overlap(0) = 1.0_r8

do nt=1,Lt
   do nx=0,L-1
   vrel(nx,nt)=vrel(nx,nt-1)
   vrel(nx,nt)=vrel(nx,nt) &
      +1.0_r8/m*vrel(mod(nx-1+L,L),nt-1) &
      +1.0_r8/m*vrel(mod(nx+1,L),nt-1) &
      -2.0_r8/m*vrel(nx,nt-1)
   vrel(nx,nt)=vrel(nx,nt) &
      +1.0_r8/(4.0_r8*m**2.0_r8)*vrel(mod(nx-2+L,L),nt-1) &
      +1.0_r8/(4.0_r8*m**2.0_r8)*vrel(mod(nx+2,L),nt-1) &
      -1.0_r8/m**2.0_r8*vrel(mod(nx-1+L,L),nt-1) &
      -1.0_r8/m**2.0_r8*vrel(mod(nx-1+L,L),nt-1) &
      +6.0_r8/(4.0_r8*m**2.0_r8)*vrel(nx,nt-1)
   enddo
   vrel(0,nt) = vrel(0,nt) - c*alphat*vrel(0,nt-1)
   overlap(nt) = 0.0_r8
   do nx = 0,L-1
      overlap(nt) = overlap(nt) + vrel(nx,nt)/sqrt(1.0_r8*L)
   enddo
enddo

do nt = 1,Lt
   write(*,*)'nt = ',nt,'energy = ', -log(overlap(nt)/overlap(nt-1))/at,'(MeV)'
enddo
write(*,*)
write(*,*) 'Now try with Monte Carlo'
end program day2
