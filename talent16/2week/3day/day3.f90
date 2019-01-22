program day3
use random
use step
implicit none

integer, parameter :: i4=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(15)
integer, parameter :: r8=selected_real_kind(15,9)
integer(kind=i4), parameter :: nstepseq=100000
integer(kind=i4), parameter :: nsteps=1000000
integer(kind=i8), allocatable :: x(:)
integer(kind=i4) :: s,n,i
real(kind=r8) :: ave,err,f1count,f2count,f
irn=17
M=reshape((/1.1,0.1,0.1,0.1,0.8,0.1,0.1,0.1,0.8/),shape(M))
Lt=20
call setstep(M,Lt)
allocate(x(0:Lt))

!initialize walkers
x(0)=0.0_r8
x(Lt)=0.0_r8
do n=1,Lt-1
   x(n)=nint(ran(irn)*3_i8-0.5_r8)
enddo

f=1.0_r8
do i=1,Lt
   f=f*M(x(i),x(i-1)) !initial configuration
enddo

do s=1,nstepseq
   call metropolis(x,f,f1count,f2count)
enddo

f1count=0_r8
f2count=0_r8

do s=1,nsteps
   call metropolis(x,f,f1count,f2count)
   if (f2count.gt.0) then
      ave=f1count/f2count
      write(*,*) 's,ave=',s,ave
   endif
enddo
ave=nsteps/f2count
err=f1count/f2count

write(*,*) ave,err

end program day3
