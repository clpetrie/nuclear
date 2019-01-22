program day1
use mymod
implicit none
integer, parameter :: np=10
integer, parameter :: steps=1000 !steps/block
integer, parameter :: method=2 ! 1=rejection, 2=metropolis
integer, parameter :: nb=100 !number of blocks
integer :: test
real, dimension(np) :: x,xsave
real, dimension(nb) :: v1,v2
real :: v1exp,v1err,v2exp,v2err
integer :: i,j,accepted,k

call random_seed

do k=1,nb
   do i=1,10
      call random_number(x(i))
   enddo

   v1(k)=0
   v2(k)=0
   accepted=0
   do i=1,steps
      if (method.eq.1) then
         call rejection(x)
      else if (method.eq.2) then
         xsave=x
         call metropolis(x)
         if (any(x.ne.xsave)) accepted=accepted+1
      endif
      v1(k)=v1(k)+pot1(x)
      v2(k)=v2(k)+pot2(x)
   enddo
   v1(k)=v1(k)/steps
   v2(k)=v2(k)/steps
   !write(*,*) 'acceptance ratio=',float(accepted)/steps
enddo
v1exp=sum(v1)/nb
v1err=sqrt((sum((v1-v1exp)**2)/nb)/nb)
v2exp=sum(v2)/nb
v2err=sqrt((sum((v2-v2exp)**2)/nb)/nb)
write(*,'(A4,E11.5,A5,E11.5)') 'V1 = ',v1exp,' +/- ',v1err
write(*,'(A4,E11.5,A5,E11.5)') 'V2 = ',v2exp,' +/- ',v2err

end program day1
