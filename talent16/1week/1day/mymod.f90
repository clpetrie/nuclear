module mymod
implicit none
real, parameter :: beta=1.0
real, parameter :: b=0.1 !something small
real, parameter :: V0=1.0
contains
   function pot1(x)
   real, parameter :: pi=3.141592654
   real, dimension(:) :: x
   real :: pot1
   integer :: i
   pot1=0
   do i=1,size(x)
      pot1=pot1+1+cos(20*pi*x(i))
   enddo
   end function pot1

   function pot2(x)
   real, dimension(:) :: x
   real :: pot2
   integer :: i,j
   pot2=0
   do i=1,size(x)-1
      do j=i+1,size(x)
         if (abs(x(i)-x(j)).lt.0.05) pot2=pot2+V0
      enddo
   enddo
  end function pot2

   function W(x) !this is what you sample from
   real, dimension(:) :: x
   real :: W
   W=exp(-beta*(pot1(x)+pot2(x)))
   end function W

   subroutine rejection(x)
   real :: vmax
   real, dimension(:) :: x
   real :: ratio,rand1
   integer :: i
   logical :: accept
   vmax=2+0.5*size(x)*(size(x)-1)
   accept=.false.
   do while (.not.accept)
      do i=1,10
         call random_number(x(i))
      enddo
      ratio=W(x)/vmax
      call random_number(rand1)
      if (rand1.lt.ratio) accept=.true.
   enddo
   end subroutine rejection

   subroutine metropolis(x)
   real, dimension(:) :: x
   real, dimension(size(x)) :: xnew
   real :: rand1,ratio
   integer :: i
   do i=1,size(x)
      call random_number(rand1)
      xnew(i)=x(i)+b*(rand1-0.5)
      if (xnew(i)<0) xnew(i)=xnew(i)+1 !periodic bc
      if (xnew(i)>1) xnew(i)=xnew(i)-1 !periodic bc
   enddo
   ratio=W(xnew)/W(x)
   call random_number(rand1)
   if (rand1.lt.ratio) x=xnew
   end subroutine metropolis
end module mymod
