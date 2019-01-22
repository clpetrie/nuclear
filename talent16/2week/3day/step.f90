module step
   use random
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), parameter :: delta=0.1
   integer(kind=i8) :: irn
   integer(kind=i4) :: Lt
   real(kind=r8), dimension(0:2,0:2) :: M
contains
   subroutine setstep(MMin,Ltin)
   real(kind=r8), dimension(:,:) :: MMin
   integer(kind=i4) :: Ltin
   M=MMin
   Lt=Ltin
   end subroutine

   subroutine metropolis(x,f,f1count,f2count)
   integer(kind=i8), dimension(0:Lt) :: x
   integer(kind=i8) :: xnew
   real(kind=r8) :: ratio,acceptance,rand1
   integer(kind=i4) :: n,nt
   real(kind=r8) :: fnew,f
   real(kind=r8) :: f1count,f2count
   nt = mod(int(ran(irn)*(Lt-1)),Lt-1)+1 !random from 1 to 19
   if (ran(irn) .lt. 0.5D0) then
      xnew = mod(x(nt)+1,3)
   else
      xnew = mod(x(nt)+2,3)
   endif
   fnew = f/(M(x(nt+1),x(nt))*M(x(nt),x(nt-1))) &
      *M(x(nt+1),xnew)*M(xnew,x(nt-1))
   ratio=fnew/f
   if (ran(irn).lt.ratio) then
      x(nt)=xnew
      f=fnew
   endif
   f1count=f1count+1.0_r8
   if (x(Lt-1).eq.0.0_r8 .and. x(Lt-2).eq.0.0_r8) then
      f2count=f2count+1.0_r8/(M(0,0)*M(0,0))
   endif
   end subroutine metropolis
end module step
