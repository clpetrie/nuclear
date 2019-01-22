module mymod
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8) :: delta,alpha,omega,dt,etrial
   integer(kind=i8) :: irnin2
contains
   subroutine setmod(deltain,alphain,omegain,dtin,etrialin)
   real(kind=r8) :: deltain,alphain,omegain,dtin,etrialin
   delta=deltain
   alpha=alphain
   omega=omegain
   dt=dtin
   etrial=etrialin
   end subroutine setmod

   function psit(x)
   real(kind=r8) :: x
   real(kind=r8) :: psit
   psit=exp(-0.5*alpha*x**2)
   end function psit

   function eloc(x)
   real(kind=r8) :: x
   real(kind=r8) :: eloc
   eloc=0.5*alpha*(1-alpha*x**2)+0.5*omega**2*x**2
   end function eloc

   subroutine metropolis(x,acceptance)
   real(kind=r8), dimension(:) :: x
   real(kind=r8), allocatable :: xnew(:)
   real(kind=r8) :: rand1,ratio,acceptance
   integer(kind=i4) :: n,accepted,nwalk
   nwalk=size(x)
   allocate(xnew(nwalk))
   accepted=0
   do n=1,nwalk
!      call random_number(rand1)
      rand1=ran(irnin2)
      xnew(n)=x(n)+delta*(rand1-0.5)
      ratio=(psit(xnew(n))/psit(x(n)))**2
!      call random_number(rand1)
      rand1=ran(irnin2)
      if (rand1.lt.ratio) then
         x(n)=xnew(n)
         accepted=accepted+1
      endif
   enddo
   acceptance=float(accepted)/nwalk
   end subroutine metropolis

   subroutine step(x,tau)
   real(kind=r8) :: x,randg,tau
!   call rgauss(randg)
   randg=rgauss(irnin2)
   x=x-dt*alpha*x+sqrt(dt)*randg
   tau=tau+dt
   end subroutine step

   subroutine branch(x,n,nwalk,nwalk0)
   real(kind=r8), dimension(:) :: x
   integer(kind=i4) :: n,nwalk,nw,i,nwalk0
   real(kind=r8) :: weight,rand1
   if (nwalk.gt.1.1*nwalk0 .or. nwalk.lt.0.9*nwalk0) then
      weight=exp(-dt*(eloc(x(n))-etrial))*nwalk0/nwalk !weight with population control
   else
      weight=exp(-dt*(eloc(x(n))-etrial)) !weight
   endif
!   call random_number(rand1)
   rand1=ran(irnin2)
   nw=weight+rand1
   if (nw.eq.0) then
      do i=n,nwalk-1
         x(i)=x(i+1)
      enddo
      nwalk=nwalk-1
      elseif (nw.gt.1) then
      do i=1,nw
         x(nwalk+i)=x(n)
         nwalk=nwalk+1
      enddo
   endif
   end subroutine branch

   subroutine stat(obs,ave,err)
   real(kind=r8), dimension(:) :: obs
   real(kind=r8), intent(out) :: ave,err
   integer(kind=i4) :: num
   num=size(obs)
   ave=sum(obs)/num
   err=sqrt(sum((obs-ave)**2)/num)
   end subroutine stat

!   subroutine rgauss(rout)
!   implicit none
!   real(kind=r8) :: pi,rand1,rand2,rout
!   pi=4.0*atan(1.0)
!   call random_number(rand1)
!   call random_number(rand2)
!   rout=sqrt(-2.0*log(rand1))*cos(2.0*pi*rand2)
!   end subroutine rgauss

   function ran(irn)
   integer(kind=i8),  parameter :: mask24 = ishft(1_i8,24)-1
   integer(kind=i8),  parameter :: mask48 = ishft(1_i8,48_i8)-1_i8
   real(kind=r8),  parameter :: twom48=2.0_r8**(-48)
   integer(kind=i8),  parameter :: mult1 = 44485709377909_i8
   integer(kind=i8),  parameter :: m11 = iand(mult1,mask24)
   integer(kind=i8),  parameter :: m12 = iand(ishft(mult1,-24),mask24)
   integer(kind=i8),  parameter :: iadd1 = 96309754297_i8
   integer(kind=i8) :: irn
   real(kind=r8) :: ran
   integer(kind=i8) :: is1,is2
   is2=iand(ishft(irn,-24),mask24)
   is1=iand(irn,mask24)
   irn=iand(ishft(iand(is1*m12+is2*m11,mask24),24)+is1*m11+iadd1,mask48)
   ran=ior(irn,1_i8)*twom48
   return
   end function ran

   function rgauss(irn)
   real(kind=r8) :: pi
   real(kind=r8) :: rgauss,x1,x2
   integer(kind=i8) :: irn
   x1=ran(irn)
   x2=ran(irn)
   pi=4.0_r8*atan(1.0_r8)
   rgauss=sqrt(-2.0_r8*log(x1))*cos(2.0_r8*pi*x2)
   return
   end function rgauss
end module mymod
