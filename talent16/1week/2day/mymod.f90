module mymod
   use random
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8) :: delta,alpha,omega,dt,etrial
   integer(kind=i4) :: npart
   integer(kind=i8) :: irn
contains
   subroutine setmod(deltain,alphain,omegain,dtin,etrialin,npartin)
   real(kind=r8) :: deltain,alphain,omegain,dtin,etrialin
   integer(kind=i4) :: npartin
   delta=deltain
   alpha=alphain
   omega=omegain
   dt=dtin
   etrial=etrialin
   npart=npartin
   end subroutine setmod

   function psit(x)
   real(kind=r8), dimension(:) :: x
   real(kind=r8) :: psit
   integer(kind=i4) :: i
   psit=0
   do i=1,npart
      psit=psit+exp(-0.5*alpha*x(i)**2)
   enddo
   end function psit

   function eloc(x)
   real(kind=r8),dimension(:) :: x
   real(kind=r8) :: eloc
   integer(kind=i4) :: i
   eloc=0
   do i=1,npart
      eloc=eloc+0.5*alpha*(1-alpha*x(i)**2)+0.5*omega**2*x(i)**2
   enddo
   end function eloc

   function egrowth(x,nwalk)
   real(kind=r8), dimension(:,:) :: x
   real(kind=r8) :: egrowth,aw,ew
   real(kind=r8), allocatable :: nweight(:)
   integer(kind=i4) :: n,nwalk
   real(kind=r8), dimension(3) :: en
   allocate(nweight(nwalk))
   do n=1,nwalk
      nweight(n)=exp(-dt*(eloc(x(n,:))-etrial))
   enddo
   call stat(nweight,aw,ew)
   egrowth=etrial-1.0_r8/dt*log(aw)
   end function egrowth

   subroutine metropolis(x,nwalk,acceptance)
   real(kind=r8), dimension(:,:) :: x
   real(kind=r8), allocatable :: xnew(:,:)
   real(kind=r8) :: rand1,ratio,acceptance
   integer(kind=i4) :: i,n,accepted,nwalk
   allocate(xnew(nwalk,npart))
   accepted=0
   do n=1,nwalk
      do i=1,npart
         xnew(n,i)=x(n,i)+delta*(ran(irn)-0.5)
      enddo
      ratio=(psit(xnew(n,:))/psit(x(n,:)))**2
      if (ran(irn).lt.ratio) then
         x(n,:)=xnew(n,:)
         accepted=accepted+1
      endif
   enddo
   acceptance=float(accepted)/nwalk
   end subroutine metropolis

   subroutine step(x)
   real(kind=r8), dimension(:) :: x
   real(kind=r8) :: randg
   integer(kind=i4) :: i
   do i=1,npart
      x(i)=x(i)-dt*alpha*x(i)+sqrt(dt)*rgauss(irn)
   enddo
   end subroutine step

   subroutine branch(x,nwalkin,nwalk0)
   real(kind=r8), dimension(:,:) :: x
   real(kind=r8), allocatable :: x2(:,:)
   integer(kind=i4) :: n,nwalk,nw,iw,nwalk0,nt,nwalkin
   real(kind=r8) :: weight,rand1
   allocate(x2(10*nwalk0,npart))
   nt=1
   do n=1,nwalkin
      if (nwalk.gt.1.1*nwalk0 .or. nwalk.lt.0.9*nwalk0) then
         weight=exp(-dt*(eloc(x(n,:))-etrial))*nwalk0/nwalk !weight with population control
      else
         weight=exp(-dt*(eloc(x(n,:))-etrial))
      endif
      nw=weight+ran(irn)
      nwalk=nwalk+nw-1
      do iw=nt,nt+nw
            x2(iw,:)=x(n,:)
            nt=nt+1
      enddo
      x=x2
   enddo
   end subroutine branch

   subroutine stat(obs,ave,err)
   real(kind=r8), dimension(:) :: obs
   real(kind=r8), intent(out) :: ave,err
   real(kind=r8) :: ave2
   integer(kind=i4) :: num
   num=size(obs)
   ave=sum(obs)/num
   ave2=0
   ave2=sum(obs**2)/num
   err=sqrt(abs(ave2-ave**2)/num)
   end subroutine stat
end module mymod
