module mymod
   use random
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), parameter :: dx=0.0005_r8
   real(kind=r8), parameter :: rsmall=1.0E-10_r8
   real(kind=r8) :: delta,dt,etrial,h2o2m,l
   integer(kind=i4) :: npart
   integer(kind=i8) :: irn
contains
   subroutine setmod(deltain,dtin,h2o2min,lin,etrialin,npartin)
   real(kind=r8) :: deltain,dtin,etrialin,h2o2min,lin
   integer(kind=i4) :: npartin
   delta=deltain
   dt=dtin
   h2o2m=h2o2min
   l=lin
   etrial=etrialin
   npart=npartin
   end subroutine setmod

!   function psit(x)
!   real(kind=r8), parameter :: b=2.99
!   real(kind=r8), dimension(:,:) :: x
!   real(kind=r8) :: psit,r
!   integer(kind=i4) :: i,j
!   psit=1
!   do i=1,npart-1
!      do j=i+1,npart
!         r=sqrt(sum((x(:,i)-x(:,j))**2.0_r8))
!         r=sqrt(sum((x(:,i)-x(:,j)-nint((x(:,i)-x(:,j))/l)*l)**2.0_r8))
!         psit=psit*exp(-0.5_r8*(b/r)**5.0_r8)
!         psit=psit*exp(-0.5_r8*(b/r)**5.0_r8-0.5_r8*(b/(l-r))**5.0_r8-0.5_r8*(b/(0.5_r8*l))**5.0_r8)
!      enddo
!   enddo
!   end function psit

   function lpsit(x)
   real(kind=r8), parameter :: b=2.99_r8
   real(kind=r8), dimension(:,:) :: x
   real(kind=r8) :: lpsit,r
   integer(kind=i4) :: i,j
   lpsit=0.0_r8
   do i=1,npart-1
      do j=i+1,npart
         r=sqrt(sum((x(:,i)-x(:,j)-nint((x(:,i)-x(:,j))/l)*l)**2.0_r8))
         if (r.lt.rsmall) r=rsmall
         if (r.gt.0.5_r8*l) cycle
!         lpsit=lpsit+(b/r)**5.0_r8
         lpsit=lpsit+(b/r)**5.0_r8+(b/(l-r))**5.0_r8-2.0_r8*(b/(0.5_r8*l))**5.0_r8
      enddo
   enddo
   lpsit=-lpsit*0.5_r8
   end function lpsit

   function ddlpsit(x,i)
   real(kind=r8), dimension(:,:) :: x
   real(kind=r8), allocatable :: xsave(:,:)
   real(kind=r8) :: ddlpsit
   integer(kind=i4) :: i,ic
   allocate(xsave(3,npart))
   ddlpsit=0.0_r8
   xsave=x
   do ic=1,3
      x(ic,i)=x(ic,i)+dx
      ddlpsit=ddlpsit+lpsit(x)
      x=xsave
      ddlpsit=ddlpsit-2.0_r8*lpsit(x)
      x(ic,i)=x(ic,i)-dx
      ddlpsit=ddlpsit+lpsit(x)
      x=xsave
   enddo
   ddlpsit=ddlpsit/dx**2.0_r8
   end function ddlpsit

   function dlpsit(x,i)
   real(kind=r8), dimension(:,:) :: x
   real(kind=r8), allocatable :: xsave(:,:)
   real(kind=r8),dimension(3) :: dlpsit
   integer(kind=i4) :: i,ic
   allocate(xsave(3,npart))
   dlpsit=0.0_r8
   xsave=x
   do ic=1,3
      x(ic,i)=x(ic,i)+dx
      dlpsit(ic)=dlpsit(ic)+lpsit(x)
      x=xsave
      x(ic,i)=x(ic,i)-dx
      dlpsit(ic)=dlpsit(ic)-lpsit(x)
      x=xsave
   enddo
   dlpsit=dlpsit/(2.0_r8*dx)
   end function dlpsit

   function eloc(x)
   real(kind=r8), parameter :: eps=10.22_r8
   real(kind=r8), parameter :: sigma=2.556_r8
   real(kind=r8), dimension(:,:) :: x
   real(kind=r8), dimension(3) :: eloc !total, kinetic, potential energies
   real(kind=r8) :: r
   integer(kind=i4) :: i,j
   real(kind=r8),dimension(3,npart) :: grad
   eloc=0.0_r8
   do i=1,npart
      grad(:,i)=dlpsit(x,i)
      eloc(2)=eloc(2)-h2o2m*(ddlpsit(x,i)+dot_product(grad(:,i),grad(:,i)))
   enddo
   do i=1,npart-1
      do j=i+1,npart
         r=sqrt(sum((x(:,i)-x(:,j)-nint((x(:,i)-x(:,j))/l)*l)**2.0_r8))
         if (r.lt.rsmall) r=rsmall
         eloc(3)=eloc(3)+4.0_r8*eps*((sigma/r)**12.0_r8-(sigma/r)**6.0_r8)
      enddo
   enddo
   eloc(1)=eloc(2)+eloc(3)
   end function eloc

   function egrowth(x,nwalk)
   real(kind=r8), dimension(:,:,:) :: x
   real(kind=r8) :: egrowth,aw,ew
   real(kind=r8), allocatable :: nweight(:)
   integer(kind=i4) :: n,nwalk
   real(kind=r8), dimension(3) :: en
   allocate(nweight(nwalk))
   do n=1,nwalk
      en=eloc(x(:,:,n))
      nweight(n)=exp(-dt*(en(1)-etrial))
   enddo
   call stat(nweight,aw,ew)
   egrowth=etrial-1.0_r8/dt*log(aw)
   end function egrowth

   subroutine metropolis(x,nwalk,acceptance)
   real(kind=r8), dimension(:,:,:) :: x
   real(kind=r8), allocatable :: xnew(:,:,:)
   real(kind=r8) :: ratio,acceptance
   integer(kind=i4) :: i,ic,n,accepted,nwalk
   allocate(xnew(3,npart,nwalk))
   accepted=0.0_r8
   xnew=0.0_r8
   do n=1,nwalk
      do i=1,npart
         do ic=1,3
            xnew(ic,i,n)=x(ic,i,n)+delta*(ran(irn)-0.5_r8)
            !boundary conditions
            if (xnew(ic,i,n).gt.0.5_r8*l) xnew(ic,i,n)=xnew(ic,i,n)-l
            if (xnew(ic,i,n).lt.-0.5_r8*l) xnew(ic,i,n)=xnew(ic,i,n)+l
         enddo
      enddo
      ratio=exp(2.0_r8*(lpsit(xnew(:,:,n))-lpsit(x(:,:,n))))
      if (ran(irn).lt.ratio) then
         x(:,:,n)=xnew(:,:,n)
         accepted=accepted+1
      endif
   enddo
   acceptance=float(accepted)/nwalk
   end subroutine metropolis

   subroutine step(x)
   real(kind=r8), dimension(:,:) :: x
   integer(kind=i4) :: i,ic
   real(kind=r8), dimension(3,npart) :: grad
   do i=1,npart
      grad(:,i)=dlpsit(x,i)
      do ic=1,3
         x(ic,i)=x(ic,i)+2.0_r8*h2o2m*dt*grad(ic,i)+rgauss(irn)*sqrt(2.0_r8*h2o2m*dt)
         !boundary conditions
         if (x(ic,i).gt.0.5_r8*l) x(ic,i)=x(ic,i)-l
         if (x(ic,i).lt.-0.5_r8*l) x(ic,i)=x(ic,i)+l
      enddo
   enddo
   end subroutine step

   subroutine branch(x,nwalkin,nwalk0)
   real(kind=r8), dimension(:,:,:) :: x
   real(kind=r8), allocatable :: x2(:,:,:)
   integer(kind=i4) :: n,nt,nwalkin,nwalk,nwalk0,nw,iw
   real(kind=r8) :: weight
   real(kind=r8),dimension(3) :: en
   allocate(x2(3,npart,10*nwalk0))
   nt=1
   do n=1,nwalkin
      en=eloc(x(:,:,n))
      if (nwalkin.gt.1.2*nwalk0 .or. nwalkin.lt.0.8*nwalk0) then
         weight=exp(-dt*(en(1)-etrial))*nwalk0/nwalkin !weight with population control
      else
         weight=exp(-dt*(en(1)-etrial))
      endif
      nw=weight+ran(irn)
      nwalk=nwalk+nw-1
      do iw=nt,nt+nw
         x2(:,:,iw)=x(:,:,n)
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
   ave2=0.0_r8
   ave2=sum(obs**2.0_r8)/num
   err=sqrt(abs(ave2-ave**2.0_r8)/num)
   !err=sqrt(sum((obs-ave)**2.0_r8)/num)
   end subroutine stat
end module mymod
