module mymod
   use stack
   use random
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), parameter :: dx=0.0005_r8
   real(kind=r8) :: delta,dt,etrial,alpha,omega
   integer(kind=i4) :: npart
   integer(kind=i8) :: irn
contains
   subroutine setmod(deltain,dtin,alphain,omegain,etrialin,npartin)
   real(kind=r8) :: deltain,dtin,etrialin,alphain,omegain
   integer(kind=i4) :: npartin
   delta=deltain
   dt=dtin
   alpha=alphain
   omega=omegain
   etrial=etrialin
   npart=npartin
   end subroutine setmod

   function lpsit(x)
   real(kind=r8), dimension(:,:) :: x
   real(kind=r8) :: lpsit
   integer(kind=i4) :: i
   lpsit=0.0_r8
   do i=1,npart
      lpsit=lpsit-0.5_r8*alpha*sum(x(:,i)**2.0_r8)
   enddo
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
   real(kind=r8), dimension(:,:) :: x
   real(kind=r8), dimension(3) :: eloc !total, kinetic, potential energies
   integer(kind=i4) :: i
   real(kind=r8),dimension(3,npart) :: grad
   eloc=0.0_r8
   do i=1,npart
      grad(:,i)=dlpsit(x,i)
      eloc(2)=eloc(2)-0.5_r8*(ddlpsit(x,i)+dot_product(grad(:,i),grad(:,i)))
      eloc(3)=eloc(3)+0.5_r8*omega**2.0_r8*sum(x(:,i)**2.0_r8)
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
!            !boundary conditions
!            if (xnew(ic,i,n).gt.0.5_r8*l) xnew(ic,i,n)=xnew(ic,i,n)-l
!            if (xnew(ic,i,n).lt.-0.5_r8*l) xnew(ic,i,n)=xnew(ic,i,n)+l
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
         x(ic,i)=x(ic,i)+dt*grad(ic,i)+rgauss(irn)*sqrt(dt)
!         !boundary conditions
!         if (x(ic,i).gt.0.5_r8*l) x(ic,i)=x(ic,i)-l
!         if (x(ic,i).lt.-0.5_r8*l) x(ic,i)=x(ic,i)+l
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
