program qmc
use random
use mymod
implicit none
integer, parameter :: i4=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(15)
integer, parameter :: r8=selected_real_kind(15,9)
integer(kind=i4) :: idmc,neq,nav,nstep,nwalk,nwalk0
real(kind=r8), allocatable :: x(:,:,:)
real(kind=r8), allocatable :: eblock(:,:),errblock(:,:),elocal(:,:)
integer(kind=i4) :: b,s,n,i,ic
real(kind=r8),dimension(3) :: energy,error
real(kind=r8) :: tau,acceptance,rho
logical :: sff
read(5,*) irn     ! random seed
read(5,*) idmc    ! 1=vmc, 2=dmc
read(5,*) neq     ! number of blocks for equilibration
read(5,*) nav     ! number of blocks for average
read(5,*) nstep   ! number of steps per block
read(5,*) nwalk   ! number of walkers
read(5,*) sff     ! start from file?
read(5,*) npart   ! number of particles
read(5,*) delta   ! controls vmc step
read(5,*) rho     ! density L=(N/rho)**(1/3)
l=(npart/rho)**(1.0_r8/3.0_r8)
read(5,*) h2o2m   ! hbar^2/(2m)
read(5,*) dt      ! time step
read(5,*) etrial  ! trial energy
nwalk0=nwalk
call setmod(delta,dt,h2o2m,l,etrial,npart)
allocate(x(3,npart,10*nwalk0))
allocate(eblock(3,nav),errblock(3,nav),elocal(3,10*nwalk))
tau=0.0_r8
eblock=0.0_r8

open(unit=6,file='data.out') !DELETE

!initialize walkers
rewind 9
if (sff) then
   read(9,'(i10,f15.8,i20)') nwalk,tau,irn
   do n=1,nwalk
      read(9,'(100e15.7)') x(:,:,n)
   enddo
else
   do n=1,nwalk
      do i=1,npart
         do ic=1,3
            x(ic,i,n)=l*(ran(irn)-0.5_r8)
         enddo
      enddo
   enddo
endif

!DELETE TOTAL
!do i=1,npart
!   do ic=1,3
!      x(ic,i,:)=-0.5_r8*l+i*(1.0_r8/(npart+1.0_r8))
!   enddo
!enddo
!energy=eloc(x(:,:,1))
!write(*,*) 'lpsit=',lpsit(x(:,:,1))
!write(*,*) 'T=',energy(2)
!write(*,*) 'V=',energy(3)

!x(1,1,1)=0.01_r8
!x(2,1,1)=-0.3_r8
!x(3,1,1)=0.8_r8
!x(1,2,1)=-0.5_r8
!x(2,2,1)=0.3_r8
!x(3,2,1)=-0.1_r8
!write(*,*) 'r start=',sqrt(sum((x(:,1,1)-x(:,2,1)-nint((x(:,1,1)-x(:,2,1))/l)*l)**2.0_r8))
!write(*,*) 'r start=',sqrt(sum((x(:,1,1)-x(:,2,1))**2.0_r8))
!write(*,*) dlpsit(x(:,:,1),1)
!write(*,*) ddlpsit(x(:,:,1),1)
!write(*,*) eloc(x(:,:,1))
!write(*,*) lpsit(x(:,:,1))
!DELETE TILL HERE

select case (idmc)
   case(1)
      write(6,'(A)') 'Variational Monte Carlo'
   case(2)
      write(6,'(A)') 'Diffusion Monte Carlo'
end select
write(6,'(''Number of blocks for equilibration = '',t50,i10)') neq
write(6,'(''Number of blocks for statistics = '',t50,i10)') nav
write(6,'(''Number of steps for each block = '',t50,i10)') nstep
write(6,'(''Number of particles = '',t50,i10)') npart
write(6,'(''Initial number of walkers = '',t50,i10)') nwalk
write(6,'(''Density = '',t50,f10.5)') rho
select case(idmc)
   case(1)
      write(6,'(''Delta for step size = '',t50,f10.5)') delta
   case(2)
      write(6,'(''Time step = '',t50,f10.5)') dt
end select
write(6,'(''Trial energy = '',t50,f10.5)') etrial

select case (idmc)
   case(1)
      do b=1,neq+nav
         if (b.eq.neq+1) write(6,*)
         if (b.eq.neq+1) write(6,'(''equilibium done!'',t30)')
         if (b.eq.neq+1) write(6,*)
         if (b.eq.neq+1) write(6,'(''***********************************************'',t50)')
         if (b.eq.neq+1) write(6,*)
         do s=1,nstep
            call metropolis(x,nwalk,acceptance)
         enddo
         if(b.gt.neq) then
            do n=1,nwalk
               elocal(:,n)=eloc(x(:,:,n))/npart !energy/particle
            enddo
            call stat(elocal(1,1:nwalk),eblock(1,b-neq),errblock(1,b-neq))
            call stat(elocal(2,1:nwalk),eblock(2,b-neq),errblock(2,b-neq))
            call stat(elocal(3,1:nwalk),eblock(3,b-neq),errblock(3,b-neq))
            write(6,'(''block ='',t25,i6)') b-neq
            write(6,'(a30,e11.5,a5,e11.5)') 'total energy/part for block = ',eblock(1,b-neq),' +/- ',errblock(1,b-neq)
            write(6,'(a32,e11.5,a5,e11.5)') 'kinetic energy/part for block = ',eblock(2,b-neq),' +/- ',errblock(2,b-neq)
            write(6,'(a34,e11.5,a5,e11.5)') 'potential energy/part for block = ',eblock(3,b-neq),' +/- ',errblock(3,b-neq)
            write(6,'(''acceptance ='',t25,e11.5)') acceptance
            write(6,*)
         endif
      enddo

   case(2)
   do b=1,nav
      write(6,*)
      write(6,'(''block ='',t30,i14)') b
      write(6,'(''walker number = '',t30,i14)') nwalk
      do s=1,nstep
         do n=1,nwalk
            call step(x(:,:,n))
         enddo
         call branch(x,nwalk,nwalk0)
         tau=tau+dt
      enddo
      do n=1,nwalk
         elocal(:,n)=eloc(x(:,:,n))/npart ! energy/particle
      enddo
      call stat(elocal(1,1:nwalk),eblock(1,b),errblock(1,b))
      call stat(elocal(2,1:nwalk),eblock(2,b),errblock(2,b))
      call stat(elocal(3,1:nwalk),eblock(3,b),errblock(3,b))
      write(6,'(''tau ='',t25,e11.5)') tau
      write(6,'(a30,e11.5,a5,e11.5)') 'total energy/part for block = ',eblock(1,b),' +/- ',errblock(1,b)
      write(6,'(a32,e11.5,a5,e11.5)') 'kinetic energy/part for block = ',eblock(2,b),' +/- ',errblock(2,b)
      write(6,'(a34,e11.5,a5,e11.5)') 'potential energy/part for block = ',eblock(3,b),' +/- ',errblock(3,b)
      write(6,'(a31,e11.5,a5,e11.5)') 'growth energy/part for block = ',egrowth(x)
   enddo
end select
call stat(eblock(1,:),energy(1),error(1))
call stat(eblock(2,:),energy(2),error(2))
call stat(eblock(3,:),energy(3),error(3))

!write walker configuration to file fort.9
rewind 9
write(9,'(i10,f15.8,i20)') nwalk,tau,irn
do n=1,nwalk
   write(9,'(100e15.7)') x(:,:,n)
enddo
close(9)

write(6,*)
write(6,'(a)') 'Finished!'
write(6,'(a13,e11.5,a5,e11.5)') 'E/particle = ',energy(1),' +/- ',error(1)
write(6,'(a13,e11.5,a5,e11.5)') 'T/particle = ',energy(2),' +/- ',error(2)
write(6,'(a13,e11.5,a5,e11.5)') 'V/particle = ',energy(3),' +/- ',error(3)
write(6,*)

close(6) !DELETE

end program qmc
