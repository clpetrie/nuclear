program sho
use random
use mymod
implicit none
integer, parameter :: i4=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(15)
integer, parameter :: r8=selected_real_kind(15,9)
integer(kind=i4) :: idmc,neq,nav,nstep,nwalk,nwalk0
real(kind=r8), allocatable :: x(:,:)
real(kind=r8), allocatable :: eblock(:),errblock(:),elocal(:)
integer(kind=i4) :: b,s,n,i
real(kind=r8) :: energy,error,tau,acceptance
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
read(5,*) alpha   ! alpha parameter of the trial wave function, psi=exp(-alpha*x**2/2)
read(5,*) omega   ! strenght of the harmonic potential
read(5,*) dt      ! time step
read(5,*) etrial  ! trial energy
nwalk0=nwalk
call setmod(delta,alpha,omega,dt,etrial,npart)
allocate(x(10*nwalk0,npart))
allocate(eblock(nav),errblock(nav),elocal(10*nwalk))
tau=0

!initialize walkers
rewind 9
if (sff) then
   read(9,'(i10,f15.8,i20)') nwalk,tau,irn
   do n=1,nwalk
      read(9,'(100e15.7)') x(n,:)
   enddo
else
   do n=1,nwalk
      do i=1,npart
         x(n,i)=ran(irn)*2-1
      enddo
   enddo
endif

select case (idmc)
   case(1)
      write(6,'(A)') 'Variational Monte Carlo'
   case(2)
      write(6,'(A)') 'Diffusion Monte Carlo'
end select
write(6,'(''Number of blocks for equilibration = '',t50,i10)') neq
write(6,'(''Number of blocks for statistics = '',t50,i10)') nav
write(6,'(''Number of steps for each block = '',t50,i10)') nstep
write(6,'(''Walkers = '',t50,i10)') nwalk
write(6,'(''Alpha parameter of the trial wave function = '',t50,f10.5)') alpha
write(6,'(''Omega strength of the potential = '',t50,f10.5)') omega
write(6,'(''Time step = '',t50,f10.5)') dt
write(6,'(''Trial energy = '',t50,f10.5)') etrial

eblock=0
select case (idmc)
   case(1)
      do b=1,neq+nav
         if (b.eq.neq) write(6,*)
         if (b.eq.neq) write(6,'(''equilibium done!'',t30)')
         if (b.eq.neq) write(6,*)
         if (b.eq.neq) write(6,'(''***********************************************'',t50)')
         if (b.eq.neq) write(6,*)
        do s=1,nstep
            call metropolis(x,nwalk,acceptance)
         enddo
         if(b.gt.neq) then
            do n=1,nwalk
               elocal(n)=eloc(x(n,:))
            enddo
            call stat(elocal(1:nwalk),eblock(b-neq),errblock(b-neq))
            write(6,'(''block ='',t25,i6)') b-neq
            write(6,'(a19,e11.5,a5,e11.5)') 'energy for block = ',eblock(b-neq),' +/- ',errblock(b-neq)
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
            call step(x(n,:))
         enddo
         call branch(x,nwalk,nwalk0)
         tau=tau+dt
      enddo
      do n=1,nwalk
         elocal(n)=eloc(x(n,:))
      enddo
      call stat(elocal(1:nwalk),eblock(b),errblock(b))
      write(6,'(''tau ='',t25,e11.5)') tau
      write(6,'(a19,e11.5,a5,e11.5)') 'energy for block = ',eblock(b),' +/- ',errblock(b)
      write(6,'(a26,e11.5,a5,e11.5)') 'growth energy for block = ',egrowth(x,nwalk)
      write(6,*)
   enddo
end select
call stat(eblock,energy,error)

!write walker configuration to file fort.9
rewind 9
write(9,'(i10,f15.8,i20)') nwalk,tau,irn
do n=1,nwalk
   write(9,'(100e15.7)') x(n,:)
enddo
close(9)

write(6,*)
write(6,'(a)') 'Finished!'
write(6,'(a4,e11.5,a5,e11.5)') 'E = ',energy,' +/- ',error
write(6,*)

end program sho
