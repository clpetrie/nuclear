   program harmonicoscillator
   use stack
   use step
   use random
   use wavefunction
   use estimator
   use operators
   use mympi
   implicit none
   integer, parameter :: i4=selected_int_kind(9)
   integer, parameter :: i8=selected_int_kind(15)
   integer, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8) :: omega,dt,etrial
   integer(kind=i4) :: nwalk,neq,nav,nstep,nt,nwalk1
   integer(kind=i4) :: i,j,istin,istout,it
   integer(kind=i4) :: idmc,ito
   real(kind=r8) :: csi
   integer(kind=i8) :: irn
   real(kind=r8) :: alpha
   real(kind=r8) :: tau
   real(kind=r8) :: eg,egp,egerr,egave,egnow,egavee
   type(walker) :: w
   character(len=90), dimension(:), allocatable :: answer
   integer(kind=i4), allocatable :: numfig(:)
   call mpiinit0
   if (myrank().eq.0) then
      read(5,*) irn   ! number used to initialize the random number generator
      read(5,*) omega ! strenght of the harmonic potential
      read(5,*) alpha ! alpha parameter of the trial wave function, psi=exp(-alpha*x**2/2)
      read(5,*) idmc
      read(5,*) dt    ! time step
      read(5,*) etrial ! trial energy
      read(5,*) nwalk ! number of walkers
      read(5,*) neq   ! number of blocks for equilibration
      read(5,*) nav   ! number of blocks for average
      read(5,*) nstep ! number of steps per block
      write(6,'(''One dimensional harmonic oscillator'')')
      select case (idmc)
         case(1) 
            write(6,'(''Variational Monte Carlo'')')
         case(2)
            write(6,'(''Diffusion Monte Carlo'')')
         case default
            write(6,'(''Invalid idmc!!!'')')
            stop
      end select
      write(6,'(''Random number seed = '',t40,i20)') irn
      write(6,'(''Omega strenght of the potential = '',t50,f10.5)') omega
      write(6,'(''Alpha parameter of the trial wave function = '',t50,f10.5)') alpha
      write(6,'(''Time step = '',t50,f10.5)') dt
      write(6,'(''Trial energy = '',t50,f10.5)') etrial
      write(6,'(''Walkers = '',t50,i10)') nwalk
      write(6,'(''Number of blocks for equilibration = '',t50,i10)') neq
      write(6,'(''Number of blocks for statistics = '',t50,i10)') nav
      write(6,'(''Number of steps for each block = '',t50,i10)') nstep
      write(6,*) 
   endif
   call bcast(omega)
   call bcast(alpha)
   call bcast(idmc)
   call bcast(dt)
   call bcast(etrial)
   call bcast(nwalk)
   call bcast(neq)
   call bcast(nav)
   call bcast(nstep)
   call setstep(dt,etrial,idmc)
   call setpsi(alpha,omega)
   call mpi_bcast(i,size(i),mpi_integer,0,mpi_comm_world,ierror)
   istin=1
   istout=2
   select case (idmc)
      case(1)
         call setestnum(5)
         allocate(answer(5))
         call addest(1,'total energy')
         call addest(2,'kinetic energy')
         call addest(3,'potential energy')
         call addest(4,'jackson-feenberg energy')
         call addest(5,'acceptance')
      case(2)
         call setestnum(4)
         allocate(answer(4))
         call addest(1,'total energy')
         call addest(2,'kinetic energy')
         call addest(3,'potential energy')
         call addest(4,'growth energy')
   end select
   call mpiinit1(5.0_r8,idmc.eq.1)
   call create(2,2*nwalk)
   istin=1
   istout=2
   allocate(numfig(0:nproc()-1))
   do i=1,nwalk
      if (myrank().eq.0) then
         csi=ran(irn)
         w%x=1.0_r8/omega*(0.5-csi)
         call hpsi(w)
         w%irn=irn
         call push(istin,w)
      endif
      ito=mod(i-1,nproc())
      call movewalkers(1,1,0,ito,i)
   enddo
   nwalk1=numstack(istin)
   tau=0.0_r8
   it=0
   call zerest
   do i=1,nav+neq
      if (i.eq.neq+1) then
         it=0
         call zerest
         if (myrank().eq.0) write(6,'(/,''Equilibration done!'')')
      endif
      it=it+1
      do j=1,nstep
         call step1(istin,istout)
         nwalk1=numstack(istin)
         call addall(nwalk1,nt)
         call bcast(nt)
         if (idmc.eq.2) tau=tau+dt
      enddo
      call checkpop(istin)
      call calcop(istin,istout,idmc)
      call update
      if (myrank().eq.0) then
         write (6,'(/,''iteration ='',t30,i14)') it
         write (6,'(''walker number = '',t30,i14)') nt
         if (idmc.eq.1) then
            answer=resstring(tau)
            write (6,'(a90)') (answer(j),j=1,size(answer))
            answer=resstringave(tau)
            write (6,'(a90)') (answer(j),j=1,size(answer))
         else
            answer=resstring(tau)
            write (6,'(a90)') (answer(j),j=1,size(answer)-1)
            call result(4,egnow,egerr,egave,egavee)
            eg=etrial-log(egnow)/dt
            egp=etrial-log(egnow+egerr)/dt
            egerr=abs(egp-eg)
            write (6,'(''blkgrowth energy '',t30,1x,1p,e15.8,e19.10,'' +- '',e19.10)') tau,eg,egerr
            answer=resstringave(tau)
            write (6,'(a90)') (answer(j),j=1,size(answer)-1)
            eg=etrial-log(egave)/dt
            egp=etrial-log(egave+egavee)/dt
            egerr=abs(egp-eg)
            write (6,'(''avegrowth energy '',t30,1x,1p,e15.8,e19.10,'' +- '',e19.10)') tau,eg,egerr
         endif
      endif
   enddo
   call done
   if (myrank().eq.0) then
      write (6,'(''Finished!'')')
   endif
   end program harmonicoscillator
