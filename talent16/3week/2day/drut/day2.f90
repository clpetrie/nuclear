program day2
use random
use step
implicit none
integer, parameter :: i4=selected_int_kind(9)
integer, parameter :: i8=selected_int_kind(15)
integer, parameter :: r8=selected_real_kind(15,9)
integer(kind=i4), parameter :: L=20_i4
integer(kind=i4), parameter :: Lt=20_i4
integer(kind=i4), parameter :: neq=10_i4
integer(kind=i4), parameter :: nav=1000_i4
integer(kind=i4), parameter :: nsteps=10_i4
real(kind=r8), parameter :: g=1.0_r8
real(kind=r8), parameter :: tau=0.05_r8
real(kind=r8), parameter :: delta=3.0_r8 !step size for metropolis
real(kind=r8), parameter :: beta=tau*Lt
real(kind=r8), parameter :: mu=0.0_r8
real(kind=r8), parameter :: c=exp(g*tau)-1.0_r8
real(kind=r8), allocatable :: s(:,:)
integer(kind=i4) :: nx,nt,n,ns,niter
real(kind=r8), dimension(0:L-1,0:L-1) :: phi
real(kind=r8), dimension(0:L-1,0:L-1) :: phiin
real(kind=r8) :: ldet,accept,attempt
irn=17
allocate(s(0:L-1,0:Lt-1))
call setstep(L,Lt,c,tau,delta,beta,mu)

!initialize sigma
do nx=0,L-1
   do nt=0,Lt-1
      s(nx,nt)=ran(irn)-0.5_r8
   enddo
enddo

!initialize s.p. waves
do n=0,L-1
   do nx=0,L-1
      phi(n,nx)=0.0_r8
      if (nx.eq.n) phi(n,nx)=1.0_r8
   enddo
enddo

phiin=phi
call calcdet(s,phiin,ldet) !get the log det of the initial config
niter=0_i4
accept=0.0_r8
attempt=0.0_r8
do ns=1,neq+nav
   niter=niter+1_i4
   if (ns.eq.neq+1) then
      write(*,*)
      write(*,*) 'acceptance=',accept/attempt
      write(*,*) 'equilibrium done!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      niter=1_i4
      accept=0.0_r8
      attempt=0.0_r8
   endif
   call metropolis(s,phi,ldet,accept,attempt)
   if (ns.gt.neq .and. mod(ns,nsteps).eq.0) then
      write(*,*)
      write(*,*) 'iteration=',niter
      write(*,*) 'acceptance=',accept/attempt
      call calcdet(s,phi,ldet)
      write(*,*) 'log(det(M))=',ldet
      accept=0.0_r8
      attempt=0.0_r8
   endif
enddo

end program day2
