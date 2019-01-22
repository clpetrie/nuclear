module step
   use random
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i8) :: irn
   integer(kind=i4), private :: L,Lt
   real(kind=r8), private :: c,tau,delta,beta,mu,z
contains
   subroutine setstep(Lin,Ltin,cin,tauin,deltain,betain,muin)
   integer(kind=i4) :: Lin,Ltin
   real(kind=r8) :: cin,tauin,deltain,betain,muin
      L=Lin
      Lt=Ltin
      c=cin
      tau=tauin
      delta=deltain
      beta=betain
      mu=muin
      z=exp(beta*mu)
   end subroutine setstep

   subroutine metropolis(s,phiin,ldetold,accept,attempt)
   integer(kind=i4) :: nx,nt,n
   real(kind=r8), dimension(0:L-1,0:Lt-1) :: s
   real(kind=r8), dimension(0:L-1,0:Lt-1) :: phi
   real(kind=r8), dimension(0:L-1,0:Lt-1) :: phiin
   real(kind=r8), allocatable :: mmat(:,:)
   real(kind=r8), allocatable :: snew(:,:)
   real(kind=r8) :: ut,ldetold,ldetnew,accept,attempt
   allocate(snew(0:L-1,0:Lt-1))
   allocate(mmat(0:L-1,0:L-1))
   phi=phiin
   !propose move
   snew=s
   do nx=0,L-1
      phi=phiin
      do nt=0,Lt-1
         snew(nx,nt)=s(nx,nt)+delta*(ran(irn)-0.5_r8)
         call calcdet(snew,phi,ldetnew)
         attempt=attempt+1.0_r8
         if (ran(irn).lt.exp(ldetnew-ldetold)) then
            accept=accept+1.0_r8
            s(nx,nt)=snew(nx,nt)
            ldetold=ldetnew
         endif
      enddo
   enddo
   end subroutine metropolis

   subroutine calcdet(s,phi,ldet)
   real(kind=r8), dimension(0:L-1,0:Lt-1) :: s
   real(kind=r8), dimension(0:L-1,0:L-1) :: phi
   real(kind=r8), dimension(0:L-1,0:L-1) :: M
   real(kind=r8), dimension(0:L-1,0:L-1) :: Mdiag
   integer(kind=i4) :: n,nx,nt,i
   real(kind=r8) :: ldet
   integer :: info
   integer, allocatable :: ipiv(:)
   allocate(ipiv(0:L-1))
   do n=0,L-1
      do nx=0,L-1
         M(n,nx)=0.0_r8
         if(n.eq.nx) M(n,nx)=1.0_r8 !identity matrix
         do nt=0,Lt-1
!            phi(n,nx)=0.25_r8*tau*(phi(n,mod(nx+1,L))+phi(n,mod(nx-1+L,L)))+(1.0_r8-0.5_r8*tau)*phi(n,nx)
            phi(n,nx)=(1.0_r8+0.25_r8*tau*(phi(n,mod(nx+1,L))+phi(n,mod(nx-1+L,L))-2.0_r8*phi(n,nx)))*phi(n,nx)
            phi(n,nx)=(1.0_r8+sqrt(c)*s(nx,nt))*phi(nx,nt)
            phi(n,nx)=(1.0_r8+0.25_r8*tau*(phi(n,mod(nx+1,L))+phi(n,mod(nx-1+L,L))-2.0_r8*phi(n,nx)))*phi(n,nx)
!            phi(n,nx)=0.25_r8*tau*(phi(n,mod(nx+1,L))+phi(n,mod(nx-1+L,L)))+(1.0_r8-0.5_r8*tau)*phi(n,nx)
         enddo
         M(n,nx)=M(n,nx)+z*phi(n,nx)
      enddo
   enddo
   Mdiag=M
   call dgetrf(L,L,Mdiag,L,ipiv,info)
   ldet=0.0_r8
   do i=0,L-1 
   !   if (ipiv(i).ne.i) then 
      if (Mdiag(i,i).lt.0) then 
         ldet=ldet+log(-Mdiag(i,i)) 
      else 
         ldet=ldet+log(Mdiag(i,i)) 
      endif 
   enddo
   end subroutine calcdet
end module step
