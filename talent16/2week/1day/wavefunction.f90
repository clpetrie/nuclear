module wavefunction
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, save :: alpha,omega
contains
   subroutine setpsi(alphain,omegain)
   real(kind=r8) :: alphain,omegain
   alpha=alphain
   omega=omegain
   end subroutine setpsi

   subroutine hpsi(w)
   use stack
   type(walker) :: w
   w%psi=exp(-0.5_r8*alpha*w%x**2)
   w%dpsi=-alpha*w%x
   w%d2psi=alpha**2*w%x**2-alpha
   w%v=0.5_r8*omega**2*w%x**2
   end subroutine hpsi
end module wavefunction
