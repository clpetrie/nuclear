module operators
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
contains
   subroutine calcop(istin,istout,idmc)
   use stack
   use estimator
   use wavefunction
   integer(kind=i4) :: istin,istout,idmc
   real(kind=r8) :: ek,ev,en,ejf
   logical :: empty
   type(walker) :: w
   do while (.true.)
      call pop(istin,w,empty)
      if (empty) exit
      call hpsi(w)
      ek=-0.5_r8*w%d2psi
      ev=w%v
      en=ek+ev
      call addval(1,en,w%weight)
      call addval(2,ek,w%weight)
      call addval(3,ev,w%weight)
      if (idmc.eq.1) then
         ejf=-0.25_r8*(w%d2psi-w%dpsi*w%dpsi)
         call addval(4,ejf,w%weight)
      endif
      call push(istout,w)
   enddo
   istout=istin
   istin=3-istin
   end subroutine calcop
end module operators
