module step
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, save ::  dt,sigma,etrial
   integer(kind=i4), private, save :: idmc
contains
   subroutine setstep(dtin,etrialin,idmcin)
   real(kind=r8) :: dtin,etrialin
   integer(kind=i4) :: idmcin
   dt=dtin
   etrial=etrialin
   idmc=idmcin
   sigma=sqrt(dt)
   end subroutine setstep

   subroutine step1(istin,istout)
   use stack
   use random
   use wavefunction
   use estimator
   integer(kind=i4) :: istin,istout
   logical :: empty
   real(kind=r8) :: csi,prob,rn,eloc
   integer(kind=i4) :: iwt,i
   type(walker) :: w1,w2
   do while (.true.)
      call pop(istin,w1,empty)
      if (empty) exit
      select case (idmc)
         case (1) ! VMC step
            csi=ran(w1%irn)
            w2%x=w1%x+sigma*(csi-0.5_r8)
            call hpsi(w2)
            prob=(w2%psi/w1%psi)**2
            rn=ran(w1%irn)
            if (rn.lt.prob) then ! the move is accepted!
               call addval(5,1.0_r8,1.0_r8)
            else ! the move has been rejected
               w2=w1
               call addval(5,0.0_r8,1.0_r8)
            endif
            w2%irn=w1%irn
            w2%weight=1.0_r8
            call push(istout,w2)
         case(2) ! DMC step
            w2%irn=w1%irn
            w2%x=w1%x+sigma**2*w1%dpsi+sigma*rgauss(w2%irn)
            call hpsi(w2)
            eloc=-0.5_r8*w2%d2psi+w2%v
            w2%weight=exp(-dt*(eloc-etrial))
            call addval(4,w2%weight,1.0_r8)
            iwt=w2%weight+ran(w2%irn)
            w2%weight=1.0_r8
            do i=1,iwt
               csi=ran2(w2%irn)
               call push(istout,w2)
            enddo
      end select
   enddo
   istout=istin
   istin=3-istin
   end subroutine step1
end module step
