module step
   use random
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i8) :: irn
   integer(kind=i4), private :: L,Lt
   real(kind=r8), private :: m,a,at
contains
   subroutine setstep(Lin,Ltin,m_in,ain,atin)
   integer(kind=i4) :: Lin,Ltin
   real(kind=r8) :: m_in,ain,atin
      L=Lin
      Lt=Ltin
      m=m_in
      a=ain
      at=atin
   end subroutine setstep
end module step
