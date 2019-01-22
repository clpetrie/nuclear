module random
    implicit none
    integer, private, parameter :: i4=selected_int_kind(9)
    integer, private, parameter :: i8=selected_int_kind(15)
    integer, private, parameter :: r8=selected_real_kind(15,9)
contains
    function ran(irn)
    integer(kind=i8),  parameter :: mask24 = ishft(1_i8,24)-1
    integer(kind=i8),  parameter :: mask48 = ishft(1_i8,48_i8)-1_i8
    real(kind=r8),  parameter :: twom48=2.0_r8**(-48)
    integer(kind=i8),  parameter :: mult1 = 44485709377909_i8
    integer(kind=i8),  parameter :: m11 = iand(mult1,mask24)
    integer(kind=i8),  parameter :: m12 = iand(ishft(mult1,-24),mask24)
    integer(kind=i8),  parameter :: iadd1 = 96309754297_i8
    integer(kind=i8) :: irn
    real(kind=r8) :: ran
    integer(kind=i8) :: is1,is2
    is2=iand(ishft(irn,-24),mask24)
    is1=iand(irn,mask24)
    irn=iand(ishft(iand(is1*m12+is2*m11,mask24),24)+is1*m11+iadd1,mask48)
    ran=ior(irn,1_i8)*twom48
    return
    end function ran

   function ran2(irn)
   integer(kind=i8),  parameter :: mask24 = ishft(1_i8,24)-1
   integer(kind=i8), parameter :: mask48 = ishft(1_i8,48)-1
   real(kind=r8), parameter :: twom48=2.0_r8**(-48)
   integer(kind=i8), parameter :: mult2 = 34522712143931_i8
   integer(kind=i8), parameter :: m21 = iand(mult2,mask24)
   integer(kind=i8), parameter :: m22 = iand(ishft(mult2,-24),mask24)
   integer(kind=i8), parameter :: iadd2 = 55789347517_i8
   integer(kind=i8) :: irn
   real(kind=r8) :: ran2
   integer(kind=i8) :: is1,is2
   is2 = iand(ishft(irn,-24),mask24)
   is1 = iand(irn,mask24)
   irn = iand(ishft(iand(is1*m22+is2*m21,mask24),24)+is1*m21+iadd2,mask48)
   ran2 = ior(irn,1_i8)*twom48
   return
   end function ran2

    function rgauss(irn)
    real(kind=r8) :: pi
    real(kind=r8) :: rgauss,x1,x2
    integer(kind=i8) :: irn
    x1=ran(irn)
    x2=ran(irn)
    pi=4.0_r8*atan(1.0_r8)
    rgauss=sqrt(-2.0_r8*log(x1))*cos(2.0_r8*pi*x2)
    return
    end function rgauss
end module random
