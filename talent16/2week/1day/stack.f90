module stack   
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   type :: walker
      real(kind=r8) :: x,psi,dpsi,d2psi
      real(kind=r8) :: v,weight
      integer(kind=i8) :: irn
   end type
   type (walker), private, allocatable, save :: s(:,:)
   integer(kind=i4), private, allocatable, save :: ist(:) 
   integer(kind=i4), private, save :: nstack,nwalk
   interface assignment (=)
      module procedure copywalker
   end interface
contains
   subroutine copywalker(wl,wr)
   type (walker), intent(inout) :: wl
   type (walker), intent(in) :: wr
   wl%x=wr%x
   wl%psi=wr%psi
   wl%dpsi=wr%dpsi
   wl%d2psi=wr%d2psi
   wl%v=wr%v
   wl%weight=wr%weight
   wl%irn=wr%irn
   end subroutine copywalker

   subroutine create(nst,nw)
   integer(kind=i4) :: nst,nw
   nstack=nst
   nwalk=nw
   allocate(s(nwalk,nstack),ist(nstack))
   ist=0
   end subroutine create

   subroutine push(i,w)
   integer(kind=i4) :: i
   type(walker) :: w
   if (ist(i).ge.nwalk) then
      write (6,'(1x,'' stack overflow'',i10)') ist(i)
      stop
   endif
   ist(i)=ist(i)+1
   s(ist(i),i)=w
   end subroutine push

   subroutine pop(i,w,empty)
   integer(kind=i4) :: i
   type(walker) :: w
   logical :: empty
   empty=ist(i).eq.0
   if (.not.empty) then
      w=s(ist(i),i)
      ist(i)=ist(i)-1
   endif
   end subroutine pop

   function numstack(i)
   integer(kind=i4) :: i,numstack
   numstack=ist(i)
   end function numstack
end module stack
