module estimator
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   character(len=30), private, dimension(:), allocatable, save :: label
   real(kind=r8), private, dimension(:), allocatable, save :: valblk,valnow,wtblk
   real(kind=r8), private, dimension(:), allocatable, save :: valblk2,errnow
   real(kind=r8), private, dimension(:), allocatable, save :: valsum,valsum2
   integer(kind=i4), private, save :: nest
   integer(kind=i4), private, save, allocatable :: nval(:),nvalsum(:)
contains
   subroutine setestnum(n)
   integer(kind=i4) :: n
   nest=n
   allocate(valblk(n),valnow(n))
   allocate(valblk2(n),errnow(n))
   allocate(valsum(n),valsum2(n))
   allocate(wtblk(n))
   allocate(label(n))
   allocate(nval(n),nvalsum(n))
   end subroutine setestnum

   subroutine zerest
   integer(kind=i4) :: i
   do i=1,nest
      valblk(i)=0.0_r8
      valblk2(i)=0.0_r8
      valnow(i)=0.0_r8
      errnow(i)=0.0_r8
      valsum(i)=0.0_r8
      valsum2(i)=0.0_r8
      wtblk(i)=0.0_r8
      nval(i)=0
      nvalsum(i)=0
   enddo
   end subroutine zerest

   subroutine addest(i,l)
   integer(kind=i4) :: i
   character(len=*) :: l
   integer(kind=i4) :: ln,j
   ln=min(len(l),len(label(i)))
   label(i)(1:ln)=l(1:ln)
   do j=ln+1,len(label(i))
      label(i)(j:j)=" "
   enddo
   return
   end subroutine addest

   subroutine addval(i,val,wt)
   integer(kind=i4) :: i
   real(kind=r8) :: val,wt
   valblk(i)=valblk(i)+val*wt
   valblk2(i)=valblk2(i)+val*val*wt
   wtblk(i)=wtblk(i)+wt
   nval(i)=nval(i)+1
   return
   end subroutine addval

   subroutine update
   use mympi
   integer(kind=i4) :: i
   real(kind=r8) :: valblkall(nest),valblk2all(nest),wtblkall(nest)
   integer(kind=i4) :: nvalall(nest)
   call addall(valblk,valblkall)
   call addall(valblk2,valblk2all)
   call addall(wtblk,wtblkall)
   call addall(nval,nvalall)
   if (myrank().eq.0) then
      do i=1,nest
         valnow(i)=valblkall(i)/wtblkall(i)
         valblk2(i)=valblk2all(i)/wtblkall(i)
         errnow(i)=sqrt(abs(valblk2(i)-valnow(i)**2)/nvalall(i))
         valsum(i)=valsum(i)+valnow(i)
         valsum2(i)=valsum2(i)+valnow(i)**2
         nvalsum(i)=nvalsum(i)+1
      enddo
   endif
   valblk=0.0_r8
   valblk2=0.0_r8
   wtblk=0.0_r8
   nval=0
   end subroutine update
 
   subroutine result(i,val,err,vals,vale)
   integer(kind=i4) :: i
   real(kind=r8) :: val,err,vals,vale
   val=valnow(i)
   err=errnow(i)
   vals=valsum(i)/nvalsum(i)
   vale=sqrt(abs(valsum2(i)/nvalsum(i)-vals**2)/nvalsum(i))
   return
   end subroutine result

   function resstring(tau)
   character(len=90), dimension(nest) :: resstring
   integer(kind=i4) :: i
   real(kind=r8) :: val,err,vals,vale
   real(kind=r8) :: tau
   do i=1,nest
      call result(i,val,err,vals,vale)
      write (resstring(i),'(a30,t30,1x,1p,e15.8,e19.10,'' +- '',e19.10)') 'blk'//label(i),tau,val,err
   enddo
   return
   end function resstring

   function resstringave(tau)
   character(len=90), dimension(nest) :: resstringave
   integer(kind=i4) :: i
   real(kind=r8) :: val,err,vals,vale
   real(kind=r8) :: tau
   do i=1,nest
      call result(i,val,err,vals,vale)
      write (resstringave(i),'(a30,t30,1x,1p,e15.8,e19.10,'' +- '',e19.10)') 'ave'//label(i),tau,vals,vale
   enddo
   return
   end function resstringave
end module estimator
