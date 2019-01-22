program block
implicit none
integer, parameter :: sbin=3
integer :: num, inum, eof, iter, i, bin, v
real, dimension(10000) :: rin, np0in, np0errin, np1in, np1errin, ppin, pperrin, nnin, nnerrin
real, dimension(9,10000) :: inv
real, allocatable :: outv(:,:)
real :: norm, dr
!real, allocatable :: ro(:), np0o(:), np0erro(:), np1o(:), np1erro(:), ppo(:), pperro(:), nno(:), nnerro(:)

open(8,file='gnp.dmc')
iter=0
do
!   read(8,*,iostat=eof) rin(iter+1), np0in(iter+1), np0errin(iter+1), np1in(iter+1), np1errin(iter+1), ppin(iter+1), pperrin(iter+1), nnin(iter+1), nnerrin(iter+1)
   read(8,*,iostat=eof) inv(1,iter+1), inv(2,iter+1), inv(3,iter+1), inv(4,iter+1), inv(5,iter+1), inv(6,iter+1), &
      inv(7,iter+1), inv(8,iter+1), inv(9,iter+1)
   if(eof.lt.0) exit
   iter=iter+1
enddo
close(8)
iter=iter-1 !remove the empty one
num=iter/sbin
allocate(outv(9,num))

outv=0 !I get weird NaN's if I don't initialize it here.
!bin data
do v=1,9
   if (v.eq.3 .or. v.eq.5 .or. v.eq.7 .or. v.eq.9) cycle
   bin=0
   inum=1
   do i=1,iter
      bin=bin+1
      if (bin.gt.sbin) then
         outv(v,inum)=outv(v,inum)/sbin
         inum=inum+1
         if(inum.gt.num) cycle
         outv(v,inum)=outv(v,inum)+inv(v,i)
         bin=1
      else
         outv(v,inum)=outv(v,inum)+inv(v,i)
      endif
   enddo
enddo

!bin errors
do v=1,9
   if (v.eq.1 .or. v.eq.2 .or. v.eq.4 .or. v.eq.6 .or. v.eq.8) cycle
   bin=0
   inum=1
   do i=1,iter
      bin=bin+1
      if (bin.gt.sbin) then
         outv(v,inum)=sqrt(outv(v,inum)/sbin)
         inum=inum+1
         if(inum.gt.num) cycle
         outv(v,inum)=outv(v,inum)+inv(v,i)**2
         bin=1
      else
         outv(v,inum)=outv(v,inum)+inv(v,i)**2
      endif
   enddo
enddo


!normalize to 1
dr=outv(1,2)-outv(1,1)
do v=1,9
   if (v.eq.1 .or. v.eq.3 .or. v.eq.5 .or. v.eq.7 .or. v.eq.9) cycle
   norm=0
   do i=1,num
      norm=norm+outv(v,i)
   enddo
   norm=1.0/(norm*dr)
   outv(v,:)=outv(v,:)*norm
   outv(v+1,:)=outv(v+1,:)*norm
enddo

open(10,file='gnp.dmc.block')
do i=1,num
   write(10,*) outv(:,i)
enddo
close(10)

end program block
