   program crunch
   implicit none
   integer, parameter :: i4=selected_int_kind(9)
   integer, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4) :: nobs,ndata,n,n1,n2
   integer(kind=i4), parameter :: ndatamax=1000000
   character(len=25), allocatable :: label(:),label2(:,:)
   character(len=80) :: word,readw,filename,filename1,filename2,filename3
   character(len=90) :: line
   character(len=500) :: longline
   character(len=2) :: cc
   character(len=17) :: cdummy
   real(kind=r8), allocatable :: obstmp(:,:),errtmp(:,:),imtimetmp(:)
   real(kind=r8), allocatable :: obs(:,:),err(:,:),imtime(:)
   real(kind=r8), allocatable :: aveobs(:,:),errobs(:,:),errobsblk(:,:),errblk(:)
   real(kind=r8), allocatable :: corr(:,:),weight(:),r(:)
   real(kind=r8), allocatable :: rho(:,:,:),rhoerr(:,:,:)
   real(kind=r8), allocatable :: averho(:,:),errrho(:,:)
   integer(kind=i4) :: nline,i,idx,imin,step,nprnt,nblk,ic,idummy,iobs(2),ivmc,idmc
   integer(kind=i4) :: ntab,j,k,obsidx,ndat1,ndat2
   real(kind=r8) :: ref,dummy,extobs,exterr
   real(kind=r8), allocatable :: extrho(:,:),extrhoerr(:,:)
   logical :: domix,domix2,mcfile,endofline,endoftable,endoffile
   domix=.false.
   domix2=.false.
   if (command_argument_count().eq.1) then
      call get_command_argument(1,filename)
      filename=adjustl(filename)
      filename=filename(1:index(filename,' ')-1)
   else if (command_argument_count().eq.2) then
      call get_command_argument(1,filename1)
      filename1=adjustl(filename1)
      filename1=filename1(1:index(filename1,' ')-1)
      call get_command_argument(2,filename2)
      filename2=adjustl(filename2)
      filename2=filename2(1:index(filename2,' ')-1)
      domix=.true.
   else if (command_argument_count().eq.3) then
      call get_command_argument(1,filename1)
      filename1=adjustl(filename1)
      filename1=filename1(1:index(filename1,' ')-1)
      call get_command_argument(2,filename2)
      filename2=adjustl(filename2)
      filename2=filename2(1:index(filename2,' ')-1)
      call get_command_argument(3,filename3)
      filename3=adjustl(filename3)
      filename3=filename3(1:index(filename3,' ')-1)
      domix2=.true.
   else 
      write(6,'(''Error: you should use "exec filename" or "exec filevmc filedmc"'')')
      stop
   endif
   if (.not.domix.and..not.domix2) then
      open(unit=10,file=filename)
      write(6,'(''Crunch program''/)')
      write(6,'(''Name of the file to be analized: '',a80)') filename
! try to understand if this is a vmc/dmc output or some other observable
      mcfile=.false.
      call findw('Variational',nline)
      if (nline.ne.-1) then 
         write(6,'(''Found VMC output file'')')
         mcfile=.true.
      endif
      rewind 10
      if (nline.eq.-1) call findw('Diffusion',nline)
      if (nline.ne.-1.and..not.mcfile) then 
         write(6,'(''Found DMC output file'')')
         mcfile=.true.
      endif
      rewind 10
      if (nline.eq.-1) call findw('Calculate',nline)
      if (nline.ne.-1.and..not.mcfile) then 
         write(6,'(''Found output file analizing existing configurations'')')
         mcfile=.true.
      endif
      rewind 10
      if (nline.eq.-1) then
         write(6,'(''This seems to be a file with observables'')')
      endif
! first try to estimate the number of observables and their name
      if (mcfile) then
         nline=1
         call findw('equilibration done',nline)
         call findw('walker number',nline)
         nobs=1
         call findw('Time for',nobs)
         nobs=nobs-2 ! number of observables
!        nobs=nobs-2 ! skip the last two values, complex energy and complex weight
         write(6,'(''Number of observables = '',i5)') nobs
         allocate(label(nobs))
         allocate(obstmp(nobs,ndatamax),errtmp(nobs,ndatamax),imtimetmp(ndatamax))
         rewind 10  ! restart from the beginning of the file
         call findw('equilibration done',nline)
         call findw('walker number',nline)
         write(6,'(/''Observables read from input:'')')
         do i=1,nobs
            read(10,'(a90)') line
            label(i)=line(1:25)
         enddo
         call listobs(nobs,label)
         rewind 10  ! restart from the beginning of the file
         call findw('equilibration done',nline)
         word='Never!'
         readw=''
         idx=1
         do while (readw.ne.word)
            call findw('walker number',nline)
            if (nline.eq.-1) exit
            do i=1,nobs
               read(10,'(a90)') line
               line=line(index(line,'  '):len(line))
               line=adjustl(line)
               read(line,*) imtimetmp(idx),obstmp(i,idx),cc,errtmp(i,idx)
            enddo
            if (idx.ge.2) then ! vmc case, where imtime is always zero
               if (imtimetmp(idx).eq.imtimetmp(idx-1)) then
                  imtimetmp(idx-1)=idx-1
                  imtimetmp(idx)=idx
               endif
            endif
            idx=idx+1
         enddo
         ndata=idx-1
         close(10)
         filename=' '
         write(6,'(/''Number of data = '',i10)') ndata
         allocate(obs(nobs,ndata),err(nobs,ndata),imtime(ndata))
         obs(:,:)=obstmp(:,1:ndata)
         err(:,:)=errtmp(:,1:ndata)
         imtime(:)=imtimetmp(1:ndata)
         deallocate(obstmp,errtmp,imtimetmp)
         allocate(aveobs(nobs,ndata),errobs(nobs,ndata),errobsblk(nobs,ndata))
         allocate(errblk(nobs))
         allocate(corr(nobs,ndata))
         imin=1
         call avecalc(obs,imin,ndata,aveobs,errobs,nobs,1,errblk)
         write(6,'(''Average of data calculated!'')')
         write(6,'(''Setup done!'')')
         write(6,'(''Type help if you dont know what to do'')')
         do while(.true.)
            write(6,'(''---> '',$)')
            read(5,'(a90)') line
            if (line(1:3).eq.'ave') then
               line=line(index(line,' '):len(line))
               line=adjustl(line)
               read(line,*,err=10,end=10) nblk
               call avecalc(obs,imin,ndata,aveobs,errobs,nobs,nblk,errblk)
               call showobs(nobs,label,aveobs(:,ndata),errobs(:,ndata),errblk,filename)
               cycle
               10 write(6,'(''Use: "ave nblk"'')')
            else if (line(1:3).eq.'obs') then
               call listobs(nobs,label)
            else if (line(1:4).eq.'corr') then
               line=line(index(line,' '):len(line))
               line=adjustl(line)
               read(line,*,err=20,end=20) ref
               call autocorr(nobs,obs,imin,ndata,ref,corr,label)
               cycle
               20 write(6,'(''Use: "corr limit"'')')
            else if (line(1:5).eq.'print') then
               line=line(index(line,' '):len(line))
               line=adjustl(line)
               if (line(1:3).eq.'ave') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=30,end=30) step
                  call prnt(filename,nobs,label,aveobs,errobs,imin,ndata,step,imtime)
                  cycle
                  30 write(6,'(''Use: "print ave step"'')')
               else if (line(1:4).eq.'corr') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=40,end=40) nprnt
                  call prntcorr(filename,corr,nprnt,ndata)
                  cycle
                  40 write(6,'(''Use: "print corr observable"'')')
               else if (line(1:3).eq.'err') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=50,end=50) n
                  call prnterr(filename,errobsblk,n,ndata)
                  cycle
                  50 write(6,'(''Use: "print err observable"'')')
               else
                  write(6,'(''Use: "print ave/corr/err"'')')
               endif
            else if (line(1:4).eq.'list') then
               call showobs(nobs,label,aveobs(:,ndata),errobs(:,ndata),errblk,filename)
               write(6,'(/''Note: calculate averages before list if you have not done yet!'')')
               write(6,'(''List written to file '',a)') filename
            else if (line(1:6).eq.'errors') then
               line=line(index(line,' '):len(line))
               line=adjustl(line)
               call errors(obs,imin,ndata,errobsblk,nobs)
               cycle
               60 write(6,'(''Use: "errors"'')')
            else if (line(1:4).eq.'plot') then
               line=line(index(line,' '):len(line))
               line=adjustl(line)
               if (line(1:7).eq.'obs ave'.or.line(1:7).eq.'ave obs') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=70,end=70) n
                  call plot2xydy(imtime(imin:ndata),obs(n,imin:ndata),err(n,imin:ndata), &
                         & aveobs(n,imin:ndata),errobs(n,imin:ndata))
                  cycle
                  70 write(6,'(''Use: "plot obs ave observable"'')')
               else if (line(1:7).eq.'obs obs') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=112,end=112) n1,n2
                  call plot2xydy(imtime(imin:ndata),obs(n1,imin:ndata),err(n1,imin:ndata), &
                         & obs(n2,imin:ndata),err(n2,imin:ndata))
                  cycle
                  112 write(6,'(''Use: "plot obs obs observable1 observable2"'')')
               else if (line(1:7).eq.'ave ave') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=114,end=114) n1,n2
                  call plot2xydy(imtime(imin:ndata),aveobs(n1,imin:ndata),errobs(n1,imin:ndata), &
                         & aveobs(n2,imin:ndata),errobs(n2,imin:ndata))
                  cycle
                  114 write(6,'(''Use: "plot ave ave observable1 observable2"'')')
               else if (line(1:3).eq.'obs') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=80,end=80) n
                  call plotxydy(imtime(imin:ndata),obs(n,imin:ndata),err(n,imin:ndata))
                  cycle
                  80 write(6,'(''Use: "plot obs observable"'')')
               else if (line(1:3).eq.'ave') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=90,end=90) n
                  call plotxydy(imtime(imin:ndata),aveobs(n,imin:ndata),errobs(n,imin:ndata))
                  cycle
                  90 write(6,'(''Use: "plot ave observable"'')')
               else if (line(1:4).eq.'corr') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=100,end=100) n
                  call ploty(corr(n,:))
                  cycle
                  100 write(6,'(''Use: "plot corr observable"'')')
               else if (line(1:6).eq.'errors') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=110,end=110) n
                  call ploty(errobsblk(n,:))
                  cycle
                  110 write(6,'(''Use: "plot errors observable"'')')
               else
                  write(6,'(''Use: "plot ave/obs/errors"'')')
               endif
            else if (line(1:7).eq.'analize') then
               line=line(index(line,' '):len(line))
               line=adjustl(line)
               read(line,*,err=120,end=120) n
               call avecalc(obs,imin,ndata,aveobs,errobs,nobs,1,errblk)
               call autocorr(nobs,obs,imin,ndata,0.2_r8,corr,label)
               call errors(obs,imin,ndata,errobsblk,nobs)
               call plotall(imtime(imin:ndata),obs(n,imin:ndata),err(n,imin:ndata), &
                 &   aveobs(n,imin:ndata),errobs(n,imin:ndata),corr(n,:), &
                 &   errobsblk(n,:))
               cycle
               120 write(6,'(''Use: "analize obs"'')')
            else if (line(1:3).eq.'set') then
               line=line(index(line,' '):len(line))
               line=adjustl(line)
               if (line(1:4).eq.'imin') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=130,end=130) n
                  imin=max(1,n)
                  if (imin.gt.ndata) then
                     write(6,'(''imin cannot be larger than ndata!'')')
                     imin=ndata
                  endif
                  write(6,'(''imin set to '',i10)') imin
                  write(6,'(''Number of data '',i10)') ndata-imin+1
                  cycle
                  130 write(6,'(''Use: set imin val"'')')
               else if (line(1:8).eq.'filename') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=140,end=140) filename
                  write(6,'(''Output file name: '',a)') filename
                  cycle
                  140 write(6,'(''Use: set filename name"'')')
               else
                  write(6,'(''Use: "set imin/filename"'')')
               endif
            else if (line(1:4).eq.'exit') then
               exit
            else if (line(1:2).ne.'') then
               write(6,'(''Available commands are:'')')
               write(6,'(''"set [arg]": set various things'')')
               write(6,'(''"analize obs": plot analysis of an observable'')')
               write(6,'(''"ave nblk": calculate averages of data'')')
               write(6,'(''"corr limit": calculate autocorrelations of data'')')
               write(6,'(''"list": list all the observables'')')
               write(6,'(''"errors": calculate errors for different blocks length'')')
               write(6,'(''print [arg]: save data in a file'')')
               write(6,'(''plot [arg]: plot data'')')
               write(6,'(''exit: pretty obvious'')')
            endif
         enddo
      else
! try to estimate the number of observables per line, number of points, and number of data
         rewind 10
         ndata=0
         ntab=0
         nobs=0
         endoffile=.false.
         endoftable=.false.
         endofline=.false.
         do while (.not.endoffile)
            read(10,'(a90)',err=1,end=1) line
            read(10,'(a90)',err=1,end=1) line
            read(10,'(a90)',err=1,end=1) line
            endoftable=.false.
            ntab=0
            do while(.not.endoftable)
               read(10,'(a500)') longline
               read(longline,*,err=2,end=2) dummy ! I assume the first number is r or q or something similar
               longline=longline(index(longline,' '):len(longline))
               longline=adjustl(longline)
               do while (.not.endofline)
                  read(longline,*,err=3,end=3) dummy
                  longline=longline(index(longline,' '):len(longline))
                  longline=adjustl(longline)
                  nobs=nobs+1
                  cycle
                  3 endofline=.true.
               enddo
               ntab=ntab+1
               cycle
               2 endoftable=.true.
            enddo
            read(10,'(a90)',err=1,end=1) line
            ndata=ndata+1
            cycle
            1 endoffile=.true.
         enddo
         nobs=(nobs-1)/2 ! first line is r or q, and then each observable is followed by an error
         write(6,'(''Number of observables per line = '',i10)') nobs
         write(6,'(''Number of points = '',i10)') ntab
         write(6,'(''Number of data = '',i10)') ndata
         rewind 10
         allocate(label(0:nobs),weight(ndata))
         allocate(r(ntab),rho(ndata,ntab,nobs),rhoerr(ndata,ntab,nobs))
         allocate(averho(ntab,nobs),errrho(ntab,nobs))
         do i=1,ndata
            read(10,'(a90)') line
            read(line,'(a1)') cc
            line=line(index(line,' '):len(line))
            line=adjustl(line)
            read(line,'(a2)') label(0)
            line=line(index(line,' '):len(line))
            line=adjustl(line)
            do j=1,nobs
               read(line,*) label(j)
               line=line(index(line,' '):len(line))
               line=adjustl(line)
               line=line(index(line,' '):len(line))
               line=adjustl(line)
            enddo   
            read(10,'(a90)') line
            line=line(index(line,'=')+1:len(line))
            line=adjustl(line)
            read(line,*) weight(i)
            read(10,'(a90)') line
            do j=1,ntab
               read(10,'(a500)') longline
               read(longline,*) r(j),(rho(i,j,k),rhoerr(i,j,k),k=1,nobs)
            enddo
            read(10,'(a90)') line
            read(10,'(a90)') line
         enddo
         close(10)
         filename=' '
         imin=1
         obsidx=1
         write(6,'(''Setup done!'')')
         write(6,'(''Type help if you dont know what to do'')')
         do while(.true.)
            write(6,'(''---> '',$)')
            read(5,'(a90)') line
            if (line(1:3).eq.'ave') then
               line=line(index(line,' '):len(line))
               line=adjustl(line)
               read(line,*,err=300,end=300) nblk
               call avecalc2(rho,weight,imin,ndata,ntab,nobs,averho,errrho,nblk)
               cycle
               300 write(6,'(''Use: "ave nblk"'')')
            else if (line(1:4).eq.'plot') then
               line=line(index(line,' '):len(line))
               line=adjustl(line)
               if (line(1:3).eq.'obs') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=310,end=310) ndat1,ndat2
                  call plot2xydy(r,rho(ndat1,:,obsidx),rhoerr(ndat1,:,obsidx),rho(ndat2,:,obsidx),rhoerr(ndat2,:,obsidx))
                  cycle
                  310 write(6,'(''Use: "plot obs dataset1 dataset2"'')')
               else
                  call plotxydy(r,averho(:,obsidx),errrho(:,obsidx))
               endif
            else if (line(1:3).eq.'set') then
               line=line(index(line,' '):len(line))
               line=adjustl(line)
               if (line(1:4).eq.'imin') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=350,end=350) n
                  imin=max(1,n)
                  if (imin.gt.ndata) then
                     write(6,'(''imin cannot be larger than ndata!'')')
                     imin=ndata
                  endif
                  write(6,'(''imin set to '',i10)') imin
                  write(6,'(''Number of data '',i10)') ndata-imin+1
                  cycle
                  350 write(6,'(''Use: set imin val"'')')
               else if (line(1:8).eq.'filename') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=360,end=360) filename
                  write(6,'(''Output file name: '',a)') filename
                  cycle
                  360 write(6,'(''Use: set filename name"'')')
               else if (line(1:3).eq.'obs') then
                  line=line(index(line,' '):len(line))
                  line=adjustl(line)
                  read(line,*,err=370,end=370) obsidx
                  if (obsidx.lt.1) obsidx=1
                  if (obsidx.gt.nobs) obsidx=nobs
                  write(6,'(''Analizing in datafile observable number: '',i10)') obsidx
                  cycle
                  370 write(6,'(''Use: set obs val'')') obsidx
               else
                  write(6,'(''Use: "set imin/filename"'')')
               endif
            else if (line(1:5).eq.'print') then
               call showobs2(nobs,ntab,label,r,averho(:,:),errrho(:,:),filename)
               write(6,'(/''Note: calculate averages before list if you have not done yet!'')')
               write(6,'(''List written to file '',a)') filename
            else if (line(1:4).eq.'exit') then
               exit
            else if (line(1:2).ne.'') then
               write(6,'(''Available commands are:'')')
               write(6,'(''"set [arg]": set various things'')')
               write(6,'(''"ave nblk": calculate averages of data'')')
               write(6,'(''"print": ave all the observables in a file'')')
               write(6,'(''plot [arg]: plot data'')')
               write(6,'(''exit: pretty obvious'')')
            endif
         enddo
      endif
   else if (domix) then
      open(unit=10,file=filename1)
      open(unit=11,file=filename2)
      write(6,'(''Crunch program''/)')
      write(6,'(''Name of the VMC file: '',a80,/)') filename1
      write(6,'(''Name of the DMC file: '',a80,/)') filename2
! read observables until the end of the two files
! note that the might have different observables
      nobs=64 ! this is the maximum, if you need more, increase this number
      allocate(label2(2,nobs))
      allocate(obs(2,nobs),err(2,nobs))
      iobs=0
      do ic=1,2
         do i=1,nobs
            read(10+ic-1,'(a90)',err=500,end=500) line
            read(line,*) idummy
            line=line(index(line,' '):len(line))
            line=adjustl(line)
            label2(ic,i)=line(1:32)
            line=line(index(line,'  '):len(line))
            line=adjustl(line)
            read(line,*) obs(ic,i),dummy,err(ic,i)
            iobs(ic)=iobs(ic)+1
            cycle
            500 exit
         enddo
      enddo
      write(6,'(''Reading files done!''/)')
      write(6,'(''Observables of VMC file:'',i10)') iobs(1)
      do i=1,iobs(1)
         write(6,'(a25,2e19.10)') label2(1,i),obs(1,i),err(1,i)
      enddo
      write(6,'(/''Observables of DMC file:'',i10)') iobs(2)
      do i=1,iobs(2)
         write(6,'(a25,2e19.10)') label2(2,i),obs(2,i),err(2,i)
      enddo
      write(6,'(/''VMC, DMC, and extrapolated observables 2*DMC-vmc:''/)')
      do ivmc=1,iobs(1)
         do idmc=1,iobs(2)
            if (label2(2,idmc).eq.label2(1,ivmc)) then
               extobs=2.0_r8*obs(2,idmc)-obs(1,ivmc)
               exterr=err(2,idmc)+err(1,ivmc)
               write(6,'(a25,6f13.5)') label2(1,ivmc),obs(1,ivmc),err(1,ivmc),obs(2,idmc),err(2,idmc),extobs,exterr
            endif
         enddo
      enddo
   else if (domix2) then
      open(unit=10,file=filename1)
      open(unit=11,file=filename2)
      open(unit=12,file=filename3)
      write(6,'(''Crunch program''/)')
      write(6,'(''Name of the VMC file: '',a80,/)') filename1
      write(6,'(''Name of the DMC file: '',a80,/)') filename2
      write(6,'(''Name of the output file: '',a80,/)') filename3
! read observables until the end of the two files
! note that they must have the same observables and same number of points!
      ntab=0
      nobs=0
      endoftable=.false.
      endofline=.false.
      read(10,'(a500)') longline
      do while(.not.endoftable)
         read(10,'(a500)',err=600,end=600) longline
         read(longline,*,err=600,end=600) dummy ! I assume the first number is r or q or something similar
         longline=longline(index(longline,' '):len(longline))
         longline=adjustl(longline)
         do while (.not.endofline)
            read(longline,*,err=620,end=620) dummy
            longline=longline(index(longline,'  '):len(longline))
            longline=adjustl(longline)
            nobs=nobs+1
            cycle
            620 endofline=.true.
         enddo
         ntab=ntab+1
         cycle
         600 endoftable=.true.
      enddo
      nobs=(nobs-1)/2 ! first line is r or q, and then each observable is followed by an error
      write(6,'(''Number of observables per line = '',i10)') nobs
      write(6,'(''Number of points = '',i10)') ntab
      allocate(r(ntab),rho(2,ntab,nobs),rhoerr(2,ntab,nobs))
      allocate(extrho(ntab,nobs),extrhoerr(ntab,nobs))
      rewind 10
      do ic=1,2
         read(10+ic-1,'(a500)') longline
         do j=1,ntab
            read(10+ic-1,'(a500)') longline
            read(longline,*) r(j),(rho(ic,j,k),rhoerr(ic,j,k),k=1,nobs)
         enddo
      enddo
      write(6,'(''Reading files done!''/)')
      write(6,'(/''VMC, DMC, and extrapolated observables 2*DMC-vmc:''/)')
     extrho(:,:)=2.0_r8*rho(2,:,:)-rho(1,:,:)
!      extrho(:,:)=rho(2,:,:)**2/(rho(1,:,:)+1.0e-12)
      extrhoerr(:,:)=sqrt(rhoerr(1,:,:)**2+4.0_r8*(rhoerr(2,:,:)**2))
!extrhoerr(:,:)=sqrt((rhoerr(1,:,:)**2)*(rho(2,:,:)**4)/(rho(1,:,:)**4)+ &
!4.0_r8*(rho(2,:,:)**2)*(rhoerr(2,:,:)**2)/(rho(1,:,:)**2))
      rewind 10
      read(10,'(a500)') longline
      write(12,*) trim(longline)
      do i=1,ntab
         write(12,'(100e19.10)') r(i),(extrho(i,j),extrhoerr(i,j),j=1,nobs)
!        write(12,'(100e19.10)') r(i),(rho(1,i,j),rhoerr(1,i,j),j=1,nobs)
      enddo
      close(10)
      close(11)
      close(12)
   endif
   write(6,'(''Finish!'')')
   contains
      subroutine findw(word,n)
      implicit none
      character(len=*) :: word
      character(len=30) :: readw
      integer(kind=i4) :: n
      readw=''
      do while (readw(1:len(word)).ne.word)
         read(10,'(a30)',end=10,err=10) readw
         n=n+1
      enddo
      return
      10 continue
      n=-1 
      end subroutine findw

      subroutine avecalc(obs,imin,imax,aveobs,errobs,nobs,nblk,errblk)
      implicit none
      integer(kind=i4) :: nobs
      real(kind=r8) :: obs(:,:),aveobs(:,:),errobs(:,:),errblk(:)
      real(kind=r8) :: ave(nobs),err(nobs),oblk(nobs)
      integer(kind=i4) :: imin,imax,i,j,idx,nblk,nblks,idx2,ndata
      write(6,'(''avecalc: min, max, nblk values are = '',3i7)') imin,imax,nblk
      aveobs=0.0_r8
      errobs=0.0_r8
      ave=0.0_r8
      err=0.0_r8
      idx=1
      do i=imin,imax
         ave(:)=ave(:)+obs(:,i)
         err(:)=err(:)+obs(:,i)**2
         aveobs(:,i)=ave/idx
         errobs(:,i)=err/idx
         errobs(:,i)=sqrt(abs(aveobs(:,i)**2-errobs(:,i))/idx)
         idx=idx+1
      enddo
! calculate errors using blocks of data with size nblk
! some data could be not included
      ndata=imax-imin+1
      nblks=ndata/nblk
      idx2=imin
      ave=0.0_r8
      err=0.0_r8
      do i=1,nblks
         oblk=0.0_r8
         do j=1,nblk
            oblk=oblk+obs(:,idx2)
            idx2=idx2+1
         enddo
         oblk=oblk/nblk
         ave=ave+oblk
         err=err+oblk**2
      enddo
      ave=ave/nblks
      errblk(:)=err/nblks
      errblk(:)=sqrt(abs(ave(:)**2-errblk(:))/nblks)
      end subroutine avecalc

      subroutine avecalc2(rho,wt,imin,ndata,ntab,nobs,averho,errrho,nblk)
      implicit none
      integer(kind=i4) :: ndata,ntab,nobs,nblk
      real(kind=r8) :: rho(:,:,:),averho(:,:),errrho(:,:),wt(:)
      real(kind=r8) :: ave(ntab,nobs),err(ntab,nobs)
      real(kind=r8) :: obl(ntab,nobs),wtsum
      integer(kind=i4) :: imin,i,j,idx,nblks
      averho=0.0_r8
      errrho=0.0_r8
      wtsum=0.0_r8
      do i=imin,ndata
         do j=1,ntab
            averho(j,:)=averho(j,:)+wt(i)*rho(i,j,:)
         enddo
         wtsum=wtsum+wt(i)
      enddo
      averho=averho/wtsum
! calculate errors using blocks of data with size nblk
! some data could be not included
      nblks=(ndata-imin+1)/nblk
      idx=imin
      ave=0.0_r8
      err=0.0_r8
      do i=1,nblks
         obl=0.0_r8
         wtsum=0.0_r8
         do j=1,nblk
            obl=obl+wt(idx)*rho(idx,:,:)
            wtsum=wtsum+wt(idx)
            idx=idx+1
         enddo
         obl=obl/wtsum
         ave=ave+obl
         err=err+obl**2
      enddo
      ave=ave/nblks
      err=err/nblks
      errrho=sqrt(abs(ave**2-err)/nblks)
      end subroutine avecalc2

      subroutine prnt(filename,nobs,label,aveobs,errobs,imin,imax,step,imtime)
      implicit none
      integer(kind=i4) :: unit,nobs,i,j,imin,imax,step
      character(len=25) :: label(:)
      character(len=*) :: filename
      real(kind=r8) :: aveobs(:,:),errobs(:,:),imtime(:)
      open(unit=10,file=filename)
      do i=imin,imax,step
         do j=1,nobs
            write(10,'(a25,1x,1p,3e19.10)') label(j),imtime(i),aveobs(j,i),errobs(j,i)
         enddo
      enddo
      close(10)
      end subroutine prnt

      subroutine listobs(nobs,label)
      implicit none
      integer(kind=i4) :: nobs,i
      character(len=25) :: label(:)
      do i=1,nobs
         write(6,'(i5,3x,a25)') i,label(i)
      enddo
      end subroutine listobs

      subroutine autocorr(nobs,obs,imin,imax,ref,corr,label)
      implicit none
      real(kind=r8) :: obs(:,:),ref
      character(len=25) :: label(:)
      integer(kind=i4) :: nobs,imin,imax,i,j,ndata
      real(kind=r8) :: xi(nobs),xxi(nobs),xj(nobs),xxj(nobs),xij(nobs)
      real(kind=r8) :: corr(:,:)
      logical :: idxfound(nobs)
      integer(kind=i4) :: idx(nobs)
      idxfound=.false.
      corr=0.0_r8
      ndata=imax-imin+1
      idx=1
      do i=1,ndata/2 ! this is arbitrary, but should be large enough
         xi=0.0_r8
         xj=0.0_r8
         xij=0.0_r8
         xxi=0.0_r8
         xxj=0.0_r8
         do j=imin,imax-i+1
            xi=xi+obs(:,j)
            xxi=xxi+obs(:,j)**2
            xij=xij+obs(:,j)*obs(:,j+i-1)
            xj=xj+obs(:,j+i-1)
            xxj=xxj+obs(:,j+i-1)**2
         enddo
         xi=xi/(ndata-i+1)
         xxi=xxi/(ndata-i+1)
         xij=xij/(ndata-i+1)
         xj=xj/(ndata-i+1)
         xxj=xxj/(ndata-i+1)
         corr(:,i)=(xij-xi*xj)/sqrt((xxi-xi**2)*(xxj-xj**2))
         do j=1,nobs
            if (.not.idxfound(j)) then
               if (corr(j,i).gt.ref) then
                  idx(j)=idx(j)+1
               else
                  idxfound(j)=.true.
               endif
            endif
         enddo
      enddo
      write(6,'(''autocorrelations: reference is = '',f5.2)') ref
      do i=1,nobs
         write(6,'(''observable '',i3,3x,a25,'' index = '',i5)') i,label(i),idx(i)
      enddo
      end subroutine autocorr

      subroutine prntcorr(filename,corr,nobs,ndata)
      implicit none
      character(len=*) :: filename
      real(kind=r8) :: corr(:,:)
      integer(kind=i4) :: unit,nobs,i,ndata
      open(unit=10,file=filename)
      do i=1,ndata
         write(10,*) corr(nobs,i)
      enddo
      close(10)
      end subroutine prntcorr

      subroutine showobs(nobs,label,aveobs,errobs,errblk,filename)
      implicit none
      integer(kind=i4) :: nobs,i
      character(len=25) :: label(:)
      character(len=*) :: filename
      real(kind=r8) :: aveobs(:),errobs(:),errblk(:)
      logical :: filex
      filex=.false.
      if (filename(1:2).ne.' ') then
         filex=.true.
         open(unit=10,file=filename)
      endif
      write(6,'(''observable, average, bad error, good error'')') 
      do i=1,nobs
         write(6,'(i5,1x,a25,1x,1p,3e19.10)') i,label(i),aveobs(i),errobs(i),errblk(i)
         if (filex) write(10,'(i5,1x,a25,1x,1p,3e19.10)') i,label(i),aveobs(i),errobs(i),errblk(i)
      enddo
      if (filex) close(10)
      end subroutine showobs

      subroutine showobs2(nobs,ntab,label,r,rho,err,filename)
      implicit none
      integer(kind=i4) :: nobs,ntab,i,j
      character(len=25) :: label(0:nobs)
      character(len=*) :: filename
      real(kind=r8) :: r(:),rho(:,:),err(:,:)
      logical :: filex
      filex=.false.
      if (filename(1:2).ne.' ') then
         filex=.true.
         open(unit=10,file=filename)
         write(10,'(a,100a10)') "#",(label(j),j=0,nobs)
         do i=1,ntab
            write(10,'(100e19.10)') r(i),(rho(i,j),err(i,j),j=1,nobs)
         enddo
         close(10)
      else
         write(6,'(''Error, showobs2: a filename should be set first!'')')
      endif
      end subroutine showobs2

      subroutine errors(obs,imin,imax,errobs,nobs)
      implicit none
      integer(kind=i4) :: nobs,ndata
      real(kind=r8) :: obs(:,:),errobs(:,:),errblk(nobs)
      real(kind=r8) :: ave(nobs),err(nobs),oblk(nobs)
      integer(kind=i4) :: imin,imax,i,j,idx,nblk,nblks
! calculate errors using blocks of data with size nblk
! some data could be not included
      ndata=imax-imin+1
      errobs=0.0_r8
      do nblk=1,ndata/2
         nblks=ndata/nblk
         idx=imin
         ave=0.0_r8
         err=0.0_r8
         do i=1,nblks
            oblk=0.0_r8
            do j=1,nblk
               oblk=oblk+obs(:,idx)
               idx=idx+1
            enddo
            oblk=oblk/nblk
            ave=ave+oblk
            err=err+oblk**2
         enddo
         ave=ave/nblks
         errblk(:)=err/nblks
         errblk(:)=sqrt(abs(ave(:)**2-errblk(:))/nblks)
         errobs(:,nblk)=errblk(:)
      enddo
      end subroutine errors

      subroutine prnterr(filename,errobs,nobs,ndata)
      implicit none
      character(len=*) :: filename
      integer(kind=i4) :: unit,nobs,ndata,i
      real(kind=r8) :: errobs(:,:)
      open(unit=10,file=filename)
      do i=1,ndata/10
         write(10,'(i10,e15.8)') i,errobs(nobs,i)
      enddo
      close(10)
      end subroutine prnterr

      subroutine ploty(y)
      implicit none
      real(kind=8) :: y(:)
      integer(kind=i4) :: i
      open (unit=99,file='tmp.dat',status='replace')
      do i=1,size(y)
         if (y(i).ne.0.0_r8) write (99,'(3e15.7)') y(i)
      enddo
      close(99)
      open(unit=99,file='tmp.plt',status='replace')
      write(99,'(a)') 'set terminal x11'//' title  "Gnuplot"' 
      write(99,'(a)') 'unset key'
      write(99,'(a)') 'plot "tmp.dat" w lp'
      write(99,'(a)') 'pause -1'
      write(99,'(a)') 'q'
      close(99)
      call gnuplot('tmp.plt') 
      end subroutine ploty

      subroutine plotxy(x,y)
      implicit none
      real(kind=8) :: x(:),y(:)
      integer(kind=i4) :: i
      open (unit=99,file='tmp.dat',status='replace')
      do i=1,size(x)
         write (99,'(3e15.7)') x(i),y(i)
      enddo
      close(99)
      open(unit=99,file='tmp.plt',status='replace')
      write(99,'(a)') 'set terminal x11'//' title  "Gnuplot"' 
      write(99,'(a)') 'unset key'
      write(99,'(a)') 'plot "tmp.dat" u 1:2 w lp'
      write(99,'(a)') 'pause -1'
      write(99,'(a)') 'q'
      close(99)
      call gnuplot('tmp.plt') 
      end subroutine plotxy

      subroutine plotxydy(x,y,err)
      implicit none
      real(kind=8) :: x(:),y(:),err(:)
      integer(kind=i4) :: i
      open (unit=99,file='tmp.dat',status='replace')
      do i=1,size(x)
         write (99,'(3e15.7)') x(i),y(i),err(i)
      enddo
      close(99)
      open(unit=99,file='tmp.plt',status='replace')
      write(99,'(a)') 'set terminal x11'//' title  "Gnuplot"' 
      write(99,'(a)') 'unset key'
      write(99,'(a)') 'plot "tmp.dat" u 1:2:3 w err'
      write(99,'(a)') 'pause -1'
      write(99,'(a)') 'q'
      close(99)
      call gnuplot('tmp.plt') 
      end subroutine plotxydy

      subroutine plot2xydy(x,y1,err1,y2,err2)
      implicit none
      real(kind=8) :: x(:),y1(:),err1(:),y2(:),err2(:)
      integer(kind=i4) :: i
      open (unit=99,file='tmp.dat',status='replace')
      do i=1,size(x)
         write (99,'(5e15.7)') x(i),y1(i),err1(i),y2(i),err2(i)
      enddo
      close(99)
      open(unit=99,file='tmp.plt',status='replace')
      write(99,'(a)') 'set terminal x11'//' title  "Gnuplot"' 
      write(99,'(a)') 'plot "tmp.dat" u 1:2:3 w err'
      write(99,'(a)') 'replot "tmp.dat" u 1:4:5 w err'
      write(99,'(a)') 'pause -1'
      write(99,'(a)') 'q'
      close(99)
      call gnuplot('tmp.plt') 
      end subroutine plot2xydy

      subroutine plotall(x,y1,err1,y2,err2,y3,y4)
      implicit none
      real(kind=8) :: x(:),y1(:),err1(:),y2(:),err2(:),y3(:),y4(:)
      integer(kind=i4) :: i
      open (unit=99,file='tmp.dat',status='replace')
      do i=1,size(x)
         write(99,'(5e15.7)') x(i),y1(i),err1(i),y2(i),err2(i)
      enddo
      write(99,'(//)')
      do i=1,size(y3)
         if (y3(i).ne.0.0_r8) write(99,'(e15.7)') y3(i)
      enddo
      write(99,'(//)')
      do i=1,size(y4)
         if (y4(i).ne.0.0_r8) write(99,'(e15.7)') y4(i)
      enddo
      close(99)
      ! plot data into x11 terminal
      open(unit=99,file='tmp.plt',status='replace')
      write(99,'(a)') 'set terminal x11 0'//' title "observable and average"'
      write(99,'(a)') 'plot "tmp.dat" i 0 u :2:3 w err' !1
      write(99,'(a)') 'replot "tmp.dat" i 0 u :4:5 w err' !1
      write(99,'(a)') 'set terminal x11 1'//' title "autocorrelations"'
      write(99,'(a)') 'plot "tmp.dat" i 1 w lp'
      write(99,'(a)') 'set terminal x11 2'//' title "errors as a funcion of block size"'
      write(99,'(a)') 'plot "tmp.dat" i 2 w lp'
      write(99,'(a)') 'pause -1'
      write(99,'(a)') 'q'
      close(99)
      call gnuplot('tmp.plt')
      ! plot data into plot.pdf file
      open(unit=8,file='tmp2.plt',status='replace')
      write(8,'(a)') "set terminal pngcairo size 1024,2160"
      write(8,'(a)') "set output 'plots.png'"
      write(8,'(a)') "set lmargin at screen 0.15"
      write(8,'(a)') "set rmargin at screen 0.95"
      write(8,'(a)') "set multiplot layout 3,1 scale 1,0.333"
      write(8,'(a)') "set size 1,0.33333"
      write(8,'(a)') "set origin 0,0.66666"
      write(8,'(a)') "set title 'Observable and average'"
      write(8,'(a)') "plot 'tmp.dat' i 0 u :2:3 w err t 'Observable', \"
      write(8,'(a)') "     'tmp.dat' i 0 u :4:5 w err t 'Average'"
      write(8,'(a)') "set size 1,0.33333"
      write(8,'(a)') "set origin 0,0.33333"
      write(8,'(a)') "set title 'Autocorrelations'"
      write(8,'(a)') "unset key"
      write(8,'(a)') "plot 'tmp.dat' i 1 w lp t 'Autocorrelations'"
      write(8,'(a)') "set size 1,0.33333"
      write(8,'(a)') "set origin 0,0"
      write(8,'(a)') "set title 'Errors as a function of block size'"
      write(8,'(a)') "plot 'tmp.dat' i 2 w lp t 'Errors as a function of block size'"
      close(8)
      call gnuplot('tmp2.plt')
      call system('convert plots.png -quality 100 plots.pdf')
      call system('rm plots.png')
      end subroutine plotall

      subroutine gnuplot(cmd)
      implicit none
      character(len=100) :: command
      character(len=*) :: cmd
      integer :: status,system
      write(command,*) '/usr/bin/gnuplot '//trim(cmd)
      status=system(trim(command))
      if (status.ne.0) write(6,'(''Error running Gnuplot!'')')
      end subroutine gnuplot
   end program crunch
