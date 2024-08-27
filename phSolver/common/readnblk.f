c  readnblk.f (pronounce "Reed and Block Dot Eff") contains:
c
c    module readarrays ("Red Arrays") -- contains the arrays that
c     are read in from binary files but not immediately blocked 
c     through pointers.
c
c    subroutine readnblk ("Reed and Block") -- allocates space for
c     and reads data to be contained in module readarrays.  Reads
c     all remaining data and blocks them with pointers.
c


      module readarrays
      
      real*8, allocatable :: point2x(:,:)
      real*8, allocatable :: qold(:,:)
      real*8, allocatable :: uold(:,:)
      real*8, allocatable :: acold(:,:)
      integer, allocatable :: iBCtmp(:)
      real*8, allocatable :: BCinp(:,:)

      !integer*8, allocatable :: fncorp(:) !MB
      !integer, allocatable :: ltg(:)  !MB

      integer, allocatable :: point2ilwork(:)
      integer, allocatable :: nBC(:)
      integer, allocatable :: point2iper(:)
      integer, allocatable :: point2ifath(:)
      integer, allocatable :: point2nsons(:)
      
      end module



      subroutine readnblk
c
      use readarrays
      include "common.h"
c
      real*8, allocatable :: xread(:,:), qread(:,:), acread(:,:)
      real*8, allocatable :: uread(:,:)
      real*8, allocatable :: BCinpread(:,:)
      integer, allocatable :: iperread(:), iBCtmpread(:)
      integer, allocatable :: ilworkread(:), nBCread(:)
      !integer, target, allocatable :: fncorpread(:)   ! MB
      !integer fncorpsize  ! MB
      character*10 cname2
      character*8 mach2
      character*30 fmt1
      character*255 fname1,fnamer,fnamelr
      character*255 warning
      integer igeomBAK, ibndc, irstin, ierr
      integer intfromfile(50) ! integers read from headers

!     SyncIO parameters
      integer :: descriptor, descriptorG, GPID, color, nfiles, nfields
      integer ::  numparts, nppf, nl_numstart, readCode
      integer :: ierr_io, numprocs, itmp, itmp2
      character*255 fname2, temp2
      character*64 temp1

c
c Assign number of error statistics
c
      numerr = 10

c
c.... determine the step number to start with
c
      open(unit=72,file='numstart.dat',status='old')
      read(72,*) irstart
      close(72)
      lstep=irstart ! in case restart files have no fields

      itmp=1
      if (irstart .gt. 0) itmp = int(log10(float(irstart)))+1
      write (fmt1,"('(''restart.'',i',i1,',1x)')") itmp
      write (fnamer,fmt1) irstart
      fnamer = trim(fnamer) // cname2(myrank+1)
c      fnamelr = trim(fnamelr) // cname2(myrank+1)

c
c.... open input files
c

c
c.... input the geometry parameters
c

      nfiles = nsynciofiles
      numparts = numpe !This is the common settings. Beware if you try to compute several parts per process

      color = int(myrank/(numparts/nfiles)) !Should call the color routine in SyncIO here
      itmp2 = int(log10(float(color+1)))+1
      write (temp2,"('(''geombc-dat.'',i',i1,')')") itmp2
      write (fnamer,temp2) (color+1)
      fnamer = trim(fnamer) // char(0)

      itwo=2
      ione=1
      ieleven=11
      itmp = int(log10(float(myrank+1)))+1

      call queryphmpiio(fnamer, nfields, nppf);
      if (myrank == 0) then
        write(*,*) 'Number of fields in geombc-dat: ',nfields
        write(*,*) 'Number of parts per file geombc-dat: ',nppf
      endif
      call initphmpiio( nfields, nppf, nfiles, igeom,
     & 'read' // char(0))
      call openfile( fnamer, 'read' // char(0), igeom )

      write (temp1,"('(''number of nodes@'',i',i1,',A1)')") itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2 // char(0),numnp,ione,
     & 'integer' // char(0), iotype)

      write (temp1,"('(''number of modes@'',i',i1,',A1)')") itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2 // char(0),nshg,ione,
     & 'integer' // char(0), iotype)

      write (temp1,"('(''number of interior elements@'',i',i1,',A1)')")
     &       itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2 // char(0),numel,ione,
     & 'integer' // char(0), iotype)

      write (temp1,"('(''number of boundary elements@'',i',i1,',A1)')")
     &       itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2 // char(0),numelb,ione,
     & 'integer' // char(0),iotype)

      write (temp1,
     & "('(''maximum number of element nodes@'',i',i1,',A1)')") itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2 // char(0),nen,ione,
     &'integer' // char(0),iotype)

      write (temp1,"('(''number of interior tpblocks@'',i',i1,',A1)')")
     &       itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2 // char(0),nelblk,ione,
     & 'integer' // char(0) ,iotype)

      write (temp1,"('(''number of boundary tpblocks@'',i',i1,',A1)')")
     &       itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2 // char(0),nelblb,ione,
     & 'integer' // char(0), iotype)

      write (temp1,
     & "('(''number of nodes with Dirichlet BCs@'',i',i1,',A1)')") itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2 // char(0),numpbc,ione,
     & 'integer' // char(0),iotype)

      write (temp1,"('(''number of shape functions@'',i',i1,',A1)')")
     &       itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2 // char(0),ntopsh,ione,
     & 'integer' // char(0),iotype)

c
c.... calculate the maximum number of boundary element nodes
c     
      nenb = 0
      do i = 1, melCat
         if (nen .eq. nenCat(i,nsd)) nenb = max(nenCat(i,nsd-1), nenb)
      enddo
c     
      if (myrank == master) then
         if (nenb .eq. 0) call error ('input   ','nen     ',nen)
      endif
c
c.... setup some useful constants
c
      I3nsd  = nsd / 3          ! nsd=3 integer flag
      E3nsd  = float(I3nsd)     ! nsd=3 real    flag
c    
      if(matflg(1,1).lt.0) then	! incompressible
         nflow = nsd + 1
      else			            ! compressible
         nflow = nsd + 2
      endif 
      ndof   = nsd + 2
      write(*,*) 'nsd is ', nsd
      nsclr=impl(1)/100
      write(*,*) 'nsclr is ', nsclr
      ndof=ndof+nsclr           ! number of sclr transport equations to solve
      
      ndofBC = ndof + I3nsd     ! dimension of BC array
      ndiBCB = 2                ! dimension of iBCB array
      ndBCB  = ndof + 1         ! dimension of BCB array
c     
      nsymdf = (ndof*(ndof + 1)) / 2 ! symm. d.o.f.'s
c
c.... ----------------------> Communication tasks <--------------------
c

      if(numpe > 1) then

         write (temp1,"('(''size of ilwork array@'',i',i1,',A1)')") itmp
         write (fname2,temp1) (myrank+1),'?'
         call readheader(igeom,fname2 // char(0),nlwork,ione,
     &   'integer' // char(0) ,iotype)

         write (temp1,"('(''ilwork@'',i',i1,',A1)')") itmp
         write (fname2,temp1) (myrank+1),'?'
         call readheader(igeom,fname2 //char(0) ,nlwork,ione,
     &   'integer' // char(0) , iotype)

         allocate( point2ilwork(nlwork) )
         allocate( ilworkread(nlwork) )
         call readdatablock(igeom,fname2 // char(0),ilworkread,
     &                      nlwork,'integer' // char(0) , iotype)

         point2ilwork = ilworkread
         call ctypes (point2ilwork)
! ######################################################################
        ! MB here goes the code to aptly assign the ltg arrays...
!         if(svLSFlag.eq.1) then
!            fncorpsize = nshg
!            allocate(fncorp(fncorpsize))
!            call gen_ncorp(fncorp, ilworkread, nlwork, fncorpsize)
   !
   ! the  following code finds the global range of the owned nodes
   !
!            maxowned=0
!            minowned=maxval(fncorp)
!            do i = 1,nshg      
!               if(fncorp(i).gt.0) then  ! don't consider remote copies
!                  maxowned=max(maxowned,fncorp(i))
!                  minowned=min(minowned,fncorp(i))
!               endif
!            enddo
   !
   !  end of global range code
   !
!            call commuInt(fncorp, point2ilwork, 1, 'out')
!            ncorpsize = fncorpsize 
!          endif
!          write(*,*) 'svLSFlag is ', svLSFlag
!          if(svLSFlag.eq.1) then
!           allocate(ltg(ncorpsize))
!            write(*,*) 'norpsize is ', ncorpsize
!            ltg(:)=fncorp(:)
!            write(*,*) 'ltg(:) is ', ltg
!          endif
! ######################################################################

      else
         nlwork = 1
         allocate( point2ilwork(nlwork) )
         nshg0 = nshg
         point2ilwork= 1
      endif

c     
c.... read the node coordinates
c

      itwo=2
      write (temp1,"('(''co-ordinates@'',i',i1,',A1)')") itmp
      write (fname2,temp1) (myrank+1),'?'

      call readheader(igeom,fname2 // char(0),intfromfile,itwo,
     & 'double' // char(0), iotype)
      numnp=intfromfile(1)
      allocate( point2x(numnp,nsd) )
      allocate( xread(numnp,nsd) )
      ixsiz=numnp*nsd
      call readdatablock(igeom,fname2 // char(0),xread,ixsiz,
     & 'double' // char(0), iotype)
      point2x = xread
c
c.... read in and block out the connectivity
c
      call genblk (IBKSIZ)
c
c.... read the boundary condition mapping array
c
      ione=1
      write (temp1,"('(''bc mapping array@'',i',i1,',A1)')") itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2 // char(0),nshg,ione,
     & 'integer' // char(0),iotype)
      allocate( nBC(nshg) )

      allocate( nBCread(nshg) )

      call readdatablock(igeom,fname2 // char(0),nBCread,nshg,
     & 'integer' // char(0),iotype)

      nBC=nBCread

c
c.... read the temporary iBC array
c
      ione=1
      write (temp1,"('(''bc codes array@'',i',i1,',A1)')") itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2 // char(0) ,numpbc,ione,
     & 'integer' // char(0),iotype)

      if ( numpbc > 0 ) then
        allocate( iBCtmp(numpbc) )
        allocate( iBCtmpread(numpbc) )
      else
        allocate( iBCtmp(1) )
        allocate( iBCtmpread(1) )
      endif
      call readdatablock(igeom,fname2 // char(0),iBCtmpread,numpbc,
     &  'integer' // char(0),iotype)

      if ( numpbc > 0 ) then
         iBCtmp=iBCtmpread
      else  ! sometimes a partition has no BC's
         deallocate( iBCtmpread)
         iBCtmp=0
      endif
c
c.... read boundary condition data
c
      ione=1
      write (temp1,"('(''boundary condition array@'',i',i1,',A1)')")
     &     itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2 // char(0),intfromfile,
     &     ione, 'double' // char(0), iotype)
c here intfromfile(1) contains (ndof+7)*numpbc
      if ( numpbc > 0 ) then
!         if(intfromfile(1).ne.(ndof+7)*numpbc) then
!           warning='WARNING more data in BCinp than needed: keeping 1st'
!           write(*,*) warning, ndof+7
!         endif
         allocate( BCinp(numpbc,ndof+7) )
         nsecondrank=intfromfile(1)/numpbc
         allocate( BCinpread(numpbc,nsecondrank) )
         iBCinpsiz=intfromfile(1)
      else
         allocate( BCinp(1,ndof+7) )
         allocate( BCinpread(0,0) ) !dummy
         iBCinpsiz=intfromfile(1)
      endif
      call readdatablock(igeom,fname2 // char(0),BCinpread,iBCinpsiz,
     &                      'double' // char(0) ,iotype)

      if ( numpbc > 0 ) then
         BCinp(:,1:(ndof+7))=BCinpread(:,1:(ndof+7))
      else  ! sometimes a partition has no BC's
         deallocate(BCinpread)
         BCinp=0
      endif
c
c.... read periodic boundary conditions
c

      ione=1
      write (temp1,"('(''periodic masters array@'',i',i1,',A1)')") itmp
      write (fname2,temp1) (myrank+1),'?'

      call readheader(igeom,fname2 // char(0) ,nshg,
     &     ione, 'integer' // char(0), iotype)
      allocate( point2iper(nshg) )
      allocate( iperread(nshg) )
      call readdatablock(igeom,fname2 // char(0),iperread,nshg,
     &                      'integer' // char(0),iotype)
      point2iper=iperread
c
c.... generate the boundary element blocks
c
      call genbkb (ibksiz)

c
c  Read in the nsons and ifath arrays if needed
c
c  There is a fundamental shift in the meaning of ifath based on whether
c  there exist homogenous directions in the flow.  
c
c  HOMOGENOUS DIRECTIONS EXIST:  Here nfath is the number of inhomogenous
c  points in the TOTAL mesh.  That is to say that each partition keeps a 
c  link to  ALL inhomogenous points.  This link is furthermore not to the
c  sms numbering but to the original structured grid numbering.  These 
c  inhomogenous points are thought of as fathers, with their sons being all
c  the points in the homogenous directions that have this father's 
c  inhomogeneity.  The array ifath takes as an arguement the sms numbering
c  and returns as a result the father.
c
c  In this case nsons is the number of sons that each father has and ifath
c  is an array which tells the 
c
c  NO HOMOGENOUS DIRECTIONS.  In this case the mesh would grow to rapidly
c  if we followed the above strategy since every partition would index its
c  points to the ENTIRE mesh.  Furthermore, there would never be a need
c  to average to a node off processor since there is no spatial averaging.
c  Therefore, to properly account for this case we must recognize it and
c  inerrupt certain actions (i.e. assembly of the average across partitions).
c  This case is easily identified by noting that maxval(nsons) =1 (i.e. no
c  father has any sons).  Reiterating to be clear, in this case ifath does
c  not point to a global numbering but instead just points to itself.
c
      nfath=1  ! some architectures choke on a zero or undeclared
                 ! dimension variable.  This sets it to a safe, small value.
c      if (myrank.eq.master) write(*,*) 'nohomog = ', nohomog 
      if(((iLES .lt. 20) .and. (iLES.gt.0))
     &                   .or. (itwmod.gt.0) .or. (iDNS.gt.0)  ) then ! don't forget same
                                                    ! conditional in proces.f

c           read (igeom) nfath  ! nfath already read in input.f,
                                     ! needed for alloc
         ione=1
c         call creadlist(igeom,ione,nfath)
c         fname1='keyword sonfath?'
         if(nohomog.gt.0) then

            write (temp1,"('(''number of father-nodes@'',i',i1,',A1)')")
     &             itmp
            write (fname2,temp1) (myrank+1),'?'
            call readheader(igeom,fname2 // char(0),nfath,ione,
     &      'integer' // char(0), iotype)


            write (temp1,
     & "('(''number of son-nodes for each father@'',i',i1,',A1)')") itmp
            write (fname2,temp1) (myrank+1),'?'
            call readheader(igeom,fname2 // char(0),nfath,ione,
     &      'integer' // char(0), iotype)

            allocate (point2nsons(nfath))

            call readdatablock(igeom,fname2 // char(0),point2nsons,
     &                      nfath,'integer' // char(0), iotype)

c
            write (temp1,"('(''keyword ifath@'',i',i1,',A1)')") itmp
            write (fname2,temp1) (myrank+1),'?'
            call readheader(igeom,fname2 // char(0),nshg,ione,
     &      'integer' // char(0), iotype)

            allocate (point2ifath(nshg))
            call readdatablock(igeom,fname2 // char(0),point2ifath,
     &                      nshg,'integer' // char(0) , iotype)
     
c     
            nsonmax=maxval(point2nsons)
c
         else  ! this is the case where there is no homogeneity
               ! therefore ever node is a father (too itself).  sonfath
               ! (a routine in NSpre) will set this up but this gives
               ! you an option to avoid that.
            nfath=nshg
            allocate (point2nsons(nfath))
            point2nsons=1
            allocate (point2ifath(nshg))
            do i=1,nshg
               point2ifath(i)=i
            enddo
            nsonmax=1
c
         endif
      else
         allocate (point2nsons(1))
         allocate (point2ifath(1))
      endif

      call closefile( igeom, "read" // char(0) )
      call finalizephmpiio( igeom )

c
c  renumber the master partition for SPEBC
c
c      if((myrank.eq.master).and.(irscale.ge.0)) then
c         call setSPEBC(numnp, nfath, nsonmax)
c         call renum(point2x,point2ifath,point2nsons)
c      endif
c
c.... Read restart files
c

      itmp=1
      if (irstart .gt. 0) itmp = int(log10(float(irstart+1)))+1

      write (fmt1,"('(''restart-dat.'',i',i1,',1x)')") itmp

      write (fnamer,fmt1) irstart
      fnamer = trim(fnamer) // cname2(color+1)

      call queryphmpiio(fnamer // char(0), nfields, nppf);
      if (myrank == 0) then
        write(*,*) 'Number of fields in restart-dat: ',nfields
        write(*,*) 'Number of parts per file restart-dat: ',nppf
      endif
      call initphmpiio(nfields,nppf,nfiles,descriptor,
     & 'read' // char(0))
      call openfile( fnamer // char(0) , 
     & 'read' // char(0), descriptor )

      ithree=3
c      call creadlist(irstin,ithree,nshg2,ndof2,lstep)
c
c Only read time step if user has not specified time start in solver.inp
c
      if (timestart .lt. 0.0) then

        write (temp1,"('(''TimeStamp@'',i',i1,',A1)')")
     &         itmp
        write (fname1,temp1) (myrank+1),'?'
        fname1 = trim(fname1)
        intfromfile=0
        call readheader(descriptor,fname1 // char(0),intfromfile,
     &                  ione,'integer' // char(0),iotype)

        if(intfromfile(1).ne.0) then 

          call readdatablock(descriptor,fname1 // char(0),time,ione,
     &                         'double' // char(0),iotype)

        else
          if (myrank.eq.master) then
            warning='Time is set to zero (SAFE)'
            write(*,*) warning
          end if
          time = 0.0
        endif
      else 
        open(unit=72,file='numstart.dat',status='old')
        nl_numstart = 0
        do
           read(72,*,iostat=readCode)
           if(readCode.ne.0) exit
           nl_numstart = nl_numstart + 1
        enddo
        close(72)
        if(nl_numstart.eq.1) then
           time = timestart
           if(myrank.eq.master) write(*,*)
     &     'Put simulation time in numstart.dat!'
        else
           open(unit=72,file='numstart.dat',status='old')
           read(72,*)
           read(72,*) timestart
           time = timestart
        endif
        close(72)
      endif

      itmp = int(log10(float(myrank+1)))+1
      write (temp1,"('(''solution@'',i',i1,',A1)')") itmp
      write (fname1,temp1) (myrank+1),'?'
      fname1 = trim(fname1)

      intfromfile=0
!      if(myrank.eq.master) write(*,*) "calling readheader"
      call readheader(descriptor,fname1 // char(0) ,intfromfile,
     & ithree,'integer' // char(0),iotype)
!      if (mrank.eq.master) write(*,*) "called readheader"
c
c.... read the values of primitive variables into q
c
      allocate( qold(nshg,ndof) )

      if(intfromfile(1).ne.0) then ! solution field was found
        nshg2=intfromfile(1)
        ndof2=intfromfile(2)
        lstep=intfromfile(3)
        if(ndof2.ne.ndof) then
         warning='WARNING more data in restart than needed: keeping 1st'
         write(*,*) warning , ndof, ndof2, myrank, fnamer 
        endif
c
        if (nshg2 .ne. nshg)
     &        call error ('restar  ', 'nshg   ', nshg)
         allocate( qread(nshg,ndof2) )
         iqsiz=nshg*ndof2
!         if (myrank.eq.master) write(*,*) 'calling readdatablock'
         call readdatablock(descriptor,fname1 // char(0),qread,iqsiz,
     &                         'double' // char(0),iotype)
!        if (myrank.eq.master) write(*,*) 'called readdatablock'
         qold(:,1:ndof)=qread(:,1:ndof)
         deallocate(qread)
      else
         if (myrank.eq.master) then
            if (matflg(1,1).eq.0) then ! compressible
               warning='Solution is set to zero (with p and T to one)'
            else
               warning='Solution is set to zero'
            endif
            write(*,*) warning
         endif
         qold=zero
         if (matflg(1,1).eq.0) then ! compressible
            qold(:,1)=one ! avoid zero pressure
            qold(:,nflow)=one ! avoid zero temperature
         endif
      endif

      write (temp1,"('(''time derivative of solution@'',i',i1,',A1)')")
     &       itmp
      write (fname1,temp1) (myrank+1),'?'
      fname1 = trim(fname1)
      intfromfile=0
      call readheader(descriptor,fname1 // char(0) ,intfromfile,
     & ithree,'integer' // char(0),iotype)
      allocate( acold(nshg,ndof) )
      if(intfromfile(1).ne.0) then
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         lstep=intfromfile(3)

         if (nshg2 .ne. nshg)
     &        call error ('restar  ', 'nshg   ', nshg)

         allocate( acread(nshg,ndof2) )
         acread=zero
         iacsiz=nshg*ndof2
         call readdatablock(descriptor,fname1 // char(0),acread,
     &    iacsiz, 'double' // char(0),iotype)
         acold(:,1:ndof)=acread(:,1:ndof)
         deallocate(acread)
      else
         if (myrank.eq.master) then
            warning='Time derivative of solution is set to zero (SAFE)'
            write(*,*) warning
         endif
         acold=zero
      endif

c      call creadlist(irstin,ithree,nshg2,ndisp,lstep)
      if (ideformwall.eq.1) then
          write (temp1,"('(''displacement@'',i',i1,',A1)')")
     &           itmp
          write (fname1,temp1) (myrank+1),'?'
          fname1 = trim(fname1)
          intfromfile=0
          call readheader(descriptor,fname1 // char(0),intfromfile,
     &     ithree,'integer' // char(0),iotype)

         nshg2=intfromfile(1)
         ndisp=intfromfile(2)
         lstep=intfromfile(3)
         if(ndisp.ne.nsd) then
            warning='WARNING ndisp not equal nsd'
            write(*,*) warning , ndisp
         endif
c
         if (nshg2 .ne. nshg) 
     &        call error ('restar  ', 'nshg   ', nshg)
c
c.... read the values of primitive variables into uold
c
         allocate( uold(nshg,nsd) )
         allocate( uread(nshg,nsd) )
         
         iusiz=nshg*nsd
         call readdatablock(descriptor,fname1 // char(0) ,uread,iusiz,
     &                   'double' // char(0),iotype)
         uold(:,1:nsd)=uread(:,1:nsd)
       else
         allocate( uold(nshg,nsd) )
         uold(:,1:nsd) = zero
       endif

c
c.... close c-binary files
c
      call closefile( descriptor, "read" // char(0) )
      call finalizephmpiio( descriptor )
c
      deallocate(xread)
      if ( numpbc > 0 )  then
         deallocate(bcinpread)
         deallocate(ibctmpread)
      endif
      deallocate(iperread)
      if(numpe.gt.1)
     &     deallocate(ilworkread)
      deallocate(nbcread)

      return
c
 994  call error ('input   ','opening ', igeom)
 995  call error ('input   ','opening ', igeom)
 997  call error ('input   ','end file', igeom)
 998  call error ('input   ','end file', igeom)
c
      end

c
c No longer called but kept around in case....
c
      subroutine genpzero(iBC)

      use pointer_data
c
      include "common.h"
      integer iBC(nshg)
c
c....  check to see if any of the nodes have a dirichlet pressure
c
      pzero=1
      if (any(btest(iBC,2))) pzero=0  
c
      do iblk = 1, nelblb
         npro = lcblkb(1,iblk+1)-lcblkb(1,iblk)
         do i=1, npro
            iBCB1=miBCB(iblk)%p(i,1)
c     
c.... check to see if any of the nodes have a Neumann pressure 
c     but not periodic (note that 
c     
            if(btest(iBCB1,1)) pzero=0
         enddo
c     
c.... share results with other processors
c     
         pzl=pzero
         if (numpe .gt. 1)
     &        call MPI_ALLREDUCE (pzl, pzero, 1,
     &        MPI_DOUBLE_PRECISION,MPI_MIN, MPI_COMM_WORLD,ierr)
           
      enddo
c
c.... return
c
      return
c
      end

