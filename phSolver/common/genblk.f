        subroutine genblk (IBKSZ)
c
c----------------------------------------------------------------------
c
c  This routine reads the interior elements and generates the
c  appropriate blocks.
c
c Zdenek Johan, Fall 1991.
c  Parallel IO implemented with the capacity to read different number of 
c  topologies per part
c
c Michel Rasquin, 2012
c----------------------------------------------------------------------
c
        use pointer_data
c
        include "common.h"
        include "mpif.h" !Required to determine the max for itpblk
c
        integer, allocatable :: ientp(:,:)
        integer mater(ibksz)
        integer intfromfile(50) ! integers read from headers
        character*255 fname1

cccccccccccccc New Phasta IO starts here ccccccccccccccccccccccccc

        integer :: descriptor, descriptorG, GPID, color, nfiles
        integer ::  numparts, writeLock
        integer :: ierr_io, numprocs, itmp, itmp2
!MR CHANGE
        integer :: itpblktot,ierr,iseven
!MR CHANGE END
        character*255 fnamer, fname2, temp2
        character*64 temp1, temp3
!THIS NEEDS TO BE CLEANED - MR
        nfiles = nsynciofiles
!        nfields = nsynciofieldsreadgeombc
        numparts = numpe !This is the common settings. Beware if you try to compute several parts per process

!        nppp = numparts/numpe
!        nppf = numparts/nfiles

        color = int(myrank/(numparts/nfiles)) !Should call the SyncIO routine here
        itmp2 = int(log10(float(color+1)))+1
        write (temp2,"('(''geombc-dat.'',i',i1,')')") itmp2
        temp2=trim(temp2)
        write (fnamer,temp2) (color+1)
        fnamer=trim(fnamer)

        ione=1
        itwo=2
        iseven=7
        ieleven=11
        itmp = int(log10(float(myrank+1)))+1

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
        iel=1
        itpblk=nelblk
!MR CHANGE

        ! Get the total number of different interior topologies in the whole domain. 
        ! Try to read from a field. If the field does not exist, scan the geombc file.
        itpblktot=-1
        write(temp1,
     &   "('(''total number of interior tpblocks@'',i',i1,',A1)')") itmp

        write (fname2,temp1) (myrank+1),'?'
        call readheader(igeom,fname2 // char(0) ,itpblktot,ione,
     &  'integer' // char(0),iotype) 

!        write (*,*) 'Rank: ',myrank,' interior itpblktot intermediate:',
!     &               itpblktot

        if (itpblktot == -1) then 
          ! The field 'total number of different interior tpblocks' was not found in the geombc file.
          ! Scan all the geombc file for the 'connectivity interior' fields to get this information.
          iblk=0
          neltp=0
          do while(neltp .ne. -1) 

            ! intfromfile is reinitialized to -1 every time.
            ! If connectivity interior@xxx is not found, then 
            ! readheader will return intfromfile unchanged

            intfromfile(:)=-1
            iblk = iblk+1
            write (temp1,"('connectivity interior',i1)") iblk
            temp1 = trim(temp1)
            write (temp3,"('(''@'',i',i1,',A1)')") itmp
            write (fname2, temp3) (myrank+1), '?'
            fname2 = trim(temp1)//trim(fname2)

            !write(*,*) 'rank, fname2',myrank, trim(adjustl(fname2))
            call readheader(igeom,fname2 // char(0),intfromfile,
     &       iseven,'integer' // char(0),iotype)
            neltp = intfromfile(1) ! -1 if fname2 was not found, >=0 otherwise
          end do
          itpblktot = iblk-1   
        end if

        if (myrank == 0) then
          write(*,*) 'Number of interior topologies: ',itpblktot
        endif
!        write (*,*) 'Rank: ',myrank,' interior itpblktot final:',
!     &               itpblktot

!MR CHANGE END

        nelblk=0
        mattyp = 0
        ndofl = ndof
        nsymdl = nsymdf

        do iblk = 1, itpblktot
           writeLock=0;
!MR CHANGE END
c
c           read(igeomBAK) neltp,nenl,ipordl,nshl, ijunk, ijunk, lcsyst
c           call creadlist(igeomBAK,iseven,
c     &          neltp,nenl,ipordl,nshl, ijunk, ijunk, lcsyst)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           write (temp1,"('connectivity interior',i1)") iblk
           temp1=trim(temp1)
           write (temp3,"('(''@'',i',i1,',A1)')") itmp
           write (fname2, temp3) (myrank+1), '?'
           fname2 = trim(temp1)//trim(fname2)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c           fname1='connectivity interior?'

           ! Synchronization for performance monitoring, as some parts do not include some topologies
           call MPI_Barrier(MPI_COMM_WORLD,ierr) 
           call readheader(igeom,fname2 // char(0) ,intfromfile,
     &     iseven,"integer" // char(0), iotype)
           neltp  =intfromfile(1)
           nenl   =intfromfile(2)
           ipordl =intfromfile(3)
           nshl   =intfromfile(4)
           ijunk  =intfromfile(5)
           ijunk  =intfromfile(6)
           lcsyst =intfromfile(7)
           allocate (ientp(neltp,nshl))
c           read(igeomBAK) ientp
           iientpsiz=neltp*nshl

           if (neltp==0) then
              writeLock=1;
           endif

           call readdatablock(igeom,fname2 // char(0),ientp,iientpsiz,
     &                     "integer" // char(0), iotype)

!            call closefile( igeom, "read" // char(0) )
!            call finalizephmpiio( igeom )

!MR CHANGE
           if(writeLock==0) then
!MR CHANGE

             do n=1,neltp,ibksz 
                nelblk=nelblk+1
                npro= min(IBKSZ, neltp - n + 1)
c
                lcblk(1,nelblk)  = iel
c                lcblk(2,nelblk)  = iopen ! available for later use
                lcblk(3,nelblk)  = lcsyst
                lcblk(4,nelblk)  = ipordl
                lcblk(5,nelblk)  = nenl
                lcblk(6,nelblk)  = nfacel
                lcblk(7,nelblk)  = mattyp
                lcblk(8,nelblk)  = ndofl
                lcblk(9,nelblk)  = nsymdl 
                lcblk(10,nelblk) = nshl ! # of shape functions per elt
c
c.... allocate memory for stack arrays
c
                allocate (mmat(nelblk)%p(npro))
c
                allocate (mien(nelblk)%p(npro,nshl))
                allocate (mxmudmi(nelblk)%p(npro,maxsh))
c
c.... save the element block
c
                n1=n
                n2=n+npro-1
                mater=1   ! all one material for now
                call gensav (ientp(n1:n2,1:nshl),
     &                       mater,           mien(nelblk)%p,
     &                       mmat(nelblk)%p)
                iel=iel+npro
c
             enddo
!MR CHANGE
           endif
!MR CHANGE
           deallocate(ientp)
        enddo
        lcblk(1,nelblk+1) = iel
c
c.... return
c
CAD        call timer ('Back    ')
c
        return
c
1000    format(a80,//,
     &  ' N o d a l   C o n n e c t i v i t y',//,
     &  '   Elem  ',/,
     &  '  Number  ',7x,27('Node',i2,:,2x))
1100    format(2x,i5,6x,27i8)
        end
