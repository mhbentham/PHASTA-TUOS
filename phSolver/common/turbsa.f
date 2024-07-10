c-----------------------------------------------------------------------
c
c     Spalart-Allmaras turbulence model constants 
c
c-----------------------------------------------------------------------
      module turbSA

      real*8, allocatable ::  d2wall(:), x2wall(:), y2wall(:), z2wall(:)
      real*8, allocatable ::  d2wal(:,:)
      real*8, allocatable ::  wnrm(:,:)
      integer, allocatable :: otwn(:)
      real*8, allocatable :: otwn_vec(:,:)
      real*8, allocatable :: effvisc(:)
      real*8  saCb1, saCb2, saCw1, saCw2, saCw3, saCv1, saSigma,
     &        saKappa, saKappaP2Inv, saCv2P3, saCw3P6, saSigmaInv
      integer, allocatable :: sidmapg(:) ! list of all surfID's, low to high
      
      parameter ( 
     &     saCb1        = 0.1355d0,
     &     saCb2        = 0.622d0,
     &     saCw1        = 3.239067817d0,
     &     saCw2        = 0.3d0,
     &     saCw3        = 2.0d0,
     &     saKappa      = 0.41d0,
     &     saSigma      = 0.666666666666666667d0,
     &     saCv1        = 7.1d0,
     &     saKappaP2Inv = 5.94883997620464d0,
     &     saCv1P3      = 3.579109999999999d+02,
     &     saCw3P6      = 64.0d0,
     &     saSigmaInv   = 1.50d0
     &     )

      
      end module

c-----------------------------------------------------------------------
c
c     Initialize: compute the distance to the nearest wall for 
c     each vertex in the mesh.
c
c     Michael Yaworski (fall 1998)
c
c-----------------------------------------------------------------------
      subroutine initTurb( x )
      
      use     pointer_data
      use     turbSA
      include "common.h"
      include "mpif.h"

      character*20 fname1,  fmt1
      real*8   x(numnp,nsd)
      integer  nwall(numpe),      idisp(numpe)
      real*8, allocatable :: xwi(:,:,:), xw(:,:,:)

!MR CHANGE
      integer :: ierr, numVars
      integer :: ifoundd2wall
!MR CHANGE END

      numVars = 4 ! Number of d2wall variables per node (magnitude + 3 components in x, y, z
      allocate ( d2wal(numnp,numvars) )
      allocate ( d2wall(numnp) )
      allocate ( x2wall(numnp) )
      allocate ( y2wall(numnp) )
      allocate ( z2wall(numnp) )

!MR CHANGE
      if(myrank.eq.master) then
        write (*,*) 'entering initTurb'
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      ! First we try to read the d2wall field from either the restart files or the d2wall files
      call read_d2wall(myrank,numnp,numvars,d2wal,ifoundd2wall)
!MR CHANGE END

      if(ifoundd2wall.eq.1) then   ! d2wall was not found so calculate the distance

         do i=1,numnp
           d2wall(i) = d2wal(i,1) 
           x2wall(i) = d2wal(i,2)
           y2wall(i) = d2wal(i,3)
           z2wall(i) = d2wal(i,4)
         enddo

      elseif(ifoundd2wall.eq.0) then   ! d2wall was not found so calculate the distance
        if(myrank.eq.master) then
          write (*,*) 'Computing the d2wall field'
        endif
c
c   hard code trip point until we are sure it is worth doing
c

c
c.... Count the welts (wall-elements)
c
         nwalli=0
         do iblk = 1, nelblb    ! loop over boundary elt blocks
            npro = lcblkb(1,iblk+1) - lcblkb(1,iblk)
            do j = 1, npro
               if(btest(miBCB(iblk)%p(j,1),4)) nwalli=nwalli+1
            enddo
         enddo
c
c.... Create wallnode-coord list for welts on processor
c
         if (nwalli.eq.0) nwalli=1 !  patch for mpi's lack of imagination
         allocate (xwi(nwalli,nenb+1,nsd))
         xwi = 1.0d32
         xwi(:,nenb+1,:)=zero
         nwalli = 0
         do iblk = 1, nelblb    ! loop over boundary elt blocks
c
            iel    = lcblkb(1,iblk)
            nenbl  = lcblkb(6,iblk) ! no. of vertices per bdry. face
            npro   = lcblkb(1,iblk+1) - iel 
c
            do j = 1, npro      ! loop over belts in this blk
               if(btest(miBCB(iblk)%p(j,1),4)) then
                  nwalli=nwalli+1
c assemble local coord list
                  do node = 1, nenbl
                     xwi(nwalli,node,1:3)=x(mienb(iblk)%p(j,node),:)
                  enddo
c put the centroid coordinates in the last slot
                  do node = 1, nenbl
                     xwi(nwalli,nenb+1,:)=xwi(nwalli,nenb+1,:)
     &                    +xwi(nwalli,node,:)
                  enddo
                  xwi(nwalli,nenb+1,:)=xwi(nwalli,nenb+1,:)/nenbl
c
               endif
            enddo               ! loop over belts in this blk
c
         enddo                  ! loop over boundary elt blocks
c
         if (nwalli.eq.0) xwi=1.0e32 ! fix for mpi's lack of imagination
         if (nwalli.eq.0) nwalli=1 !  patch for mpi's lack of imagination
c
c  Pool "number of welts" info from all processors
c
cMPI_ALLGATHER(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,comm) 
c[ IN sendbuf] starting address of send buffer (choice) 
c[ IN sendcount] number of elements in send buffer (integer) 
c[ IN sendtype] data type of send buffer elements (handle) 
c[ OUT recvbuf] address of receive buffer (choice) 
c[ IN recvcount] number of elements received from any process (integer) 
c[ IN recvtype] data type of receive buffer elements (handle) 
c[ IN comm] communicator (handle)
         if (numpe.gt.1) then   ! multiple processors
c write the number of wall elts on the jth processor to slot j of nwall
            call MPI_ALLGATHER(nwalli,1,MPI_INTEGER,nwall,1,
     &          MPI_INTEGER,MPI_COMM_WORLD,ierr)
c count up the total number of wall elts among all processes
            nwallt=0
            do j=1,numpe
               nwallt=nwallt+nwall(j)
            enddo
         else                   ! single processor
c the local information is the global information for single-processor
            nwall=nwalli
            nwallt=nwalli
         endif                  ! if-else for multiple processors
c
c  Make all-processor wallnode-coord collage
c
        if(myrank.eq.master)write(*,*)'wall elem count = ',nwallt
         allocate (xw(nwallt,nenb+1,nsd))
         if (numpe.gt.1) then   ! multiple processors
c we will gather coordinates from local on-proc sets to a global set
c we will stack each processor's coordinate list atop that of the
c previous processor.  If the coordinate list for processor i is
c called xwi, then our global coordinate list xw will look like this:
c ---------------------------------------------------------------
c | xw1            | xw2                | xw3        |   ...    |
c ---------------------------------------------------------------
c  <---nwall(1)---> <-----nwall(2)-----> <-nwall(3)->
c  <------------------------nwallt-----------------------...---->
c To accomplish this with MPI, we use MPI_ALLGATHERV, summarized as:
cMPI_ALLGATHERV(sendbuf,sendcount,sendtype,recvbuf,recvcount,disp,recvtype,comm) 
c[ IN sendbuf] starting address of send buffer (choice) 
c[ IN sendcount] number of elements in send buffer (integer) 
c[ IN sendtype] data type of send buffer elements (handle) 
c[ OUT recvbuf] address of receive buffer (choice) 
c[ IN recvcount] number of elements received from any process (integer) 
c[ IN disp] displacement array
c[ IN recvtype] data type of receive buffer elements (handle) 
c[ IN comm] communicator (handle)
c The displacement array disp is an array of integers in which the jth
c entry indicates which slot of xw marks beginning of xwj
c So, first we will build this displacement array
            idisp(:)=0 ! starting with zero, since MPI likes C-numbering
            do j=2,numpe
               idisp(j)=idisp(j-1)+nwall(j-1) ! see diagram above
            enddo
c Now, we gather the data one slice at a time (1:nwalli)
            do j=1,nenb+1
               do k=1,nsd
                  call MPI_ALLGATHERV(xwi(:,j,k),nwalli,
     &                 MPI_DOUBLE_PRECISION,xw(:,j,k),nwall,idisp,
     &                 MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
               enddo
            enddo
         else                   ! single-processor
c global data is local data in single processor case
            xw=xwi
         endif

c
c  For each point, loop over wall nodes and calculate distance; 
c  save the distance in this node's position of d2wall if it's 
c  shorter than currently stored distance
c
         d2wall=1.0e32
         do i=1,numnp
!	  if (5000*int(i/5000).eq.i) write(*,*) 'i, numnp = ', i, numnp
            do j=1, nwallt
               do k=1,nenb+1
                  distance =  ( x(i,1) - xw(j,k,1) )**2
     &                 +( x(i,2) - xw(j,k,2) )**2
     &                 +( x(i,3) - xw(j,k,3) )**2
                  if ( d2wall(i).gt.distance ) then
                    d2wall(i) = distance
                    x2wall(i) = xw(j,k,1) - x(i,1)   ! Vector direction is TOWARDS the wall
                    y2wall(i) = xw(j,k,2) - x(i,2)
                    z2wall(i) = xw(j,k,3) - x(i,3)
                  endif
               enddo
            enddo
         enddo
         d2wall=sqrt(d2wall)
! Jun, Nov. 2014:
         do i=1,numnp
           d2wal(i,1) = d2wall(i)
           d2wal(i,2) = x2wall(i)
           d2wal(i,3) = y2wall(i)
           d2wal(i,4) = z2wall(i)
         enddo
c
c
         deallocate(xwi)
         deallocate(xw)
c
c.... write d2wall to a file so we don't have to do this again
c

          call write_d2wall(myrank,numnp,numvars,d2wal) !See new_interface.c

      endif

!MR CHANGE
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myrank.eq.master) then
        write (*,*) 'leaving initTurb'
      endif
!MR CHANGE

      return
995     call error ('d2wall  ','opening ', 72)
996     call error ('d2wall  ','opening ', 72)

      end subroutine
