c-----------------------------------------------------------------------
c
c  This module conveys temporal BC data.  Below functions read in the data
c  and interpolate it to the current time level. 
c
c-----------------------------------------------------------------------
      module specialBC

      real*8, allocatable ::  BCt(:,:,:), acs(:,:), spamp(:)
      real*8, allocatable ::  ytarget(:,:)
      integer, allocatable :: nBCt(:), numBCt(:)
     

      integer ntv,nptsmax
c$$$      integer itvn
      end module

c-----------------------------------------------------------------------
c
c  This module conveys flow rate history for the different impedance outlets
c  over one period. Below functions read in the data and store it for the
c  current time level. 
c
c-----------------------------------------------------------------------
      module convolImpFlow

      real*8, allocatable ::  QHistImp(:,:), ValueImpt(:,:,:)
      real*8, allocatable ::  ValueListImp(:,:), ConvCoef(:,:)
      real*8, allocatable ::  ImpConvCoef(:,:), pold(:)
      integer ntimeptpT,numDataImp
      integer, allocatable :: nImpt(:), numImpt(:)
      integer nptsImpmax


      end module

c-----------------------------------------------------------------------
c
c     Initialize:
c
c-----------------------------------------------------------------------
      subroutine initSponge( y,x)
      
      use     specialBC
      include "common.h"
      
      real*8   y(nshg,nflow), x(numnp,3)
      allocate (ytarget(nshg,nflow))  
      
      if(matflg(5,1).eq.5) then
         write(*,*) 'calculating IC sponge'
         ytarget = y
      else
         write(*,*) 'calculating Analytic sponge'

c
c OLD STyle sponge pushed onto target.  You need to be sure that your
c solver.inp entries for start and stop of sponge match as well as the
c growth rates
c
      vcl=datmat(1,5,1)         ! velocity on centerline
      rslc=datmat(2,5,1)        ! shear layer center radius
      bfz=datmat(3,5,1)
      we=3.0*29./682.
      rsteep=3.0
      zstart=30.0
      radst=10.0
      radsts=radst*radst
      do id=1,numnp
         radsqr=x(id,2)**2+x(id,1)**2
c         if((x(id,3).gt. zstart) .or. (radsqr.gt.radsts))  then
            rad=sqrt(radsqr)
            radc=max(rad,radst)
            zval=max(x(id,3),zstart)
            utarget=(tanh(rsteep*(rslc-rad))+one)/two*
     &                    (vcl-we) + we
            Ttarget  = press/(ro*Rgas)
            ptarget= press
            ytarget(id,1) = zero
            ytarget(id,2) = zero
            ytarget(id,3) = utarget
            ytarget(id,4) = ptarget
            ytarget(id,5) = Ttarget            
c         endif
      enddo
      endif
      return
      end


c-----------------------------------------------------------------------
c
c     Initialize:time varying boundary condition
c
c-----------------------------------------------------------------------
      subroutine initBCt( x, iBC, BC )
      
      use     specialBC
      include "common.h"
      
      real*8   x(numnp,nsd), BC(nshg,ndofBC), rj1,rj2,rj3,rj4,distd,epsd
      real*8 rj_temp
      allocatable :: rj_temp(:)
      integer  iBC(numnp)
      character*80 card
      real*8 distds
      real*8 dd
c
c  This one should be used for boundary layer meshes where bct.dat must
c  be given to greater precision than is currently being generated.
c
c      epsd=1.0d-12    ! this is distance SQUARED to save square root
      epsd=1.0d-9            ! this is distance SQUARED to save square root

      ic=0                      !count the number on this processor
     
      if(any(ibits(iBC,3,3).eq.7)) then
         if(myrank.eq.master) write(*,*) 'opening bct.dat'
c         open(unit=567, file='bct.dat',status='old')
         open(unit=567, file='bct.dat',ACTION='READ',STATUS='old')
c         read (567,'(a80)') card
c reading the #of nodal points - ntv, and #of time series  - ntpts.
c           read (card,*) ntv, nptsmax
        read(567,*) ntv,nptsmax
         allocate (nBCt(numnp))  
         allocate (numBCt(ntv))  
        if (tvbcswitch.eq.0) then
         allocate (BCt(ntv,nptsmax,4))  
	else 
	  allocate (BCt(ntv,nptsmax,ndof+1)) 
          allocate (rj_temp(ndof+1)) !* just for replace rj1,rj2..Azat 10/18/04

	endif
	if (myrank.eq.master) write(*,*) 'Loop over nodes'
        do k=1,ntv !* this is loop over the boundary nodes (global)
            read(567,*) x1,x2,x3,ntpts
c
c Find the point on the boundary (if it is on this processor)
c that matches this point
c
            do i=1,numnp
               if(ibits(ibc(i),3,3) .eq.7) then
                  dd= distds(x1,x2,x3,x(i,1),x(i,2),x(i,3))
                  if(dd.lt.epsd) then
                     ic=ic+1
                     nBCt(ic)=i ! the pointer to this point
                     numBCt(ic)=ntpts ! the number of time series
                     do j=1,ntpts
  	               if (tvbcswitch.eq.0) then
                        read(567,*) (BCt(ic,j,n),n=1,4)
	               else
c************Now the structure of the BCT.dat is changed: previuos structure was : u v w t
c***********it is replaced by follows: t p u v w sc1 sc2 ..
	                read(567,*) (BCt(ic,j,n),n=1,ndof+1)
	               endif
                     enddo
                     exit
                  endif
               endif
            enddo
            if(i.eq.numnp+1) then
c
c  if we get here the point was not found.  It must be on another
c  processor so we read past this record and move on
c
               do j=1,ntpts
              if (tvbcswitch.eq.0) then
                  read(567,*) rj1,rj2,rj3,rj4
	      else
                read(567,*) (rj_temp(jj),jj=1,ndof+1)
	      endif
               enddo
            endif
         enddo                  ! end of the loop over ntv
         BCt(:,:,4)=BCt(:,:,4)*bcttimescale
      endif                     ! any 3 component nodes
      itvn=ic
      close(567)
      if (ic.gt.0) then
         write(*,*)'myrank=',myrank,' and I found ',ic,' nodes.'
c      else
c         deallocate(nBCt)
c         deallocate(numBCt)
c         deallocate(BCt)
      endif

      return
      end


      subroutine BCint(timel,shp,shgl,shpb,shglb,x,BC,iBC)

      use     specialBC ! brings in itvn,nbct, bct, numbct, nptsmax

      include "common.h"

      real*8   BC(nshg,ndofBC), timel,t
      real*8   x(numnp,nsd),   
     &         shp(MAXTOP,maxsh,MAXQPT),
     &         shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &         shpb(MAXTOP,maxsh,MAXQPT),
     &         shglb(MAXTOP,nsd,maxsh,MAXQPT)

      integer  iBC(numnp),nlast,i,j,nper 

      do i =1,itvn ! itvn is the number of varying nodes on this proc 

         nlast=numBCt(i)     ! number of time series to interpolate from
	if (tvbcswitch.eq.0) then !it is for keeping previous structure
         nper=timel/BCt(i,nlast,4)! number of periods completed to shift off


         t=timel-nper*BCt(i,nlast,4)  ! now time in periodic domain

         do j=2,nlast   !loop to find the interval that we are in

            if(BCt(i,j,4).gt.t) then  ! this is upper bound, j-1 is lower

               wr=(t-BCt(i,j-1,4))/(BCt(i,j,4)-BCt(i,j-1,4))
               BC(nbct(i),3:5)= BCt(i,j-1,1:3)*(one-wr) 
     &                        + BCt(i,j,1:3)*wr
               exit

            endif
	  enddo !* end of loop over nlast
	else
          nper=timel/BCt(i,nlast,1) !* Now the structure is changed. see initBCt
          t=timel-nper*BCt(i,nlast,1)  ! now time in periodic domain
	  do j=2,nlast   !loop to find the interval that we are in   
            if(BCt(i,j,1).gt.t) then  ! this is upper bound, j-4 is lower
              wr=(t-BCt(i,j-1,1))/(BCt(i,j,1)-BCt(i,j-1,1))
              BC(nbct(i),1:5)= BCt(i,j-1,2:6)*(one-wr)
     &                       + BCt(i,j,2:6)*wr
              BC(nbct(i),7:8)= BCt(i,j-1,7:8)*(one-wr) !* I am not sure that is correct, because I do not understand in itrBCsclr the meaning of the BC(:,6) 
     &                       + BCt(i,j,7:8)*wr
              exit
            endif !* for j-1
	  enddo !* end of loop over nlast 
	endif !* for switch between old structure and new
      enddo !* end of loop over itvn
      return
      end

      function distds(x1,y1,z1,x2,y2,z2)
      real*8 distds 
      real*8 x1,y1,z1,x2,y2,z2,x,y,z
      x=x1-x2
      y=y1-y2
      z=z1-z2
      distds=x*x+y*y+z*z
      return
      end
c-----------------------------------------------------------------------
c   initialize the impedance boundary condition:
c   read the data in initImpt
c   interpolate the data to match the process time step in Impint
c-----------------------------------------------------------------------
      subroutine initImpt()
      
      use convolImpFlow
      include "common.h"

      open(unit=817, file='impt.dat',status='old')
         read (817,*) nptsImpmax
         allocate (numImpt(numImpSrfs))  
         allocate (ValueImpt(nptsImpmax,2,numImpSrfs))
         ValueImpt=0
         do k=1,numImpSrfs
            read (817,*) numDataImp
            numImpt(k) = numDataImp
            do j=1,numDataImp
               read(817,*) (ValueImpt(j,n,k),n=1,2) ! n=1 time, 2 value
            enddo
         enddo
      close(817)
      
      allocate (ValueListImp(ntimeptpT+1,numImpSrfs))
      ValueListImp = zero 
      ValueListImp(ntimeptpT+1,:) = ValueImpt(1,2,:) !Z(time=0), last entry
      ValueListImp(1,:) = ValueImpt(1,2,:) !Z(time=0)=Z(time=T)
      return
      end
      
      
      
      subroutine Impint(ctime,jstep)
      
      use convolImpFlow
      include "common.h"
      
      real*8 ctime, ptime
      integer nlast, nper, k, j , jstep
      
         
      do k =1,numImpSrfs
         nlast=numImpt(k)     ! number of time series to interpolate from
         nper=ctime/ValueImpt(nlast,1,k)!number of periods completed to shift off
         ptime = ctime-nper*ValueImpt(nlast,1,k)  ! now time in periodic domain
            
         do j=2,nlast   !loop to find the interval that we are in

            if(ValueImpt(j,1,k).gt.ptime) then  ! this is upper bound, j-1 is lower
               wr=(ptime-ValueImpt(j-1,1,k))
     &             / ( ValueImpt(j,1,k)-ValueImpt(j-1,1,k) )
               ValueListImp(jstep,k)= ValueImpt(j-1,2,k)*(one-wr) 
     &                        + ValueImpt(j,2,k)*wr
               exit
            endif

         enddo
      enddo
      return
      end
      
      
c----------------------------------------------------------------------- 
c     returns in pold the history dependent part of the pressure in the
c     impedance/flow rate convolution for the impedance BC      
c-----------------------------------------------------------------------      
      subroutine pHist(pressHist,QHist,betas,nTimePoint,nSrfs)

      include "common.h"
      
      integer  nTimePoint,nSrfs
      real*8   pressHist(0:MAXSURF)
      real*8   QHist(nTimePoint+1,nSrfs),betas(nTimePoint+2,nSrfs)!don't need here betas(ntimePoint+2)
      !but pb of array passing if cut at nTimePoint+1
      pressHist=zero
      do k=1,nSrfs
        do j=1,nTimePoint+1
            pressHist(k) = pressHist(k) + QHist(j,k)*betas(j,k)
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
c
c     Initialize arrays for time varying inlet boundary condition
c	where nodal values are set by user function (not bct.dat)
c
c-----------------------------------------------------------------------
      subroutine initsetBC(iBC)

      use     specialBC ! brings in itvn, nbct

      include "common.h"

      integer  iBC(numnp),i,j
c
c Only set nodal values for prescribed boundary conditions
c
      if(any(ibits(iBC,3,3).eq.7)) then
        ntv=0                      !count the number on this processor
        do i=1,numnp
          if(ibits(ibc(i),3,3) .eq.7) ntv = ntv+1
        enddo
      endif
c
c Allocate arrays and store node numbers
c
      allocate (nBCt(ntv)) 
      ic = 0 
      do i=1,numnp
        if(ibits(ibc(i),3,3) .eq.7) then
          ic=ic+1
          nBCt(ic)=i ! the pointer to this point
        endif
      enddo
c
      itvbc = ic
c
      return
      end


c-----------------------------------------------------------------------
c
c     Initialize arrays for time varying inlet boundary condition
c	where nodal values are set by user function (not bct.dat)
c
c-----------------------------------------------------------------------
      subroutine setBC(time1,x,iBC,BC)

      use     specialBC ! brings in itvn,nbct, bct, numbct, nptsmax

      include "common.h"

      real*8   BC(nshg,ndofBC), time1, t
      real*8   x(numnp,nsd), rx(nsd) 
      real*8   rvel(nsd), rpress, rsclr  
      integer  iBC(numnp),nlast,i,j,nper 
c
      do i =1,itvbc ! itvbc is the number of varying nodes on this proc
        inode = nBCt(i)
        rx = x(inode,:)
c set velocity
        if (ibits(ibc(inode),3,3) .eq.7) then
          call usrbc('velocity',rx,time1,rvel,nsd)
          BC(inode,3:5)= rvel
c          write(*,1001) istep,time1,inode,rx,rvel
 1001 format(" DEBUG: velocity ",i8,2x,e12.6,i8,6(2x,e12.6))
        endif
c set pressure
        if (ibits(ibc(inode),2,1) .eq.1) then
          call usrbc('pressure',rx,time1,rpress,1)
          BC(inode,1)= rpress
c          write(*,1004) istep,time1,inode,rx,rpress
 1004 format(" DEBUG: pressure ",i8,2x,e12.6,i8,4(2x,e12.6))
        endif
c set temp
        if (ibits(ibc(inode),1,1) .eq.1) then
          call usrbc('  temper',rx,time1,rtemp,1)
          BC(inode,2)= rtemp
c          write(*,1005) istep,time1,inode,rx,rpress
 1005 format(" DEBUG: temperature ",i8,2x,e12.6,i8,4(2x,e12.6))
        endif
c set scalar 1
        if (ibits(ibc(inode),6,1) .eq.1) then
          call usrbc('  scalar',rx,time1,rsclr,1)
          BC(inode,7)= rsclr
c          write(*,1000) istep,time1,inode,rx,rsclr 
 1000 format(" DEBUG: scalar1 ",i8,2x,e12.6,i8,4(2x,e12.6))
        endif
c set scalar 2
        if (ibits(ibc(inode),7,1) .eq.1) then
          call usrbc('  scalar',rx,time1,rsclr,1)
          BC(inode,8)= rsclr
c          write(*,1003) istep,time1,inode,rx,rsclr
 1003 format(" DEBUG: scalar2 ",i8,2x,e12.6,i8,4(2x,e12.6))
        endif
      enddo
c
      return
      end



c-----------------------------------------------------------------------
c
c     User function that defines nodal values
c
c-----------------------------------------------------------------------
      subroutine usrbc(cflag,rx,time1,xval,nsize)
      include "common.h"
c
      character*8 cflag
      integer nsize
      real*8   time1, rx(nsd), xval(nsize)
c
c Velocity
c
      if (cflag .eq. "velocity") then
         if (rx(2) .gt. 0.0) then
           xval(1) = 1.0
         else
           xval(1) = 0.9
         endif
         xval(2) = 0.0
         xval(3) = 0.0
      endif
c
c Pressure
c
      if (cflag .eq. "pressure") then
         xval(1) = 0.0
      endif
c
c Temperature
c
      if (cflag .eq. "  temper") then
         xval(1) = 0.0
      endif
c
c Scalar
c
      if (cflag .eq. "  scalar") then
         xval(1) = -rx(2)
      endif
c
      return
      end

