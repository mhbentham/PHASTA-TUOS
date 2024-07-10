        subroutine genini (iBC, BC, y, ac, banma, iper, ilwork,
     &               ifath, velbar,  nsons,x,
     &               shp,     shgl,    shpb,    shglb ) 
c
c----------------------------------------------------------------------
c This routine reads the initial values in primitive form (density,
c velocity and temperature), satisfies the boundary conditions and 
c converts them to Y-variables.
c
c input:
c  iBC    (nshg)               : boundary condition code
c  BC     (nshg,ndofBC)        : boundary condition constrain data
c   x     (numnp,nsd)	       : locations of nodes, numnp-> # of node
c				  nsd-> space dimension, 1=x, 2=y, 3=z 
C			
c
c output:
c  y      (nshg,ndof)          : initial values of Y variables
c
c
c Farzin Shakib, Winter 1986.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use specialBC   ! gets itvn from here
        use convolImpFlow ! brings in ntimeptpT and other variables
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        dimension iBC(nshg),                iper(nshg),
     &            BC(nshg,ndofBC),          y(nshg,ndof),
     &            ac(nshg,ndof),            x(numnp,nsd),
     &            shp(MAXTOP,maxsh,MAXQPT),
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)

c
        dimension banma(nshg,1)
c
        dimension ilwork(nlwork)
c
c  stuff for dynamic model s.w.avg and wall model
c
        dimension ifath(numnp),    velbar(nfath,nflow),
     &           nsons(nfath) 

        character*20 fname1
        character*10 cname2
        character*5 cname
c
c.... -------------------------->  Restart  <---------------------------
c
c.... read q from [RESTAR.INP], reset LSTEP
c
        call restar ('in  ',  y,  ac)

c.... copy marker field from array y and restore sclr 2, Jun, Jul 2014
       if (ndof.eq. 7) then
          banma(:,1) = y(:,7)
          y(:,7) = y(:,6)
       endif

c
c	if (myrank.eq.master) write(*,*) 'y1 = ', y(1:10,1)

        if((itwmod.gt.0) 
     &     .or. (nsonmax.eq.1 .and. (iLES.gt.0.or.iDNS.ne.0)) ) then 
           call rwvelb('in  ',velbar,ifail)
c
c if the read failed calculate velbar
c
           if(ifail.eq.1) then
              call getvel (y,     ilwork, iBC,
     &                     nsons, ifath, velbar)
           endif
 
        endif   ! for the itwmod or irscale

c        if (myrank.eq.master) write(*,*) 'y2 = ', y(1:10,1)
c        if (myrank.eq.master) write(*,*) 'velbar = ', velbar(1:10,2)


c
c.... time varying boundary conditions as set from file bct.dat and impt.dat 
c     (see function for format in file bctint.f)
c
        if (itvn .gt. 0 ) then !for inlet velocities
           call initBCt( x, iBC, BC)
           call BCint(time,shp,shgl,shpb,shglb,x, BC, iBC)
        endif
        if (impfile .gt. 0 ) then !for impedance BC
           if(myrank.eq.master) then 
              write(*,*) 'reading Qhistor.dat'
           endif
           open(unit=816, file='Qhistor.dat',status='old')
           read (816,*) ntimeptpT
           allocate (QHistImp(ntimeptpT+1,numImpSrfs)) 
           do j=1,ntimeptpT+1
              read(816,*) (QHistImp(j,n),n=1,numImpSrfs) !read flow history
           enddo
           close(816)           !vie added
                 
           call initImpt()      !read impedance data and initialize begin/end values
           do i=2,ntimeptpT
              call Impint((ntimeptpT-i+1)*Delt(1),i) !return Imp values in reverse order ZN->Z0
           enddo

           allocate (pold(0:MAXSURF)) !for pressure part that depends on the history only
        endif
        if (itvbc .gt. 0 ) then
           call initsetBC(iBC)
           call setBC( time, x, iBC, BC)
        endif
c
c
c.... satisfy the boundary conditions
c

        call itrBC (y, ac, banma,  iBC, BC, iper, ilwork)
c	if (myrank.eq.master) write(*,*) 'itrBC is done.'

        itempr=mod(impl(1),2)  ! tempr solve if impl odd
        if(itempr.eq.1) then
           isclr=0
           call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
        endif
        do isclr=1,nsclr
           call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
        enddo

        if((irscale.ge.0) .and. (myrank.eq.master)) then
         write(*,*) 'genscale...'
           call genscale(y, x, iBC)
         write(*,*) 'is done.'
        endif
c
c.... --------------------------->  Echo  <----------------------------
c
c.... echo the initial data
c
        if ((necho .lt. 0).and.(myrank.eq.master)) then
          do n = 1, nshg
            if (mod(n,50) .eq. 1) write(iecho,1000) ititle,(i,i=1,ndof)
            write (iecho,1100) n, (y(n,i),i=1,ndof)
          enddo
        endif
c
c.... return
c
        return
c
1000    format(a80,//,
     &  ' I n i t i a l   V a l u e s                        ',//,
     &  '    Node ',/,
     &  '   Number ',6x,6('dof',i1,:,10x))
1100    format(1p,2x,i5,5x,5(e12.5,2x))
c
        end
