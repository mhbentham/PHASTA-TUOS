      subroutine renum_cyl(x)
c
c----------------------------------------------------------------------
c This subroutine finds all nodes that are on the inlet plane and all
c nodes that are on the recycle plane; it also blocks elements that 
c are on recycle plane; all nodes are also stored with cylindrical 
c coordinates; find all father nodes for recycle plane (i.e. for
c theta = -theta given)
c
c input:
c  x      (numnp,nsd)           : node coordinates
c
c output:
c  xcyl   (numnp,nsd)           : node cylindrical coordinates 
c  ien2D  (npro, nshl)		: connectivity array for recycle plane
c				  assuming tethraheadral elements, i.e.
c				  triangular elements on face
c
c
c  Elaine Bohr
c  June 2002
c
c  Igor Bolotnov
c  Summer 2009
c----------------------------------------------------------------------
c
       use spebc
c       use pointer_data
       include "common.h"
       include "mpif.h"
       include "auxmpi.h"
c
        dimension x(numnp,nsd), nrin(numnp), nula(numnp),
     &		  erreur(nshl), xtmp(nsd), nrin2(numnp)

        integer  temp, etmp, s
	real*8 thint

	real*8,  allocatable :: xrtmp(:)
	real*8,  allocatable :: Radd(:)        
	real*8,  allocatable :: xtsXn(:), xtsYn(:), xtsCn(:)

c	thetag = thetag/180.0*pi


c
c .... changing to cylindrical coordinate system for nodal point
c       
	
	xcyl(:,1) = sqrt(x(:,1)*x(:,1) + x(:,2)*x(:,2))
	
	j = 0
        do i=1,numnp
	  if ((x(i,1).eq.0).and.(x(i,2).eq.0)) then
	    j = j+1
	    nula(j) = i
	  else
	    xcyl(i,2) = atan2(x(i,2),x(i,1)) 
	  endif
	enddo
        xcyl(:,3) = x(:,3)
	

c
c .... finding the minimum and maximum angles
c

	thmin = xcyl(nen1(1),2)
        thmax = xcyl(nen1(1),2)
c	thabsmin = xcyl(nen1(1),2)
        do i=2,npin
	  if((x(nen1(i),1).eq.0).and.(x(nen1(i),2).eq.0)) then
	    goto 10
	  else
            if(thmin.gt.xcyl(nen1(i),2)) thmin = xcyl(nen1(i),2)
            if(thmax.lt.xcyl(nen1(i),2)) thmax = xcyl(nen1(i),2)
c	    if(abs(thabsmin).gt.abs(xcyl(nen1(i),2))) thabsmin = xcyl(nen1(i),2) 
	  endif
 10	  continue
        enddo
	
	do i=1,j
	  xcyl(nula(i),2) = thmax
	enddo
	
c
c .... Finding nodes from inlet plane with same theta given
c
        if (abs(thmin).gt.abs(thmax)) then
		thint = thmin
	  else
		thint = thmax
	end if 

c	write(*,*) 'Printing five variables here: '
c        write(*,*) 'thmin, thmax, npin, thint = ', thmin, thmax, npin, thint  !, thabsmin

       j = 0
       do i=1,npin
         error1 = abs(xcyl(nen1(i),2) - thint)
c	 write(*,*) 'i, error1, xcyl(2) = ', i, error1, xcyl(nen1(i),2)
	 if(error1.le.0.001) then
	   j = j + 1
	   nrin(j)=nen1(i)
	 endif
       enddo
       nfin = j
c	write(*,*) 'nfin = ', nfin
	if (nfin.eq.1) 
     &      write(*,*) 'Warning: cannot process the BL mesh expansion'
c Igor Bolotnov, summer 2009.
c Consider the special case of thetag = 360.
c The radii outside of the BL has to be generated since they cannot be 
c found using the conventional procedure (there is no edge in the full cylinder case).

c We will use the spacing in the outside edge of the BL to fill in the radial coordinates near the center.

	if (thetag.gt.360-0.001) then
c Step 1: find the smallest radius (the farest from the wall)
	 nlast = nrin(1)
	  do i=2,nfin
	     if (xcyl(nrin(i),1).lt.xcyl(nlast,1)) then
		nlast = nrin(i)
	     end if  
	  end do
c Step 2: find next one to the smallest one:
         if (nlast.ne.nrin(1)) then
		nnext = nrin(1)
	    else
		nnext = nrin(2)
	 end if
          do i=1,nfin
        if (xcyl(nrin(i),1).lt.xcyl(nnext,1)
     1   .and.nrin(i).ne.nlast) then
                nnext = nrin(i)
             end if
          end do

	 write(*,*) 'rlast, rnext = ', xcyl(nlast, 1), xcyl(nnext, 1)

c Step 3: Compute the number of additional radial locations:
	 Naddrad = int(xcyl(nlast,1)/abs(xcyl(nlast,1)-xcyl(nnext,1)))+1
c         write(*,*) 'abs = ', abs(xcyl(nlast,1)-xcyl(nnext,1))
c	 write(*,*) '/ = ', xcyl(nlast,1)/abs(xcyl(nlast,1)-xcyl(nnext,1))
c	 write(*,*) 'int = ', int(xcyl(nlast,1)/abs(xcyl(nlast,1)-xcyl(nnext,1)))
	 write(*,*) 'Adding ', Naddrad, ' points in the unstructured',
     1      ' mesh region'

	 allocate(Radd(Naddrad))

c Step 4: Generate the additional locations:
c	 write(*,*) 'xcyl_last, next = ', xcyl(nlast,1), xcyl(nnext,1)
	 do i=1, Naddrad
	  Radd(Naddrad-i+1) = real(i-1)/real(Naddrad)*xcyl(nlast,1)
	 write(*,*) 'i, r = ', i, Radd(Naddrad-i+1) 
         end do

c Step 5: find a closest node on the inlet plane (npin) for each Radd and record its number
c	 do i=1, Naddrad
c	  nrin2(i) = nen1(1)
c	  do j = 2, npin
c	   if (abs(xcyl(nen1(j),1)-Radd(i)).lt.abs(xcyl(nrin2(i),1)-Radd(i))) then
c		nrin2(i) = nen1(j)
c 	   end if
c	  end do
c	  write(*,*) 'i, nrin2, r, th = ', i, nrin2(i), xcyl(nrin2(i),1:2)
c	 end do

c Step 6: Increase nfin and expand nrin:
c	 nfin = nfin + Naddrad
c	 nrin(nfin-Naddrad+1:nfin) = nrin2(1:Naddrad)

c Step 5a. Generate x,y arrays for the radial points to be used instead of the nrint array.

	 allocate(xtsX(nfin+Naddrad))
         allocate(xtsY(nfin+Naddrad))
	 allocate(xtsC(nfin+Naddrad))

c   Add the old points in:
	 do i = 1, nfin
	  xtsX(i) = x(nrin(i),1)
	  xtsY(i) = x(nrin(i),2)
	  xtsC(i) = xcyl(nrin(i),1)
	 end do

c The new points go in right after 
	 do i = 1, Naddrad
	  xtsX(nfin+i) = Radd(i)*cos(thmin)
	  xtsY(nfin+i) = Radd(i)*sin(thmin)
	  xtsC(nfin+i) = Radd(i)
	 end do

	 nfin = nfin + Naddrad

	 deallocate(Radd)

	end if

       allocate(xrtmp(nfin))
        if (thetag.gt.360-0.001) then
        do i=1,nfin
         xrtmp(i) = xtsC(i)
c         write(*,*) 'i, xrtmp = ', i, xrtmp(i)
       enddo

	 else
       do i=1,nfin
	 xrtmp(i) = xcyl(nrin(i),1)
c	 write(*,*) 'i, xrtmp = ', i, xrtmp(i)
       enddo
	end if	
       
c
c .... Ordering nrin by decreasing radius
c

	  if (thetag.le.360-0.001) then
	j = 0
	allocate (nrint(nfin))
	
	radinl = xcyl(nrin(1),1)
	do i=1,nfin
	  if (radinl.lt.xcyl(nrin(i),1)) radinl = xcyl(nrin(i),1)
	  rmaxtemp=xrtmp(1)
	  itmp = 1
	  do k=2,nfin
	    if (rmaxtemp.le.xrtmp(k)) then
	      rmaxtemp = xrtmp(k)
	      itmp = k
	    endif
	  enddo
	  if (radcyl.ge.xrtmp(itmp)) then
	    j = j + 1
	    nrint(j)=nrin(itmp)
	  endif
	  xrtmp(itmp) = -1
	enddo
	nfint = j
	 else   ! FULL CYL here
        j = 0
        allocate (nrint(nfin))
         allocate(xtsXn(nfin))
         allocate(xtsYn(nfin))
         allocate(xtsCn(nfin))

        
        radinl = xtsC(1)
        do i=1,nfin
          if (radinl.lt.xtsC(i)) radinl = xtsC(i)
          rmaxtemp=xrtmp(1)
          itmp = 1
          do k=2,nfin
            if (rmaxtemp.le.xrtmp(k)) then
              rmaxtemp = xrtmp(k)
              itmp = k
            endif
          enddo
          if (radcyl.ge.xrtmp(itmp)) then
            j = j + 1
            nrint(j) = nrin(itmp)
	    xtsXn(j) = xtsX(itmp)
            xtsYn(j) = xtsY(itmp)
            xtsCn(j) = xtsC(itmp)
          endif
          xrtmp(itmp) = -1
        enddo

	xtsX = xtsXn
	xtsY = xtsYn
	xtsC = xtsCn

	deallocate(xtsXn)
        deallocate(xtsYn)
        deallocate(xtsCn)


	 end if ! not a full cylinder/FULL CYL
	deallocate(xrtmp)

c
c .... off wall coordinate for the inlet plane
c
c	write(*,*) 'radinl = ', radinl
	do i=1,npin
	  xynin(i) = (radinl - xcyl(nen1(i),1))/sang
c	 write(*,*) 'i, xynin = ', i, xynin(i)
     	enddo

c
c .... off wall coordinate for the virtual points on recycle plane
c
        if (thetag.gt.360-0.001) then
        do i=1,nfin
          xyn(i) = (radcyl - xtsC(i))/sang
	 write(*,*) 'i, xtsC, xyn = ', i, xtsC(i), xyn(i)
        enddo
	 else
	do i=1,nfint
	  xyn(i) = (radcyl - xcyl(nrint(i),1))/sang
     	enddo
	end if
c
c .... Finding corresponding elements and local coordinates
c      on recycle plane for every arc with radius from nrint
c
        if (thetag.gt.360-0.001) nfint = nfin

       s = (thmax-thmin)*radcyl/ds
       allocate (xsinfin(nfint,s+1,nsd))
       allocate (elcnfin(nfint,s+1))
       allocate (imax(nfint))

        if (thetag.gt.360-0.001) then   ! THIS IS THE FULL CYLINDER CASE:
       do jj = 1, nfint
          xts1 = xtsX(jj)
          xts2 = xtsY(jj)
          xts3 = x(nrin(1),3) + plandist
          call elem_search(xintl, xts1, xts2, xts3,
     &                     xtmp(:), etmp, 1)
          xsinfin(jj,1,:) = xtmp(:)
          elcnfin(jj,1) = etmp
          imax(jj) = (thmax-thmin)*xtsC(jj)/ds
          do i=1,imax(jj)
            if ( xtsC(jj) .eq. radcyl) then
              xts1 = (xtsC(jj)-tolerence)
     &          *cos(1.0*i/imax(jj)*(thmax-thmin)+thmin)
              xts2 = (xtsC(jj)-tolerence)
     &          *sin(1.0*i/imax(jj)*(thmax-thmin)+thmin)
            else
c             reel=i*ds/xcyl(nrint(jj),1)
              xts1 = xtsC(jj)
     &          *cos(1.0*i/imax(jj)*(thmax-thmin)+thmin)
              xts2 = xtsC(jj)
     &          *sin(1.0*i/imax(jj)*(thmax-thmin)+thmin)
            endif
            xts3 = x(nrin(1),3) + plandist
            call elem_search(xintl, xts1, xts2, xts3,
     &                       xtmp(:), etmp, 1)
            xsinfin(jj,1+i,:) = xtmp(:)
            elcnfin(jj,1+i) = etmp
          enddo
        enddo
	else
       do jj = 1, nfint
          xts1 = x(nrint(jj),1) !*cos(xcyl(nrint(jj),2)+tolerence)
          xts2 = x(nrint(jj),2) !*sin(xcyl(nrint(jj),2)+tolerence)
          xts3 = x(nrint(jj),3) + plandist
	  call elem_search(xintl, xts1, xts2, xts3,
     &		           xtmp(:), etmp, 1)
     	  xsinfin(jj,1,:) = xtmp(:)
	  elcnfin(jj,1) = etmp
	  imax(jj) = (thmax-thmin)*xcyl(nrint(jj),1)/ds
     	  do i=1,imax(jj)
	    if ( xcyl(nrint(jj),1) .eq. radcyl) then
	      xts1 = (xcyl(nrint(jj),1)-tolerence) 
     &		*cos(1.0*i/imax(jj)*(thmax-thmin)+thmin)
	      xts2 = (xcyl(nrint(jj),1)-tolerence) 
     &		*sin(1.0*i/imax(jj)*(thmax-thmin)+thmin)
	    else
c	      reel=i*ds/xcyl(nrint(jj),1)
	      xts1 = xcyl(nrint(jj),1)
     &		*cos(1.0*i/imax(jj)*(thmax-thmin)+thmin)
	      xts2 = xcyl(nrint(jj),1)
     & 		*sin(1.0*i/imax(jj)*(thmax-thmin)+thmin)
	    endif
	    xts3 = x(nrint(jj),3) + plandist
	    call elem_search(xintl, xts1, xts2, xts3,
     &		             xtmp(:), etmp, 1)
     	    xsinfin(jj,1+i,:) = xtmp(:)
	    elcnfin(jj,1+i) = etmp
          enddo
	enddo 

	end if

        return
        end


      subroutine renum_cart(x)
c
c----------------------------------------------------------------------
c This subroutine finds all nodes that are on the inlet plane and all
c nodes that are on the recycle plane; it also blocks elements that 
c are on recycle plane; all nodes are also stored with cylindrical 
c coordinates; find all father nodes for recycle plane (i.e. for
c theta = -theta given)
c
c input:
c  x      (numnp,nsd)           : node coordinates
c
c output:
c  xcyl   (numnp,nsd)           : node cylindrical coordinates 
c  ien2D  (npro, nshl)		: connectivity array for recycle plane
c				  assuming tethraheadral elements, i.e.
c				  triangular elements on face
c
c
c  Elaine Bohr
c  July 2002
c----------------------------------------------------------------------
c
       use spebc
c       use pointer_data
       include "common.h"
       include "mpif.h"
       include "auxmpi.h"
c
        dimension x(numnp,nsd), nyin(numnp),
     &		  erreur(nshl), xtmp(nsd), yrtmp(numnp)

        integer  temp, etmp, s

c
c .... Finding nodes from inlet plane with minimal z value
c
       j = 0
       zmin = x(nen1(1),3)
       zmax = x(nen1(1),3)
       do i=2,npin
	 if(x(nen1(i),3).lt.zmin) zmin = x(nen1(i),3)
	 if(x(nen1(i),3).gt.zmax) zmax = x(nen1(i),3)
       enddo

       do i=1,npin
         if (x(nen1(i),3).eq.zmin) then
	   j = j + 1
	   nyin(j) = nen1(i)
	 endif  
       enddo
       nfin = j
       do i=1,nfin
	 yrtmp(i) = x(nyin(i),2)
       enddo	
       
c
c .... Ordering nyin by increasing y
c
	j = 0
	allocate (nrint(nfin))
	
	do i=1,nfin
	  rmintemp=yrtmp(1)
	  itmp = 1
	  do k=2,nfin
	    if (rmintemp.ge.yrtmp(k)) then
	      rmintemp = yrtmp(k)
	      itmp = k
	    endif
	  enddo
	  j = j + 1
	  nrint(j)=nyin(itmp)
	  yrtmp(itmp) = 10000000
	enddo
	nfint = j

c
c .... y coordinate for the inlet plane
c
	do i=1,npin
	  xynin(i) = x(nen1(i),2)
     	enddo

c
c .... y coordinate for the virtual points on recycle plane
c
	
	do i=1,nfint
	  xyn(i) = x(nrint(i),2)
     	enddo

c
c .... Finding corresponding elements and local coordinates
c      on recycle plane for every arc with radius from nrint
c
       s = (zmax - zmin) / ds
       allocate (xsinfin(nfint,s+1,nsd))
       allocate (elcnfin(nfint,s+1))
       allocate (imax(nfint))

       do jj = 1, nfint
            
          xts1 = x(nrint(jj),1) + plandist
          xts2 = x(nrint(jj),2) 
          xts3 = x(nrint(jj),3) 
	  call elem_search(xintl, xts1, xts2, xts3,
     &		           xtmp(:), etmp, 1)
     	  xsinfin(jj,1,:) = xtmp(:)
	  elcnfin(jj,1) = etmp
	  imax(jj) = s
     	  do i=1,imax(jj)
	    xts1 = x(nrint(jj),1) + plandist
	    xts2 = x(nrint(jj),2) 
	    xts3 = x(nrint(jj),3) + i*ds
	    call elem_search(xintl, xts1, xts2, xts3,
     &		             xtmp(:), etmp, 1)
     	    xsinfin(jj,1+i,:) = xtmp(:)
	    elcnfin(jj,1+i) = etmp
          enddo
	enddo 

        return
        end
