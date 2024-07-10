      subroutine genscale(y, x, iBC)
c
c----------------------------------------------------------------------
c This subroutine calculate the y^+ and eta at inflow and internal face.
c From these generate the scaling for the inflow data.
c
c input:
c  iBC    (numnp)               : boundary condition code
c  x      (numnp,nsd)           : node coordinates
c
c output:
c  y      (numnp,ndof)          : initial values of Y variables
c
c
c Elaine Bohr december 2001
c----------------------------------------------------------------------
c
       use spebc
c       use pointer_data
       include "common.h"
       include "mpif.h"
       include "auxmpi.h"
c
       dimension y(numnp,ndof),   iBC(numnp),  
     &           x(numnp,nsd), velbarR(nfint,nflow)
       dimension yl(npro, nshl, nflow), ien(npro,nshl)
       dimension ifath(numnp),  velbarl(nelint,nshl,nflow),
     &           v1(nfint),       ymapped(numnp,ndof),
     &           shapef(nshl),	shgradl(nshl,nsd),
     &           xsi(nsd), yintl(nelint,nshl,nflow),
     &           flucl(nelint,nshl,nflow),
     &           ubarintl(nelint,nshl,nflow),
     &           fluc1(npin,nflow), fluc2(npin,nflow),
     &           ubar1(npin,nflow), ubar2(npin,nflow)
       integer   element, dir

       real*8    ymax, displTi, displTr, correction
       real*8    freestream(nflow), Qcur
c       save deltaint

	real*8 velfs   ! should be set up in solver.inp  ! 
        real*8 velmean, Tjet  ! set this in solver.inp!!! - mean velocity based on the flow rate

c       return
c        write(*,*) 'genscale l. 42'

	if (thetag.gt.360.0-0.0001) then
	  velmean = 1.100E0
	  Tjet = plandist/velmean/3.0		! time from the SPEBC inflow to the recycling plane
c Update velmean here according to the prescribed transient flow rate:
c The following parameters should be ideally read from solver.inp/input.config
	  TTbegin =  1.665    ! Time transient begins
	  TTend   =  2.000	! -------------- ends
	  FRbegin = 1.1E0       ! Flow rate mean velocity begin value
	  FRend   = 0.0E0	! ----------------------- end -------
	  if (time.ge.TTbegin.and.time.le.TTend) then   ! Flow rate transient time window
		velmean = FRbegin + 
     1 (FRend-FRbegin)*((time-TTbegin)/(TTend-TTbegin))**0.5E0
	  end if
	  if (time.gt.TTend) then
		velmean = FRend
	  end if
	end if

        ymapped(:,2:4)=y(:,1:3)
        ymapped(:,1)=y(:,4)
        ymapped(:,5)=y(:,5)
c	 write(*,*) 'deltaint = ', deltaint
c        call localy(y,      yl,     ien,    ndofl,  'gather  ')

        ubar2 = 0.0
        fluc2 = 0.0

        ymax = xyn(nfint)
c        write(*,*) 'genscale l. 51, nelint, nflow', nelint, nflow
c	write(*,*) 'ymax, nfint', ymax, nfint
c
c .... Localizing the solution vector on virtual plane
c

        do i = 1, nelint
c	write(*,*) 'i, ien2D', i, ien2D(i,:)
          do j = 1, 3
            yintl(i,:,j+1) = y(ien2D(i,:),j)
          enddo
          yintl(i,:,1) = y(ien2D(i,:),4)
          if(nflow.gt.4) then
            do j = 5, nflow
              yintl(i,:,j) = y(ien2D(i,:),j)
            enddo
          endif
c	if (i.lt.10) then
c        write(*,*) 'y = ', y(ien2D(i,2),1:4)
c        write(*,*) 'yintl[2] = ', yintl(i,:,2)
c        write(*,*) 'yintl[3] = ', yintl(i,:,3)
c        write(*,*) 'yintl[4] = ', yintl(i,:,4)
c	end if
        enddo  

c	write(*,*) 'genscale l. 67'

c
c .... Finding averaged velocity in spanwise direction
c      for the virtual plane
c

       do i=1,nfint
        velbarR(i,:)=0
        do j=1,imax(i)+1
         call shptet(ipord,xsinfin(i,j,:),shapef(:),shgradl(:,:))
         do k=1,nshl
           velbarR(i,:)=velbarR(i,:) 
     &     + yintl(elcnfin(i,j),k,:)*shapef(k)
c	  write(*,*) 'j, k, shapef, elcnfin(i,j)', j, k, shapef(k), elcnfin(i,j)
c	  write(*,*) 'i,j,k,yintl=', i,j,k,yintl(elcnfin(i,j),k,4) 
        enddo
        enddo
        velbarR(i,:)=velbarR(i,:) / (imax(i)+1)
c	write(*,*) 'i, velbarR_z, imax = ', i, velbarR(i,4), imax(i)
       enddo
 
c        write(*,*) 'genscale l. 86'
c
c .... Label the nodes that near the BL thickness
c 

       if (thetag.eq.0.0) then
         dir = 2		! dir of the flow - x coordinate
       else
         dir = 4		! This is the direction of the flow (z-coordinate)
       endif

c Assume the parabolic profile in the initial transient:
	if (time.lt.1.0*Tjet.and.thetag.gt.360.0-0.0001) then
	 do i = 1, nfint
	   velbarR(i,dir) = 2.0E0*sqrt(1.0E0-xtsC(i)/xtsC(1))
c	   write(*,*) 'time, i, v, xts, radc =', i, velbarR(i,dir) 
	 end do
	end if
c Setting the freestream velocity here (quick fix for now)
	velfs = 10.0  ! has to be more than anything for conduit flows
c	write(*,*) 'velfs = ', velfs

       v1(1)=100.0
       do i=2,nfint+1
          v1(i)=velbarR(i-1,dir)-0.99*velfs
          if((v1(i).gt.0).and.(v1(i-1).le.0)) then
             label=i-1
             go to 200
          endif
       enddo
       label=i-1

c        write(*,*) 'genscale l. 108'


c     
c.... Find the BL thickness by means of finding the y coord
c     

 200   continue
       dv=velbarR(label,dir)-velbarR(label-1,dir)
       dy=xyn(label)-xyn(label-1)

c	write(*,*) 'dv, dy = ', dv, dy
c	write(*,*) 'label, xyn(label-1) =', label, xyn(label-1)
c	write(*,*) 'vel, velbarR = ', vel, velbarR(label-1,dir) 
c     
c .... Current calculation of bl thickness at recycle plane
c

       if(istep.ne.0) then
          dlast=deltaint
          deltaint=xyn(label-1)
     &    + dy*(0.99*vel-velbarR(label-1,dir))/dv
     
c
c .... Early transients cause jumpy delta, smooth it.
c
c	  write(*,*) 'dlast, deltaint = ', dlast, deltaint
          deltaint=min(1.05*dlast,max(deltaint,0.95*dlast))
c	  write(*,*) 'deltaint = ', deltaint
       else
          deltaint=xyn(label-1)
     &             + dy*(0.99*vel-velbarR(label-1,dir))/dv
       endif
c        write(*,*) 'genscale l. 136, deltaint = ', deltaint

c
c .... Deltaint is now the ratio of BL thickness at the interior plane
c      to the BL thickness at the inlet plane
c

       deltaint=min(two*rbltin,max(deltaint,pt5*rbltin)) 
       rdelta = deltaint/rbltin
c        write(*,*) 'genscale l. 169, deltaint, rbltin = ', deltaint, rbltin
        if (thetag.gt.360.0-0.0001) rdelta = 1.0  ! no rescaling in the full cylinder case
       
c
c .... Finding freestream solutions
c

       freestream = 0
       icount = 0
       do i=1, nfint
        if (xyn(i).ge.deltaint) then
         freestream(:) = freestream(:) + velbarR(i,:)
         icount = icount + 1 
        endif
       enddo
	if (icount.gt.0)  freestream = freestream / icount
c        write(*,*) 'icount, freestream = ', icount, freestream(1)

c
c .... Putting the freestream values into the average outside the BLT
c

       do i=1, nfint
c	 write(*,*) 'i, xyn, deltaint = ', i ,xyn(i), deltaint
        if (xyn(i).ge.deltaint) then
         velbarR(i,:) = freestream(:)
        endif
       enddo

c Compute the flow rate on the recycling plane here:
        if (thetag.gt.360.0-0.0001) then   ! FULL CYLINDER CASE
	  Qcur = 0.0E0
	do i=1, nfint-1
	  Qcur = Qcur + 2.0E0*atan(1.0e0)*(velbarR(i,dir)+velbarR(i+1,dir))*
     &      abs(xtsC(i+1)**2.0-xtsC(i)**2.0)
	end do

c To get the average velocity divide by the area:
          Qcur = Qcur / (4.0E0*atan(1.0e0)*ymax**2.0E0)    
	if (time.lt.1.0*Tjet) then 
c	  if (time.lt.5.0E-01*Tjet) then
	   correction = Velmean*(time/(1.0E0*Tjet))		! Linearly Ramping up the velocity
c	  else
c	   correction = Velmean+0.5*Velmean*sin(2.0*3.141592E0*time/(0.5*Tjet)) 	! Sending a sin wave in
c	  end if
	else
	  correction = Velmean/Qcur
	end if
	 write(*,123) Tjet, Qcur, correction*Qcur, correction
  123	format('SPEBC info: Tjet=', E12.4, '; Outflow=', E12.4, 
     1   '; Inflow=', E12.4, '; Correction=', E12.4)
	end if ! Full cylinder

c
c .... Localizing the averaged velocity found above
c
c	write(*,*) 'nelint, nshl, nfint = ', nelint, nshl, nfint
       do i=1,nelint
        do k=1,nshl
         do j=1,nfint-1
          if (thetag.eq.0.0) then
           if ((x(ien2D(i,k),2).ge.xyn(j)) .and.
     &        (x(ien2D(i,k),2).le.(xyn(j+1)+0.000001))) then
            tmp = (x(ien2D(i,k),2) - xyn(j)) /
     &            (xyn(j+1) - xyn(j))
            do l=1,nflow
              velbarl(i,k,l) = 
     &              (velbarR(j+1,l) - velbarR(j,l)) * 
     &              tmp + velbarR(j,l)
            enddo
           endif
          elseif (thetag.gt.360.0-0.0001) then   ! FULL CYLINDER CASE
          if ((xcyl(ien2D(i,k),1).ge.xtsC(j+1)) .and.
     &        (xcyl(ien2D(i,k),1).le.xtsC(j))) then
            tmp = (xcyl(ien2D(i,k),1) - xtsC(j+1)) / 
     &            (xtsC(j) - xtsC(j+1))
            do l=1,nflow
              velbarl(i,k,l) = 
     &             (velbarR(j,l) - velbarR(j+1,l)) * 
     &             tmp + velbarR(j+1,l)
            enddo
           endif
	  else
           if ((xcyl(ien2D(i,k),1).ge.xcyl(nrint(j+1),1)) .and.
     &        (xcyl(ien2D(i,k),1).le.xcyl(nrint(j),1))) then
            tmp = (xcyl(ien2D(i,k),1) - xcyl(nrint(j+1),1)) /
     &            (xcyl(nrint(j),1) - xcyl(nrint(j+1),1))
            do l=1,nflow
              velbarl(i,k,l) =
     &             (velbarR(j,l) - velbarR(j+1,l)) *
     &             tmp + velbarR(j+1,l)
 	    enddo
          endif
	 end if
         enddo
        enddo
       enddo
c        write(*,*) 'genscale l. 240'
c	do i=1, 100 
c	  write(*,*) 'i,velbarl = ', i, velbarl(i,1:nshl,4)
c	end do
c
c --- For now only Blasius is coded ---
c
       
c
c .... Calculate fluctuations on elements of internal plane
c

       flucl = yintl - velbarl

c
c .... Calculate mean values on elements of internal plane
c

       ubarintl = velbarl 

c
c .... Calculating the coordinates of the point from where the
c      solution will be projected to the inlet plane
c

c	write(*,*) 'npin = ', npin

c       write(*,*) 'nrml = ', xnrml, ynrml, znrml
c       write(*,*) 'rdelta, sang = ', rdelta, sang
c	write(*,*) 'radcyl = ', radcyl

       do i=1,npin

c
c .... Cartesian coodinate system
c

        if (thetag.eq.0.0) then
          xts1 = x(nen1(i),1) + plandist
          if (xynin(i)*rdelta.gt.ymax) then
            xts2 = ymax
          else  
            xts2 = xynin(i)*rdelta
          endif
          xts3 = x(nen1(i),3)
c	write(*,*) 'i,x,y,zi, plandist =', i,xts1,xts2,xts3, plandist
c
c .... Cylindrical coordinate system
c

        else
c	write(*,*) 'i, xynin = ', i, xynin(i)
          if (xynin(i).le.0.00001) then
            xts1 = (radcyl-xynin(i)*rdelta*sang-tolerence)
     &             * cos(xcyl(nen1(i),2))
            xts2 = (radcyl-xynin(i)*rdelta*sang-tolerence)
     &             * sin(xcyl(nen1(i),2))
            xts3 = (aR-(radcyl-xynin(i)*rdelta*sang-tolerence)
     &             * (xnrml*cos(xcyl(nen1(i),2))
     &             +  ynrml*sin(xcyl(nen1(i),2))))/znrml
          elseif (xynin(i)*rdelta.gt.ymax) then
            xts1 = (radcyl-ymax*sang)
     &             * cos(xcyl(nen1(i),2))
            xts2 = (radcyl-ymax*sang)
     &             * sin(xcyl(nen1(i),2))
            xts3 = (aR-(radcyl-ymax*sang)
     &             * (xnrml*cos(xcyl(nen1(i),2))
     &             +  ynrml*sin(xcyl(nen1(i),2))))/znrml
          else
            xts1 = (radcyl-xynin(i)*rdelta*sang)
     &             * cos(xcyl(nen1(i),2))
            xts2 = (radcyl-xynin(i)*rdelta*sang)
     &             * sin(xcyl(nen1(i),2))
            xts3 = (aR-(radcyl-xynin(i)*rdelta*sang)
     &             * (xnrml*cos(xcyl(nen1(i),2))
     &             +  ynrml*sin(xcyl(nen1(i),2))))/znrml
          endif
        endif
  
c
c .... Searching for the appropriate element
c
c        write(*,*) 'genscale l. 278, i(npin), xts = '
c     &     , i, xts1, xts2, xts3

        call elem_search(xintl, xts1, xts2, xts3,
     &       xsi(:), element, 2)
        call shptet(ipord,xsi(:),shapef(:),shgradl(:,:))  

c        write(*,*) 'genscale l. 284, element = ', element

c
c .... Calculating the average velocity and fluctuations
c      for the inlet plane
c

        do k=1,nshl
          fluc2(i,:)= fluc2(i,:) + flucl(element,k,:)*shapef(k)
          ubar2(i,:)=ubar2(i,:) + ubarintl(element,k,:)*shapef(k)
c	write(*,*) 'element,k,ubarintl=', element,k,ubarintl(element,k,2:4)
        enddo
c       if (xcyl(nen1(i),1).lt.0.0002) then
c	write(*,*) 'genscale, i, xcyl, ubar2(i,1:3) = '
c     &    , i, xcyl(nen1(i),1), ubar2(i,2:4)
c	end if
       enddo  
        
c        write(*,*) 'genscale l. 297, fluc2(1,4) = ', fluc2(1, 4)

c$$$c
c$$$c keep freestream values set through averages
c$$$c
c$$$         ubaro=0
c$$$	 tbaro=0
c$$$         icount=0
c$$$         do i=1,nfin
c$$$            if(yin(i).ge.rbltin) then
c$$$               nzl=nsons(i)  !Elaine
c$$$               nzb=ienson1(i,1)
c$$$               nze=nzb+nzl-1
c$$$	       tbaro=tbaro+2.0*ubar2(i,5)+sum(fluc2(nzb:nze,5))
c$$$               ubaro=ubaro               +sum(fluc2(nzb:nze,2))
c$$$               icount=icount+nzl
c$$$            endif
c$$$         enddo
c$$$         
c$$$c     alternative to myway
c$$$c     
c$$$         ubaro=ubaro/icount
c$$$	 tmeaninflow=0.0625097048890964
c$$$	 fact= tmeaninflow/(tbaro/icount)
c$$$	 if (fact.ge. 0.9999999 .and. fact.le.1.0000001) fact = 1.0
c$$$         
c$$$         do i=1,nfin
c$$$            if(yin(i).ge.rbltin) then
c$$$               ubar2(i,2)=1.0-ubaro
c$$$            endif
c$$$         enddo
         fact = 1.0 
         rvscal = 1.0
 
c
c .... Putting the freestream value outside the BLT into ubar2
c

        do i = 1, npin
          if (xynin(i).ge.rbltin) then
            ubar2(i,:) = freestream(:)
c	    ubar2(i,dir) = 1.0
            fluc2(i,:) = 0
          endif
        enddo
c        write(*,*) 'genscale l. 342'
c$$$c
c$$$c .... For the cylindrical case the freestream velocity needs
c$$$c      to be corrected for the blockage
c$$$c
c$$$
c$$$	if (thetag.ne.0.0) then
c$$$	  displTi = 0.0
c$$$	  displTr = 0.0
c$$$	  do i = 2, nfint
c$$$
c$$$c
c$$$c .... Displacement thickness for inlet plane
c$$$c
c$$$
c$$$	    displTi = displTi + (1 - y(nrint(i),3))
c$$$     &              * (xyn(i) - xyn(i-1)) * (radcyl - xyn(i))
c$$$
c$$$c
c$$$c .... Displacement thickness for recycle plane
c$$$c
c$$$
c$$$     	    displTr = displTr + (1 - velbarR(i,4))
c$$$     &              * (xyn(i) - xyn(i-1)) * (radcyl - xyn(i))
c$$$	  enddo
c$$$c	  displTi = radcyl - sqrt(radcyl*radcyl - displTi)
c$$$c	  displTr = radcyl - sqrt(radcyl*radcyl - displTr)
c$$$	  correction = (radcyl*radcyl - displTr)
c$$$     &               / (radcyl*radcyl - displTi)
c$$$        else
	  if (thetag.lt.360.0-0.01) then		! In case of full cylinder we already defined it
            correction = 1.0
	  end if
c$$$	endif

c     
c .... Scaled plane extraction boundary condition
c     
c        write(*,*) 'genscale l. 437, fluc2(1:3,4) = ', fluc2(1:3, 4)
c        write(*,*) 'genscale l. 437, ubar2(1:3,4) = ', ubar2(1:3, 4)


         ymapped(nen1(1:npin),1)= correction * (ubar2(:,1)+fluc2(:,1))
         ymapped(nen1(1:npin),2)= correction *
     &        (ubar2(:,2)+fluc2(:,2)) 
         ymapped(nen1(1:npin),3)= correction *
     &        (ubar2(:,3)+fluc2(:,3))*rvscal
	  if (thetag.gt.360.0-0.001.and.time.lt.1.0*Tjet) then
           ymapped(nen1(1:npin),4)= correction*ubar2(:,4)
	  else
           ymapped(nen1(1:npin),4)= 
     1 abs(correction * (ubar2(:,4)+fluc2(:,4)))  ! ensure positive inflow
	 end if
         ymapped(nen1(1:npin),5)= correction * fact*(ubar2(:,5) +
     &                            fluc2(:,5))

c	write(*,*) 'ymapped = ', ymapped(nen1(1:3),4)
c     
c .... Ready to put the solution on the inflow plane 
c
c	write(*,*) 'intpres = ', intpres
      if(intpres.eq.1) then     !interpolating pressure at inflow
         where (btest(iBC,11))
            y(:,1) = ymapped(:,2)
            y(:,2) = ymapped(:,3)
            y(:,3) = ymapped(:,4)
            y(:,4) = ymapped(:,1)
            y(:,5) = ymapped(:,5)
         endwhere
      else                      ! not interpolating pressure at inflow
c        write(*,*) 'genscale, 410; ',
c     1  'solution is assighed at inflow, y(2) = ', ymapped(1:3,2)
c        write(*,*) 'genscale, 410; ',
c     1  'solution is assighed at inflow, y(3) = ', ymapped(1:3,3)
c        write(*,*) 'genscale, 410; ',
c     1  'solution is assighed at inflow, y(4) = ', ymapped(1:3,4)
         where (btest(iBC,11))
            y(:,1) = ymapped(:,2)
            y(:,2) = ymapped(:,3)
            y(:,3) = ymapped(:,4)
            y(:,5) = ymapped(:,5)
         endwhere
      endif

c     
c     debugging variables
c     
c      if(iter.eq.nitr) then
c         write(555,556)lstep+1,deltaint,label,nfint
c         write(554,556)lstep+1,yplusi(2),ypluso(2),factu,factt,gamt
c         call flush(554)
c         call flush(555)
c      endif
c 556  format(i6,5(2x,e14.7))
c
c.... return
c

      return
      end
