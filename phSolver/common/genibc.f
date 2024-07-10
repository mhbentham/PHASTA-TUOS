        subroutine geniBC (iBC,x)
c
c----------------------------------------------------------------------
c This routine reads the boundary condition codes.
c
c output: 
c  iBC   (nshg)        : Boundary Condition code
c
c         = 1 * iBC_1 + 2 * iBC_2 + 4 * iBC_3
c              density   temperature   pressure
c
c    if nsd = 3:
c
c        +  8 * iBC_4 +  16 * iBC_5 +  32 * iBC_6
c           x1-velocity   x2-velocity   x3-velocity
c
c        + 64 * iBC_7 + 128 * iBC_8 + 256 * iBC_9 + 512 * iBC_10
c          sclr1         sclr2        sclr3         sclr4
c
c        + 1024 * iBC_11  + 2048* iBC_12 
c          perioidicity     spebc          
c
c  nBC   (nshg)        : Boundary Condition mapping array
c
c
c Farzin Shakib, Winter 1986.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
c
        use readarrays          ! used to access iBCtmp
        include "common.h"
c
c Arrays in the following 1 line are now dimensioned in readnblk
c        dimension iBCtmp(numpbc)
c
        dimension iBC(nshg)
        dimension itemp(6)

        integer, parameter :: maxp = 1000                              !Farhad
        dimension x(numnp,nsd)                                         !Farhad
        integer n, m, m1, m2, CrkEdgNodes                              !Farhad
        real*8  CrkTime, RadiusF                                       !Farhad
        real*8  Vec1(2), Vec2(2), Angle(nshg), xEdge(maxp,3)           !Farhad
        real*8  Magnitude1, Magnitude2, CrossProd, DotProd, Ang        !Farhad
	   real*8 xleft, xright, zleft, zright		! Igor
	   real*8 rx(maxp), ry(maxp), rz(maxp), Msize 
	   integer Mtype, Mnum
        logical exts                                                   !Farhad    

c.... set the iBC array
        iBC = 0

c       if(numpbc.eq.0) return  ! sometimes there are no BC's on a partition !Farhad. I commented this line out
        where (nBC(:) .ne. 0) iBC(:) = iBCtmp(nBC(:))

        inquire(file='FSHAP_R.dat', exist=exts)                         ! Igor: New format for rough channel flows
        if ((exts).and.(myrank.eq.master)) then                         !Farhad
           write(*,*) 'Roughness geometry file is available'            !Farhad
        endif                                                           !Farhad
! Igor:
	if (exts) then 
        open(unit=89, file ='FSHAP_R.dat', status='old')                  !Igor
	   do m = 1,3
		read(89, *)
	   end do
	   read(89,*) Mtype
	 if ((exts).and.(myrank.eq.master)) write(*,*) 'Mtype = ', Mtype
         if (Mtype.eq.1) then
         read(89,*) Msize, Mnum	! Type: 1 - rectangular, 2 - hemispheres; Size is height/radius, Mnum is number of items
	   do m = 1, Mnum
		read(89,*) rx(m), ry(m), rz(m)
	   end do
!	 write(*,*) 'Mnum = ', Mnum
	 else if (Mtype.eq.2) then
	   read(89,*) Msize, Mx, My, Mz
	   do m = 1, Mx
		read(89,*) rx(m)
	   end do
	   do m = 1, My
		read(89,*) ry(m)
	   end do
	   do m = 1, Mz
		read(89,*) rz(m)
	   end do
	 end if      
	end if ! exts
      if (exts) then  ! Igor
       do n = 1, nshg
       RadiusF = sqrt((x(n,3)-pi/2)**2.0 + 
     1  (x(n,2))**2.0 + (x(n,1)-pi)**2.0)    !  Sphere shaped - channel flow turbulence
         if (RadiusF .lt. 0.25) then    ! Igor, sphere shaped for channel flow turbulence
            iBC(n)=56
         end if
       end do
       if (Mtype.eq.2) then   ! hemispherical roughness
	if (myrank.eq.master) write(*,*) 'Hemis: ', Mx, My, Mz
        do n = 1, nshg                                                    
		do m = 1, Mx
		do m1 = 1, My
		do m2 = 1, Mz
            RadiusF = sqrt((x(n,1)-rx(m))**2.0 
     1 + (x(n,2)-ry(m1))**2.0 + (x(n,3)-rz(m2))**2.0)    
           if (RadiusF .le. Msize) then   
         write(*,*) 'Fixed point x, y, z: ', x(n,1), x(n,2), x(n,3)
                 iBC(n)=56                                              
          endif
        enddo  ! m2
        enddo  ! m1
        enddo  ! m
       enddo   ! n
	 elseif (Mtype.eq.1) then   ! Rectangular roughness
        do n = 1, nshg                                                    
		do m = 1, Mnum
		xleft = x(n,1) - (rx(m)-0.5E0*Msize)
		xright = (rx(m)+0.5E0*Msize) - x(n,1)
		yleft = x(n,2) - (ry(m)-Msize)
		yright = (ry(m)+Msize) - x(n,2)
!	   if (m.lt.10.and.n.lt.10) write(*,*) 'n,m,xleft,xright: ', n,m,xleft,xright,yleft,yright
		if (xleft.gt.0.and.xright.gt.0.and.yleft.gt.0.and.yright.gt.0) then   
!		  write(*,*) 'Fixed point N, x, y, z: ', n, x(n,1), x(n,2), x(n,3)
		   iBC(n)=56                                              
          endif
        enddo  ! m
       enddo   ! n
	 end if ! Mtype
	end if ! Igor, exts

c
c.... echo the input iBC array only if other than zero
c
        if (necho .lt. 3) then
          nn = 0
          do n = 1, nshg
            if (nBC(n) .ne. 0) then
              nb = nBC(n)
              nn = nn + 1
              if (mod(nn,50).eq.1) write(iecho,1000)ititle,(j,j=1,ndof)
              itemp(   1) = mod(iBCtmp(nb)   ,2) - mod(iBCtmp(nb)/ 4,2)
              itemp(   2) = mod(iBCtmp(nb)/ 8,2)
              itemp(   3) = mod(iBCtmp(nb)/16,2)
              itemp(   4) = mod(iBCtmp(nb)/32,2)
              itemp(ndof) = mod(iBCtmp(nb)/ 2,2)
              write(iecho,1100) n,(itemp(i),i=1,ndof)
            endif
          enddo
        endif
        deallocate(iBCtmp)
c
c.... return
c
        return
c
c.... end of file error handling
c
999     call error ('geniBC  ','end file',ibndc)
c
1000    format(a80,//,
     &  ' N o d a l   B o u n d a r y   C o n d i t i o n   C o d e',//,
     &  '    Node   ',13x,6('dof',i1,:,6x))
1100    format(2x,i5,10x,5i10)
c
        end
