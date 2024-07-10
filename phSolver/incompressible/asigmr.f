        subroutine AsIGMR (y,         ac,         banma,
     &                     x,         xmudmi,
     &                     shp,       shgl,       ien,     
     &                     res,       qres,
     &                     xKebe,     xGoC,       rerr,
     &                     cfl,       icflhits,   elemvol_local,
     &                     xarray,    yarray,     zarray,  
     &                     bubradius, bubradius2, coordtag)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use stats
      use rlssave  ! Use the resolved Leonard stresses at the nodes.
      use timedata    ! time series
      use turbsa      ! access to d2wall and effvisc & x2wall, y2wall, z2wall


      include "common.h"
c
        dimension y(nshg,ndofl),              ac(nshg,ndofl),
     &            x(numnp,nsd),              
     &            shp(nshl,ngauss),           shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),
     &            res(nshg,nflow),
     &            qres(nshg,idflx),           cfl(nshg),
     &            icflhits(nshg)
c
        dimension banma(nshg,1),              bml(npro,nshl,1)
c
        dimension yl(npro,nshl,ndofl),         acl(npro,nshl,ndofl),
     &            xl(npro,nenl,nsd),           dwl(npro,nenl),      
     &            rl(npro,nshl,nflow), 
     &            ql(npro,nshl,idflx),
     &            evl(npro,nenl),
     &            cfll(npro,nshl),
     &            xwl(npro,nenl),ywl(npro,nenl),zwl(npro,nenl)
c        
        dimension xKebe(npro,9,nshl,nshl), 
     &            xGoC(npro,4,nshl,nshl)
c
        dimension rlsl(npro,nshl,6) 

c
        real*8    lStsVec(npro,nshl,nResDims)
        
        dimension xmudmi(npro,ngauss)
        dimension sgn(npro,nshl)
c
        real*8 rerrl(npro,nshl,6), rerr(nshg,numerr)

c!... Matt Talley's Bubble Coalescence Control
        real*8 xarray(ibksiz), yarray(ibksiz), zarray(ibksiz)
        real*8 bubradius, bubradius2

        integer coordtag(ibksiz)
c
c!.... gather the variables
c
c
c!.... get the matrix of mode signs for the hierarchic basis functions. 
c
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
        
        call localy(y,      yl,     ien,    ndofl,  'gather  ')
        call localy(ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
        call local (qres,   ql,     ien,    idflx,  'gather  ')
        if (iRANS .eq. -2.or.((iDNS.gt.0).and.(abs(itwmod).eq.1))
     &      .or. iBT.eq.1) then ! kay-epsilon & DNS slip-vel
           call localx (d2wall,   dwl,     ien,    1,     'gather  ')
        endif
 
        if( (iLES.gt.10).and.(iLES.lt.20)) then  ! bardina 
           call local (rls, rlsl,     ien,       6, 'gather  ')  
        else
           rlsl = zero
        endif      

        if ((iDNS.gt.0).and.(itwmod.eq.-2)) then
          call local(effvisc, evl,    ien,    1,      'gather  ')
        endif

c	write(*,*) 'iLES = ', iLES

        if ((iDNS.gt.0).and.(abs(itwmod).eq.1)) then
         call local(x2wall, xwl,    ien,    1,      'gather  ')
         call local(y2wall, ywl,    ien,    1,      'gather  ')
         call local(z2wall, zwl,    ien,    1,      'gather  ')
        endif

c
c.... zero the matrices if they are being recalculated
c
        if (lhs. eq. 1)  then
           xKebe = zero
           xGoC  = zero
        endif   
c
!       Update the marker field.
       IF(iBT .eq. 1) THEN
          call banmaUpdate(xl, yl, banma, ien, bml)
       ENDIF !iBT

c.... get the element residuals, LHS matrix, and preconditioner
c
        rl     = zero
        cfll   = zero

        if(ierrcalc.eq.1) rerrl = zero

        call e3  (yl,      acl,     dwl,     shp,
     &            shgl,    xl,      rl,      
     &            ql,      xKebe,   xGoC,    xmudmi, 
     &            sgn,     rerrl,   rlsl,
     &            cfll,    evl,     xwl,     ywl,
     &            zwl,     xarray,  yarray,  zarray,
     &            elemvol_local,    bml,
     &            bubradius, bubradius2, coordtag)
c
c.... assemble the statistics residual
c
        if ( stsResFlg .eq. 1 ) then
           call e3StsRes ( xl, rl, lStsVec )
           call local( stsVec, lStsVec, ien, nResDims, 'scatter ')
        else
c
c.... assemble the residual
c
           call local (res,    rl,     ien,    nflow,  'scatter ')
           
           if ( ierrcalc .eq. 1 ) then
              call local (rerr, rerrl,  ien, 6, 'scatter ')
           endif
        endif
c
c.... sum the CFL value from IPs.  These wil be divided by the number of
c     contributors in elmgmr to get average CFL value at node
c
        call localSum (cfl, cfll, ien, icflhits, 1)
c
c.... end
c
        if (exts) then
           if ((iter.eq.1).and.(mod(lstep,freq).eq.0)) then
c	  if (myrank.eq.master) write(*,*) 'Calling the timeseries now:'
              call timeseries(yl,xl,ien,sgn,1,4)
              call timeseries(ql(:,:,10),xl,ien,sgn,5,5)
              call timeseries(ql(:,:,11:19),xl,ien,sgn,6,14)
c	Two-phase flows averaging:
	      if (numvar.ge.15) then
		call timeseries(yl(:,:,6),xl,ien,sgn,15,15)
	      endif
              if (numvar.eq.19) then
                call timeseries(ql(:,:,idflx-2:idflx),xl,ien,sgn,16,18)
              endif
           endif
        endif
        
        return
        end



C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c-----------------------------------------------------------------------
c=======================================================================


        subroutine AsIGMRSclr(y,       ac,      x,       
     &                     shp,     shgl,    ien,     
     &                     res,     qres,    xSebe, xmudmi,
     &                     cfl,     icflhits,  cflold)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use     turbSA   ! access to d2wall and effvisc
      include "common.h"
c
        dimension y(nshg,ndofl),              ac(nshg,ndofl),
     &            x(numnp,nsd),              
     &            shp(nshl,ngauss),           shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),
     &            res(nshg),                  qres(nshg,nsd),
     &            cfl(nshg),                  icflhits(nshg),
     &            cflold(nshg)

c
        real*8    yl(npro,nshl,ndofl),        acl(npro,nshl,ndofl),
     &            xl(npro,nenl,nsd),         
     &            rl(npro,nshl),              ql(npro,nshl,nsd),
     &            dwl(npro,nenl),             evl(npro,nenl),
     &            cfll(npro,nshl),            cfllold(npro,nshl)
c        
        real*8    xSebe(npro,nshl,nshl),      xmudmi(npro,ngauss) 
c
c.... gather the variables
c
        real*8 sgn(npro,nshl)
c
c.... get the matrix of mode signs for the hierarchic basis functions. 
c
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
        
        call localy(y,      yl,     ien,    ndofl,  'gather  ')
        call localy(ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
        if(iRANS.lt. 0) 
     &  call localx(d2wall, dwl,    ien,    1,      'gather  ')
        call local (qres,   ql,     ien,    nsd,    'gather  ')

        if ((iDNS.gt.0).and.(itwmod.eq.-2)) then
          call local(effvisc, evl,    ien,    1,      'gather  ')
        endif

        if (iLSet.eq.2) then
c          call local(cfl, cfll, ien,  1,  'gather  ')
c          cfllold = cfll
          call local(cflold, cfllold, ien,  1,  'gather  ')
        endif
c
c.... zero the matrices if they are being recalculated
c
        if (lhs. eq. 1)  then
           xSebe = zero
        endif   
c
c.... get the element residuals, LHS matrix, and preconditioner
c
      rl = zero
      cfll = zero
      call e3Sclr  (yl,      acl,     shp,
     &              shgl,    xl,      dwl,
     &              rl,      ql,      xSebe,   
     &              sgn,     xmudmi,  cfll,
     &              cfllold, evl)
c
c.... assemble the residual
c
        call local (res,    rl,     ien,    1,  'scatter ')
c
c.... assemble the CFL values.  cfl will contain the sum of
c     all contributing integration points.  Will divide by
c     the number of contributors to get the average CFL number.
        if (iLSet.eq.2) then
          call localSum (cfl, cfll, ien, icflhits, 1)
c TEST 
	  cfl = cflold
	  cfll = cfllold
        endif
c
c.... end
c
        return
        end
