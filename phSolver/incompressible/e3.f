        subroutine e3 (yl,      acl,     dwl,     shp,
     &                 shgl,    xl,      rl,      ql,
     &                 xKebe,   xGoC,    xmudmi,  sgn, 
     &                 rerrl,   rlsl,    cfll,    evl,
     &                 xwl,     ywl,     zwl,     xarray,
     &                 yarray,  zarray,  
     &                 elemvol_local,    bml,
     &                 bubradius, bubradius2, coordtag)
c                                                                      
c----------------------------------------------------------------------
c
c     This routine calculates the residual and tangent matrix for the 
c     UBar formulation of the incompressible Navier Stokes equations.
c
c
c input:    e    a   1..5   when we think of U^e_a  and U is 5 variables
c  yl     (npro,nshl,ndof)      : Y variables (not U)
c  acl    (npro,nshl,ndof)      : Y acceleration (Y_{,t})
c  shp    (nen,ngauss)           : element shape-functions  N_a
c  shgl   (nsd,nen,ngauss)       : element local-shape-functions N_{a,xi}
c  wght   (ngauss)               : element weight (for quadrature)
c  xl     (npro,nenl,nsd)       : nodal coordinates at current step (x^e_a)
c  ql     (npro,nshl,nsd*nsd) : diffusive flux vector (don't worry)
c  rlsl   (npro,nshl,6)       : resolved Leonard stresses
c  evl    (npro,nenl)         : effective viscosity for turb. wall function
c
c output:
c  rl     (npro,nshl,nflow)      : element RHS residual    (G^e_a)
c  rml    (npro,nshl,nflow)      : element modified residual  (G^e_a tilde)
c  xKebe  (npro,9,nshl,nshl)  : element LHS tangent mass matrix
c  xGoC   (npro,4,nshl,nshl)    : element LHS tangent mass matrix
c  cfll   (npro,nshl) 		: CFL of the element
c
c Note: This routine will calculate the element matrices for the
c        Hulbert's generalized alpha method integrator
c
c Mathematics done by:  Michel Mallet, Farzin Shakib (version 1)
c                       Farzin Shakib                (version 2)
c
c K. E. Jansen,   Winter 1999.   (advective form formulation)
c C. H. Whiting,  Winter 1999.   (advective form formulation)
c----------------------------------------------------------------------
c
        use spat_var_eps   ! use spatially-varying epl_ls
c
        include "common.h"
c
        dimension yl(npro,nshl,ndof), 
     &            acl(npro,nshl,ndof),       
     &            shp(nshl,ngauss),       shgl(nsd,nshl,ngauss),
     &            xl(npro,nenl,nsd),      dwl(npro,nenl),
     &            rl(npro,nshl,nflow),    ql(npro,nshl,idflx),
     &            cfll(npro,nshl),        evl(npro,nenl),
     &            xwl(npro,nenl),         ywl(npro,nenl),
     &            zwl(npro,nenl)
c      
        dimension bml(npro,nshl,1)
      	dimension xKebe(npro,9,nshl,nshl), xGoC(npro,4,nshl,nshl)

c
c.... Matt Talley's Bubble Coalescence Control
c
        real*8    xarray(ibksiz), yarray(ibksiz), zarray(ibksiz)
        real*8    bubradius, bubradius2

        integer   coordtag(ibksiz)

c
c.... local declarations
c
        dimension g1yi(npro,ndof),        g2yi(npro,ndof),
     &            g3yi(npro,ndof),        shg(npro,nshl,nsd),
     &            aci(npro,3),            dxidx(npro,nsd,nsd),       
     &            WdetJ(npro),            rho(npro),
     &            cp(npro),               k_T(npro),
     &            pres(npro),             u1(npro),
     &            u2(npro),               u3(npro),
     &            rLui(npro,nsd),         uBar(npro,nsd),
     &            xmudmi(npro,ngauss),     sgn(npro,nshl), 
     &            shpfun(npro,nshl),      shdrv(npro,nsd,nshl),
     &            rmu(npro),              tauC(npro),
     &            tauM(npro),             tauBar(npro),
     &            src(npro,3)

        dimension rlsl(npro,nshl,6),      rlsli(npro,6)

        real*8    rerrl(npro,nshl,6)
        integer   aa

c
c... needed for CFL calculation
c
      real*8 cfll_loc(npro)
      dimension Sclr(npro)

      real*8 epsilon_ls_tmp


c
c     
c.... local reconstruction of diffusive flux vector for quadratics
c     or greater but NOT for bflux since local mass was not mapped
c
        if ( idiff==2 .and. ires .eq. 1 ) then
           call e3ql (yl,        dwl,       shp,       shgl, 
     &                xl,        ql,        xmudmi, 
     &                sgn,       evl)
        endif
c
c.... loop through the integration points
c

        do intp = 1, ngauss

        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
c
c.... get the hierarchic shape functions at this int point
c
        call getshp(shp,          shgl,      sgn, 
     &              shpfun,       shdrv)
c
c.... get necessary fluid properties (including eddy viscosity)
c
        call getdiff(dwl,  yl, shpfun, xmudmi, xl, rmu, rho, cp, k_T,
     &               elem_local_size(lcblk(1,iblk):lcblk(1,iblk+1)-1),
     &               evl)
c
c.... calculate the integration variables
c
        call e3ivar (yl,          acl,       shpfun,
     &               shdrv,       xl,        bml,
     &               aci,         g1yi,      g2yi,    
     &               g3yi,        shg,       dxidx,   
     &               WdetJ,       rho,       pres, 
     &               u1,          u2,        u3,              
     &               ql,          rLui,      src,
     &               rerrl,       rlsl,      rlsli,
     &               dwl,         xwl,       ywl,
     &               zwl,         xarray,    yarray,
     &               zarray,      bubradius, bubradius2,
     &               coordtag,    elemvol_local) 
c
c.... compute CFL number
c
        cfll_loc = zero
        call calc_cfl(rho,          u1,       u2,
     &                u3,           dxidx,    rmu,    
     &                cfll_loc)

        do i=1,nshl
          cfll(:,i) = cfll(:,i) + shpfun(:,i)*cfll_loc
        enddo
c
c.... compute the stabilization terms
c
        call e3stab (rho,          u1,       u2,
     &               u3,           dxidx,    rLui,   
     &               rmu,          tauC,     tauM,   
     &               tauBar,       uBar )  
c
c.... compute the residual contribution at this integration point
c
       if (iLSet .eq. 2) then
          call get_phasic_vol(yl, shpfun, WdetJ)
       endif

       !if (myrank.eq.master) write (*,*) 'calling e3Res'
        call e3Res ( u1,        u2,         u3,
     &               uBar,      aci,        WdetJ,
     &               g1yi,      g2yi,       g3yi,
     &               rLui,      rmu,        rho,
     &               tauC,      tauM,       tauBar,
     &               shpfun,    shg,        src,
     &               rl,        pres,       acl,
     &               rlsli,      yl,         bml,
     &               elemvol_local)
        !if (myrank.eq.master) write (*,*) 'called e3Res'
c
c.... compute the tangent matrix contribution
c
        if (lhs .eq. 1) then
           call e3LHS ( u1,        u2,         u3,
     &                  uBar,      WdetJ,      rho,
     &                  rLui,      rmu,
     &                  tauC,      tauM,       tauBar,
     &                  shpfun,    shg,        xKebe,
     &                  xGoC )
        endif

c
c.... end of integration loop
c
      enddo
c
c.... symmetrize C
c
      if (lhs .eq. 1) then
         do ib = 1, nshl
            do iaa = 1, ib-1
               xGoC(:,4,iaa,ib) = xGoC(:,4,ib,iaa)
            enddo
         enddo
      endif
c
c.... return
c
      return
      end


c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c###################################################################


      subroutine e3Sclr (yl,      acl,     shp,
     &                   shgl,    xl,      dwl,
     &                   rl,      ql,      xSebe,   
     &                   sgn,     xmudmi,  cfll,
     &                   cfllold, evl)
c                                                                      
c----------------------------------------------------------------------
c
c     This routine calculates the residual and tangent matrix for the 
c     advection - diffusion equation for scalar.
c
c K. E. Jansen,   Winter 1999.   (advective form formulation)
c C. H. Whiting,  Winter 1999.   (advective form formulation)
c----------------------------------------------------------------------
c
      include "common.h"
c
      real*8    yl(npro,nshl,ndof),     acl(npro,nshl,ndof),       
     &          shp(nshl,ngauss),       shgl(nsd,nshl,ngauss),
     &          xl(npro,nenl,nsd),      rl(npro,nshl),          
     &          ql(npro,nshl,nsd),      xSebe(npro,nshl,nshl),
     &          dwl(npro,nshl),         cfll(npro,nshl),
     &          cfllold(npro,nshl),     evl(npro,nshl)
c
c.... local declarations
c
      real*8    gradS(npro,nsd),        shg(npro,nshl,nsd),
     &          Sdot(npro),             Sclr(npro),
     &          dxidx(npro,nsd,nsd),    WdetJ(npro),      
     &          u1(npro),     u2(npro), u3(npro),
     &          sgn(npro,nshl),         shpfun(npro,nshl),       
     &          shdrv(npro,nsd,nshl),   rLS(npro),
     &          tauS(npro),             diffus(npro),
     &          srcL(npro),             srcR(npro),
     &          gGradS(npro,nsd),       dcFct(npro),
     &          giju(npro,6)
c
c.... Source terms sometimes take the form (beta_i)*(phi,_i).  Since
c     the convective term has (u_i)*(phi,_i), it is useful to treat
c     beta_i as a "correction" to the velocity.  In calculating the
c     stabilization terms, the new "modified" velocity (u_i-beta_i) is 
c     then used in place of the pure velocity for stabilization terms,
c     and the source term sneaks into the RHS and LHS.
      real*8    uMod(npro,nsd), srcRat(npro), xmudmi(npro,ngauss)
c
      integer   aa, b

c
c... needed for CFL calculation
c
      real*8 rmu_tmp(npro), rho_tmp(npro), cfll_loc(npro)
      rmu_tmp = zero
      rho_tmp = 1.0  
c     
c.... local reconstruction of diffusive flux vector
c
        if ( idiff==2 ) then
           call e3qlSclr (yl, dwl, shp, shgl, xl, ql, sgn, evl)        
        endif
c
c.... loop through the integration points
c
        do intp = 1, ngauss

        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
c
c.... get the hierarchic shape functions at this int point
c
        call getshp(shp,          shgl,      sgn, 
     &              shpfun,        shdrv)
c
c.... get necessary fluid properties
c
        call getdiffsclr(shpfun,dwl,yl,diffus,evl)
c
c.... calculate the integration variables
c
        call e3ivarSclr(yl,          acl,       shpfun,
     &                  shdrv,       xl,        xmudmi,
     &                  Sclr,        Sdot,      gradS,
     &                  shg,         dxidx,     WdetJ,       
     &                  u1,          u2,        u3,              
     &                  ql,          rLS,       SrcR,
     &                  SrcL,        uMod,      dwl,
     &                  diffus,      srcRat,    evl,
     &                  cfllold )
c
c.... compute CFL number
c
        cfll_loc = zero
        call calc_cfl(rho_tmp,      u1,       u2,
     &                u3,        dxidx,    rmu_tmp,
     &                cfll_loc)

        do i=1,nshl
          cfll(:,i) = cfll(:,i) + shpfun(:,i)*cfll_loc
        enddo
c	if (myrank.eq.master) write(*,*) 'cfll(57) = ', cfll(57,:)
c
c.... compute the stabilization terms
c
        call e3StabSclr (uMod,    dxidx,   tauS, 
     &                   diffus,  srcR,    giju,
     &                   srcRat)
c
c... computing the DC factor for the discontinuity capturing
c
        if (idcsclr(1) .ne. 0) then
           if ((idcsclr(2).eq.1 .and. isclr.eq.1) .or. 
     &          (idcsclr(2).eq.2 .and. isclr.eq.2)) then ! scalar with dc
c
              call e3dcSclr ( gradS,    giju,     gGradS,
     &                        rLS,      tauS,     srcR,
     &                        dcFct)
           endif
        endif                   !end of idcsclr
c
c.... compute the residual contribution at this integration point
c
        call e3ResSclr ( uMod,      gGradS,
     &                   Sclr,      Sdot,       gradS,  
     &                   WdetJ,     rLS,        tauS,
     &                   shpfun,    shg,        srcR,
     &                   diffus, 
     &                   rl )
c
c.... compute the tangent matrix contribution
c
        if (lhs .eq. 1) then
           call e3LHSSclr ( uMod,      giju,       dcFct,
     &                      Sclr,      Sdot,       gradS,  
     &                      WdetJ,     rLS,        tauS,
     &                      shpfun,    shg,        srcL,
     &                      diffus,
     &                      xSebe )

        endif

c
c.... compute phasic volume for Level Set
c
       if (iLSet .eq. 2) then
          call get_phasic_vol(yl, shpfun, WdetJ)
       endif

c
c.... end of integration loop
c
      enddo

c
c.... return
c
      return
      end
