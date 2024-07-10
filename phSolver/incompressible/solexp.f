      subroutine SolSclrExp(y,    ac,         yold,      
     &                   acold,      x,          iBC,
     &                   BC,         nPermDimsS,  nTmpDimsS,  
     &                   apermS,     atempS,     iper,       
     &                   ilwork,     shp,        shgl, 
     &                   shpb,       shglb,      rowp,     
     &                   colm,       lhsS,       solinc,
     &                   cfl )
c
c----------------------------------------------------------------------
c
c This is the forward Euler explicit scalar solve for IC flow
c
c input:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c  yold   (nshg,ndof)           : Y-variables at beginning of step
c  acold  (nshg,ndof)           : Primvar. accel. variable at begng step
c  x      (numnp,nsd)            : node coordinates
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  iper   (nshg)                : periodic nodal information
c
c output:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c
c
c
c----------------------------------------------------------------------
c
      use pointer_data
        
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
c     
      real*8    y(nshg,ndof),             ac(nshg,ndof),
     &          yold(nshg,ndof),          acold(nshg,ndof),
     &          x(numnp,nsd),             BC(nshg,ndofBC),
     &          res(nshg,1),
     &          flowDiag(nshg,4),
     &          sclrDiag(nshg,1),         lhsS(nnz_tot),
     &          apermS(nshg,nPermDimsS),  atempS(nshg,nTmpDimsS)


c
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &          shpb(MAXTOP,maxsh,MAXQPT),
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
      integer   usr(100),                 eqnType,
     &          rowp(nshg*nnz),           colm(nshg+1),
     &          iBC(nshg),                ilwork(nlwork),
     &          iper(nshg)
c
      real*8    yAlpha(nshg,ndof),        acAlpha(nshg,ndof),
     &          lesP(nshg,1),             lesQ(nshg,1),
     &          solinc(nshg,1),           cfl(nshg)
c
c.... form the LHS matrices, the residual vector (at n tilde)
c
      ac(:,7)=zero
      call ElmGMRSclr(y,      ac,         x,
     &             shp,       shgl,       iBC,       
     &             BC,        shpb,       shglb,
     &             res,       iper,       ilwork,   
     &             rowp,      colm,       lhsS,
     &             cfl )

c
c.... Remove any nodes outside "interface" from
c     redistancing step
c
      if (i_focusredist .gt. 0) then
      if ((iLSet.eq.2) .and. (isclr.eq.2)) then
        where (abs(y(:,5+isclr)) .gt. epsilon_lsd)
           res(:,1) = 0.0
        endwhere
      endif
      endif

c
c.... Compute acceleration
c
      solinc(:,1)=res(:,1)/gmass(:)
c
      if (numpe > 1) then
         call commu ( solinc, ilwork, 1, 'out')
      endif
      
      nsolsc=5+isclr
      call rstaticSclr (res, y, solinc, nsolsc) ! output scalar stats
c     
c.... end
c     
      return
      end





