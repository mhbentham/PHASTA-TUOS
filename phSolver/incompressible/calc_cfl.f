      subroutine calc_cfl (rho,          u1,       u2,
     &                     u3,           dxidx,    rmu,  
     &                     cfl_loc)  
c
c----------------------------------------------------------------------
c
c This routine computes the CFL for each element at an integration point. 
c
c input:
c  u1     (npro)           : x1-velocity component
c  u2     (npro)           : x2-velocity component
c  u3     (npro)           : x3-velocity component
c  dxidx  (npro,nsd,nsd)   : inverse of deformation gradient
c
c output:
c  cfl_loc(npro) 	   : CFL of the element
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension rho(npro),                 u1(npro),
     &            u2(npro),                  u3(npro),
     &            dxidx(npro,nsd,nsd), 
     &            rmu(npro), 
     &            dt2(npro), dt3(npro)
c
        dimension gijd(npro,6),  rnu(npro),  
     &            rhoinv(npro), cfl_loc(npro)
c
c.... get the metric tensor
c      
      call e3gijd( dxidx, gijd )

      rhoinv=one/rho
      rnu=rmu*rhoinv

       dt2 = ( u1 * ( gijd(:,1) * u1
     4		    + gijd(:,4) * u2
     5		    + gijd(:,6) * u3 )
     6	     + u2 * ( gijd(:,4) * u1
     7		    + gijd(:,2) * u2
     8		    + gijd(:,5) * u3 )
     9	     + u3 * ( gijd(:,6) * u1
     a		    + gijd(:,5) * u2
     1		    + gijd(:,3) * u3 ) ) 
       dt3 = rnu ** 2
     3	          * ( gijd(:,1) ** 2
     4	            + gijd(:,2) ** 2
     5		    + gijd(:,3) ** 2
     6		    + 2.
     7		  * ( gijd(:,4) ** 2
     8		    + gijd(:,5) ** 2
     9		    + gijd(:,6) ** 2 ) )
c
         if (ires == 1) then
c            cfl_loc= cfl_loc+sqrt(max(dt2,dt3/two))/(Dtgl*two)
            cfl_loc= cfl_loc+sqrt(dt2)/(Dtgl*two) !Farhad
         endif
c     
c.... return
c
        return
        end

