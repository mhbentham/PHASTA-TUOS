      subroutine AsBForce (y,       x,       shpb,    shglb,
     &                     ienb,    materb,  iBCB,    BCB)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  boundary elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use turbSA                ! access to d2wall & effvisc
      include "common.h"
c
        dimension y(nshg,ndofl),           x(numnp,nsd),
     &            shpb(nshl,ngaussb),
     &            shglb(nsd,nshl,ngaussb),         
     &            ienb(npro,nshl), 
     &            iBCB(npro,ndiBCB),       BCB(npro,nshlb,ndBCB)       
c
        dimension yl(npro,nshl,ndofl),     xlb(npro,nenl,nsd),
     &            sgn(npro,nshl),
     &            ylotwn(npro,nshl,ndofl)

        real*8    evl(npro,nshl)
c
c.... get the matrix of mode signs for the hierarchic basis functions
c
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
c
c.... gather the variables
c
        call localy(y,      yl,     ienb,   ndofl,  'gather  ')
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')

        if ((iDNS.gt.0).and.(itwmod.eq.-2)) then
          call local(effvisc, evl,    ien,    1,      'gather  ')
        endif

c
c.... gather off-the-wall values for convective pressure boundaries
c
        ylotwn = zero
        if (any(btest(iBCB(:,1),1))) then
          call localyotwn(y, ylotwn, ienb,   ndofl,  'gather  ')
        endif
c
c.... 3D
c
        call e3bforce  (yl,      iBCB,    BCB,     shpb,    shglb,
     &                  xlb,      sgn,    ylotwn,  evl)

c     
c.... end
c
        return
        end
 
