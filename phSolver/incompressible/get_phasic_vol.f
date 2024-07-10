c..***********************************************************
c....this small routine just for calculate volume of the element
c*************************************************************
        subroutine get_phasic_vol(yl, shape, WdetJ)
c                                                                      
c----------------------------------------------------------------------
c This routine calculates characteristic size of each element 
c
c input: 
c  yl(npro,nshl,ndof)		: solution
c  shape  (npro, nshl)		: element shape-functions
c  WdetJ  (npro)		: Jacobian 
c output:
c  phvol(2)			: phasic volume (1-phase 1, 2-phase 2)
c
c----------------------------------------------------------------------
c
      include "common.h"
c....Passed arrays
      dimension shape(npro,nshl),
     &          yl(npro,nshl,ndof), WdetJ(npro)

c
c local arrays
c
      integer iel, n
      real*8 Sclr(npro)
c
c compute level set at gauss point
c
      Sclr = zero
      isc=abs(iRANS)+6
      do n = 1, nshl
        Sclr = Sclr + shape(:,n) * yl(:,n,isc)
      enddo

c
c compute volume at int. point and assing to one of the 
c Level Set phases
c
      do iel = 1, npro
        do n = 1, nshl  
          if (Sclr(iel) .gt. 0) then
            phvol(1) =  phvol(1) + abs(shape(iel,n)*WdetJ(iel))
          else
            phvol(2) =  phvol(2) + abs(shape(iel,n)*WdetJ(iel))
          endif
        enddo
      enddo
c
      return
      end

