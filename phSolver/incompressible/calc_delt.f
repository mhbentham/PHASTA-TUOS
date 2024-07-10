      subroutine calc_delt(istp)
c
c----------------------------------------------------------------------
c This routine modifies the time step based on the worst element cfl 
c number
c 
c----------------------------------------------------------------------
c
      include "common.h"
c
c modify time step if flag is set
c
      if ((iflag_cfl_dt .eq. 1) .and. (istp .gt. 1)) then
        if (CFLfl_max .le. zero) then
          write(*,*) "Zero velocity --> zero CFL - cannot modify delt"
        else
c
c compute scaling factor - not allowed to vary by more than 25% at a time
c
         factor = CFLfl_limit / CFLfl_max
         if (factor .lt. 0.75) factor = 0.75
         if (factor .gt. 1.25) factor = 1.25
         Delt(itseq) = Delt(itseq)*factor
         Dtgl = one / Delt(itseq)
         CFLfl_max = CFLfl_max * factor
        endif
      endif
c
      return
      end
      
      subroutine calc_deltau()
c
c----------------------------------------------------------------------
c This routine modifies the time step based on the worst element cfl
c number
c
c----------------------------------------------------------------------
c
      include "common.h"
      if (i_dtlset_cfl .gt. 0) then
         factor = dtlset_cfl / CFLls_max
          if (factor .lt. 0.75) factor = 0.75
          if (factor .gt. 1.25) factor = 1.25
         dtlset_new = dtlset * factor
         if (myrank.eq.master) then
	  if (factor.gt.1.01.or.factor.lt.0.99) then
           write (*,5001) dtlset, dtlset_new, CFLls_max, dtlset_cfl
          end if 
	 end if
 5001 format ("Pseudo time step for redistancing changed from ",
     &        e12.5," to ",e12.5," since max CFL of ",e12.5,
     &        " exceeds imposed limit of ",e12.5)
         dtlset = dtlset_new
         Delt(1) = dtlset ! psuedo time step for level set
         Dtgl = one / Delt(1)
c         CFLls_max = CFLls_max * factor
      endif
c
      return
      end

