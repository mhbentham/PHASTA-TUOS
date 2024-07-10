        subroutine check_proximity (y, stopjob)
c
c----------------------------------------------------------------------
c
c 
c----------------------------------------------------------------------
c
        use pointer_data  ! brings in the pointers for the blocked arrays
        use spat_var_eps   ! use spatially-varying epl_ls
c        
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        real*8        y(nshg,ndof)
	integer stopjob, itmp

        stopjob = 0             
c
c loop over element blocks to compute volume
c
           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lelCat = lcblk(2,iblk)
              lcsyst = lcblk(3,iblk)
              iorder = lcblk(4,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              nsymdl = lcblk(9,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              ngauss = nint(lcsyst)

c
c compute local element volume
c
               call elemchk (y,mien(iblk)%p,
     +              elem_local_size(lcblk(1,iblk):lcblk(1,iblk+1)-1),
     +              stopjob)

           enddo
c
c Communicate value of stopjob
c set to max value over all processors
c
        if (numpe > 1) then
          call MPI_ALLREDUCE (stopjob, itmp, 1, MPI_INTEGER,
     &         MPI_MAX, MPI_COMM_WORLD, ierr)
        else
          itmp = stopjob
        endif
        stopjob = itmp
c
      return
      end
      
      
c..***********************************************************
c....this small routine to check proximity of interface to larger elements
c*************************************************************
        subroutine elemchk (y,ien,loc_el_size,stopjob)

c                                                                      
c----------------------------------------------------------------------
c This routine 
c
c----------------------------------------------------------------------
c
        include "common.h"
c....Passed arrays
        dimension  y(nshg,ndof), ien(npro,nshl)
        real*8 loc_el_size(npro)
	integer stopjob

c
        dimension yl(npro,nshl,ndofl)
c
c*************Localizing solution*************
c	
        call localy(y,      yl,     ien,    ndofl,  'gather  ')
c
c Test if interface is near large element region
c Will test if edge of buffer region about interface exists in 
c larger elements. 
c    buffer region defined by: r_int_buf
c    large element size set by: r_int_elem_size

        i_flag = 0
c
        do ielem = 1, npro
	  if ((abs(yl(ielem,1,6)) .lt. r_int_buffer) .and. 
     &        (loc_el_size(ielem) .gt. r_int_elem_size)) then
            i_flag = 1
	  endif
	enddo
c
c Write message if interface has moved too close to large elements
c
        if (i_flag .eq. 1) then
          write(*,1000) time
 1000     format ("INTERFACE NEAR LARGE ELEMENTS: time = ",e12.6)
          stopjob = 1
c	  call MPI_BCAST(stopjob,1,MPI_INTEGER,myrank,
c     .               MPI_COMM_WORLD,0)
	endif
c
       return
       end

