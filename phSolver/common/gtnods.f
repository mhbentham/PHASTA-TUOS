      subroutine gtnods

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

!SCATTER      dimension irecvcount(numpe), numvec(numpe)
      integer*8 numvec

      if(numpe > 1) then
         irecvcount = 1
         numvec = nshg0
c         call MPI_ALLREDUCE (numvec, nshgt, irecvcount,
c     &                 MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE (numvec, nshgt, irecvcount,
     &                 MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
c         call MPI_REDUCE_SCATTER (numvec, nshgt, irecvcount,
c     &                 MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      else
         nshgt = nshg0
      endif
     
      if (myrank .eq. master) then
         write(6,*) 'Total number of modes = ',nshgt
      endif
      
      return
      end
