      subroutine genadj (colm,         rowp, icnt )
c     
      use pointer_data
c     
      include "common.h"
c     
      integer rowp(nshg*nnz),         colm(nshg+1)
      integer adjcnt(nshg),    row_fill_list(nshg,25*nnz), mloc(1)
c                                          change ^ if overflow
c                                   also change overflow check in asadj TWICE
      integer tmprdim(1), nnonzero
      real*8, allocatable, dimension(:) :: tmpr
c	write(*,*) '*********************'
      adjcnt=0
c      if (myrank.eq.master) write(*,*) ' point 1 '
 
      do iblk = 1, nelblk
c     
c.... set up the parameters
c     
         iel    = lcblk(1,iblk)
         lelCat = lcblk(2,iblk)
         lcsyst = lcblk(3,iblk)
         iorder = lcblk(4,iblk)
         nenl   = lcblk(5,iblk) ! no. of vertices per element
         nshl   = lcblk(10,iblk)
         npro   = lcblk(1,iblk+1) - iel 
         
c     
c.... compute sparse matrix data structures
c     
c      if (myrank.eq.master) write(*,*) ' point 2 ', iblk
         call Asadj (row_fill_list,                       
     &               mien(iblk)%p,  adjcnt )
         
      enddo
      
c      if (myrank.eq.master) write(*,*) ' point 3 '
      call sumgatInt ( adjcnt, nshg, nnonzero)
      if ( myrank .eq. master) then
         write (*,*) 'Number of global nonzeros ',nnonzero
      endif

c     
c     build the colm array
c     
      colm(1)=1
      do i=1,nshg
         colm(i+1)=colm(i)+adjcnt(i)
      enddo
c     
c     sort the rowp into increasing order
c     
      ibig=10*nshg
      icnt=0
      tmprdim=maxval(adjcnt)
      allocate (tmpr(tmprdim(1)))
      do i=1,nshg
         ncol=adjcnt(i)
         tmpr(1:ncol)=row_fill_list(i,1:ncol)
         do j=1,ncol
            icnt=icnt+1
            imin=minval(tmpr(1:ncol))
            mloc=minloc(tmpr(1:ncol))
            rowp(icnt)=imin
            tmpr(mloc(1))=ibig
         enddo
      enddo
      deallocate(tmpr)
      maxfill=tmprdim(1)
      if (myrank.eq.master) then
       write(*,*) 'maxfill=',maxfill
      end if
      nnza=icnt/nshg +1
      if (myrank.eq.master) then
      if(icnt.gt.nnz*nshg) then
         write(*,*) 'increase nnz in genmat to',nnza
         stop
      else
         write(*,*) 'nnz ok  nnz=',nnz,' actually needed',nnza   
         write(*,*) myrank,' is my rank and my nnz_tot is: ',icnt   
      endif
      endif
      return
      end









