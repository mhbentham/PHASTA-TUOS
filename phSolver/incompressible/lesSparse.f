c---------------------------------------------------------------------
c     
c     drvftools.f : Bundle of Fortran driver routines for ftools.f
c     
c     Each routine is to be called by les**.c
c     
c---------------------------------------------------------------------
c     
c----------------
c     drvLesPrepDiag
c----------------
c     
      subroutine drvlesPrepDiag ( flowDiag, ilwork,
     &                            iBC,      BC,      iper,
     &                            rowp,     colm,    
     &                            lhsK,     lhsP)
c     
      use pointer_data
      use pvsQbi
      use convolImpFlow !brings in the current part of convol coef for imp BC
      include "common.h"
      include "mpif.h"
c     
      dimension flowDiag(nshg,4), ilwork(nlwork)
      dimension iBC(nshg), iper(nshg), BC(nshg,ndofBC)
      real*8 lhsK(9,nnz_tot), lhsP(4,nnz_tot)
      integer rowp(nnz_tot),  colm(nshg+1)
      integer	n,	k
c
      integer sparseloc
c
c     
c.... Clear the flowdiag
c
      if((flmpl.eq.1).or.(ipord.gt.1)) then
         do n = 1, nshg
	    k = sparseloc( rowp(colm(n)), colm(n+1)-colm(n), n )
     &       + colm(n)-1
c     
	    flowdiag(n,1) = lhsK(1,k)
	    flowdiag(n,2) = lhsK(5,k)
	    flowdiag(n,3) = lhsK(9,k)
c     
	    flowdiag(n,4) = lhsP(4,k)
         enddo
      else
	flowDiag = zero
	do n = 1, nshg  ! rowsum put on the diagonal instead of diag entry
           do k=colm(n),colm(n+1)-1

c
              flowdiag(n,1) = flowdiag(n,1) + abs(lhsK(1,k)) 
c     &                          + lhsK(2,k) + lhsK(3,k)
              flowdiag(n,2) = flowdiag(n,2) + abs(lhsK(5,k)) 
c     &                          + lhsK(4,k) + lhsK(6,k)
              flowdiag(n,3) = flowdiag(n,3) + abs(lhsK(9,k)) 
c     &                          + lhsK(7,k) + lhsK(8,k)
c
              flowdiag(n,4) = flowdiag(n,4) + abs(lhsP(4,k)) 
           enddo
           flowdiag(n,:)=flowdiag(n,:)*pt33
	enddo
      endif
      if(ipvsq.ge.3) then ! for first cut only do diagonal extraction
 ! this is not yet correct for multi procs I suspect if partition
 ! boundary cuts a p=QR face
         tfact=alfi * gami * Delt(1)
         do n=1,nshg
            if(numResistSrfs.gt.zero) then
               do k = 1,numResistSrfs
                  if (nsrflistResist(k).eq.ndsurf(n)) then
                     irankCoupled=k      
                     flowdiag(n,1:3) = flowdiag(n,1:3)
     &               + tfact*ValueListResist(irankCoupled)*
     &               NABI(n,:)*NABI(n,:)
                     exit
                  endif
               enddo
            elseif(numImpSrfs.gt.zero) then
               do k = 1,numImpSrfs
                  if (nsrflistImp(k).eq.ndsurf(n)) then
                     irankCoupled=k      
                     flowdiag(n,1:3) = flowdiag(n,1:3)
     &               + tfact*ImpConvCoef(ntimeptpT+2,irankCoupled)*
     &               NABI(n,:)*NABI(n,:)
                     exit
                  endif
               enddo
            endif
         enddo
      endif
c     

c
      if(iabc==1)    !are there any axisym bc's
     &      call rotabc(flowdiag, iBC, 'in ')
c

c     
c.... communicate : add the slaves part to the master's part of flowDiag
c     
	if (numpe > 1) then 
           call commu (flowDiag, ilwork, nflow, 'in ') 
        endif
c
c.... satisfy the boundary conditions on the diagonal
c
        call bc3diag(iBC, BC,  flowDiag)
c
c     
c.... on processor periodicity was not taken care of in the setting of the 
c     boundary conditions on the matrix.  Take care of it now.
c
        call bc3per(iBC,  flowDiag, iper, ilwork, 4)
c
c... slaves and masters have the correct values
c
c     
c.... Calculate square root
c     
        do i = 1, nshg
           do j = 1, nflow
              if (flowDiag(i,j).ne.0) 
     &             flowDiag(i,j) = 1. / sqrt(abs(flowDiag(i,j)))
           enddo
        enddo
c     
        return
        end

c     
c-------------
c     drvsclrDiag
c-------------
c     
      subroutine drvsclrDiag ( sclrDiag, ilwork, iBC, BC, iper, 
     &                         rowp,     colm,   lhsS )
c     
      use pointer_data
      include "common.h"
      include "mpif.h"
c     
      integer  ilwork(nlwork),    iBC(nshg),     iper(nshg),
     &         rowp(nnz_tot),    colm(nshg+1)

      real*8   sclrDiag(nshg),    lhsS(nnz_tot), BC(nshg,ndofBC)
      integer sparseloc

      sclrDiag = zero
      do n = 1, nshg
         k = sparseloc( rowp(colm(n)), colm(n+1)-colm(n), n ) 
     &                               + colm(n)-1
c
         sclrDiag(n) = lhsS(k)
      enddo
c     
c.... communicate : add the slaves part to the master's part of sclrDiag
c     
	if (numpe > 1) then 
           call commu (sclrDiag, ilwork, 1, 'in ') 
        endif
c
c.... satisfy the boundary conditions on the diagonal
c
        call bc3SclrDiag(iBC,  sclrDiag)
c
c     
c.... on processor periodicity was not taken care of in the setting of the 
c     boundary conditions on the matrix.  Take care of it now.
c
        call bc3per(iBC,  sclrDiag, iper, ilwork, 1)
c
c... slaves and masters have the correct values
c
c     
c.... Calculate square root
c     
        do i = 1, nshg
           if (sclrDiag(i).ne.0) then
              sclrDiag(i) = 1. / sqrt(abs(sclrDiag(i)))
           endif
        enddo
c     
      return
      end

C============================================================================
C
C "fLesSparseApG":
C
C============================================================================
	subroutine fLesSparseApG(	col,	row,	pLhs,	
     &					p,	q,	nNodes,
     &                                  nnz_tot )
c
c.... Data declaration
c
	implicit none
	integer	nNodes, nnz_tot
	integer	col(nNodes+1),	row(nnz_tot)
	real*8	pLhs(4,nnz_tot),	p(nNodes),	q(nNodes,3)
c
	real*8	pisave
	integer	i,	j,	k
c
c.... clear the vector
c
	do i = 1, nNodes
	    q(i,1) = 0
	    q(i,2) = 0
	    q(i,3) = 0
	enddo
c
c.... Do an AP product
c
	do i = 1, nNodes
c
	    pisave = p(i)
cdir$ ivdep
	    do k = col(i), col(i+1)-1
		j = row(k) 
c
		q(j,1) = q(j,1) - pLhs(1,k) * pisave
		q(j,2) = q(j,2) - pLhs(2,k) * pisave
		q(j,3) = q(j,3) - pLhs(3,k) * pisave
	    enddo
	enddo
c
c.... end
c
	return
	end

C============================================================================
C
C "fLesSparseApKG":
C
C============================================================================

	subroutine fLesSparseApKG(	col,	row,	kLhs,	pLhs,
     1					p,	q,	nNodes,
     2                                  nnz_tot_hide ) 
c
c.... Data declaration
c
c	implicit none
        use pvsQbi
        include "common.h"
	integer	nNodes, nnz_tot
	integer	col(nNodes+1),	row(nnz_tot)
	real*8	kLhs(9,nnz_tot),	pLhs(4,nnz_tot)
        real*8 	p(nNodes,4),	q(nNodes,3)
c
	real*8	tmp1,	tmp2,	tmp3,	pisave
	integer	i,	j,	k
c
c.... clear the vector
c
	do i = 1, nNodes
	    q(i,1) = 0
	    q(i,2) = 0
	    q(i,3) = 0
	enddo
c
c.... Do an AP product
c
	do i = 1, nNodes
c
	    tmp1 = 0
	    tmp2 = 0
	    tmp3 = 0
	    pisave   = p(i,4)
cdir$ ivdep
	    do k = col(i), col(i+1)-1
		j = row(k) 
		tmp1 = tmp1
     1		     + kLhs(1,k) * p(j,1)
     2		     + kLhs(4,k) * p(j,2)
     3		     + kLhs(7,k) * p(j,3)
		tmp2 = tmp2
     1		     + kLhs(2,k) * p(j,1)
     2		     + kLhs(5,k) * p(j,2)
     3		     + kLhs(8,k) * p(j,3)
		tmp3 = tmp3
     1		     + kLhs(3,k) * p(j,1)
     2		     + kLhs(6,k) * p(j,2)
     3		     + kLhs(9,k) * p(j,3)
c
		q(j,1) = q(j,1) - pLhs(1,k) * pisave
		q(j,2) = q(j,2) - pLhs(2,k) * pisave
		q(j,3) = q(j,3) - pLhs(3,k) * pisave
	    enddo
	    q(i,1) = q(i,1) + tmp1
	    q(i,2) = q(i,2) + tmp2
	    q(i,3) = q(i,3) + tmp3
	enddo

        if(ipvsq.ge.2) then
         tfact=alfi * gami * Delt(1)
           call ElmpvsQ(q,p,tfact)
        endif
c
c.... end
c
	return
	end


C============================================================================
C
C "fLesSparseApNGt":
C
C============================================================================

	subroutine fLesSparseApNGt(	col,	row,	pLhs,	
     1					p,	q,	nNodes,
     2                                  nnz_tot   )
c
c.... Data declaration
c
	implicit none
	integer	nNodes, nnz_tot
	integer	col(nNodes+1),	row(nnz_tot)
	real*8	pLhs(4,nnz_tot),	p(nNodes,3),	q(nNodes)
c
	real*8	tmp
	integer	i,	j,	k
c
c.... Do an AP product
c
	do i = nNodes, 1, -1
c
	    tmp = 0
	    do k = col(i), col(i+1)-1
		j = row(k)
c
		tmp = tmp
     1		    + pLhs(1,k) * p(j,1)
     2		    + pLhs(2,k) * p(j,2)
     3		    + pLhs(3,k) * p(j,3)
	    enddo
	    q(i) = tmp
	enddo
c
c.... end
c
	return
	end

C============================================================================
C
C "fLesSparseApNGtC":
C
C============================================================================

	subroutine fLesSparseApNGtC(	col,	row,	pLhs,	
     1					p,	q,	nNodes,
     2                                  nnz_tot )
c
c.... Data declaration
c
        implicit none
        integer	nNodes, nnz_tot
        integer	col(nNodes+1),	row(nnz_tot)
        real*8	pLhs(4,nnz_tot),	p(nNodes,4),	q(nNodes)
c
	real*8	tmp
	integer	i,	j,	k
c
c.... Do an AP product
c
	do i = nNodes, 1, -1
c
	    tmp = 0
	    do k = col(i), col(i+1)-1
		j = row(k)
c
		tmp = tmp
     1		    + pLhs(1,k) * p(j,1)
     2		    + pLhs(2,k) * p(j,2)
     3		    + pLhs(3,k) * p(j,3)
     4		    + pLhs(4,k) * p(j,4)
	    enddo
	    q(i) = tmp
	enddo
c
c.... end
c
	return
	end

C============================================================================
C
C "fLesSparseApFull":
C
C============================================================================

	subroutine fLesSparseApFull(	col,	row,	kLhs,	pLhs,
     1					p,	q,	nNodes,
     2                                  nnz_tot_hide )
c
c.... Data declaration
c
c	implicit none
        use pvsQbi
        include "common.h"

	integer	nNodes, nnz_tot
	integer	col(nNodes+1),	row(nnz_tot)
	real*8	kLhs(9,nnz_tot),	pLhs(4,nnz_tot)
        real*8  p(nNodes,4),	q(nNodes,4)
c
	real*8	tmp1,	tmp2,	tmp3,	tmp4,	pisave
	integer	i,	j,	k
c
c.... clear the vector
c
	do i = 1, nNodes
	    q(i,1) = 0
	    q(i,2) = 0
	    q(i,3) = 0
	enddo
c
c.... Do an AP product
c
	do i = 1, nNodes
c
	    tmp1 = 0
	    tmp2 = 0
	    tmp3 = 0
	    tmp4 = 0
	    pisave   = p(i,4)
cdir$ ivdep
	    do k = col(i), col(i+1)-1
		j = row(k)
c
		tmp1 = tmp1
     1		     + kLhs(1,k) * p(j,1)
     2		     + kLhs(4,k) * p(j,2)
     3		     + kLhs(7,k) * p(j,3)
		tmp2 = tmp2
     1		     + kLhs(2,k) * p(j,1)
     2		     + kLhs(5,k) * p(j,2)
     3		     + kLhs(8,k) * p(j,3)
		tmp3 = tmp3
     1		     + kLhs(3,k) * p(j,1)
     2		     + kLhs(6,k) * p(j,2)
     3		     + kLhs(9,k) * p(j,3)
c
		tmp4 = tmp4
     1		     + pLhs(1,k) * p(j,1)
     2		     + pLhs(2,k) * p(j,2)
     3		     + pLhs(3,k) * p(j,3)
     4		     + pLhs(4,k) * p(j,4)
c
		q(j,1) = q(j,1) - pLhs(1,k) * pisave
		q(j,2) = q(j,2) - pLhs(2,k) * pisave
		q(j,3) = q(j,3) - pLhs(3,k) * pisave
	    enddo
	    q(i,1) = q(i,1) + tmp1
	    q(i,2) = q(i,2) + tmp2
	    q(i,3) = q(i,3) + tmp3
	    q(i,4) = tmp4
	enddo
        if(ipvsq.ge.2) then
         tfact=alfi * gami * Delt(1)
           call ElmpvsQ(q,p,tfact)
        endif
c
c.... end
c
	return
	end

C============================================================================
C
C "fLesSparseApSclr":
C
C============================================================================

	subroutine fLesSparseApSclr(	col,	row,	lhs,	
     1					p,	q,	nNodes,
     &                                  nnz_tot)
c
c.... Data declaration
c
	implicit none
	integer	nNodes, nnz_tot
	integer	col(nNodes+1),	row(nnz_tot)
	real*8	lhs(nnz_tot),	p(nNodes),	q(nNodes)
c
	real*8	tmp
	integer	i,	j,	k
c
c.... Do an AP product
c
	do i = nNodes, 1, -1
c
	    tmp = 0
	    do k = col(i), col(i+1)-1
		tmp = tmp + lhs(k) * p(row(k))
	    enddo
	    q(i) = tmp
	enddo
c
c.... end
c
	return
	end

C============================================================================
	subroutine commOut(  global,  ilwork,  n, 
     &                       iper,    iBC, BC  )
	
	include "common.h"
	
	real*8  global(nshg,n), BC(nshg,ndofBC)
	integer ilwork(nlwork), iper(nshg), iBC(nshg)
c
	if ( numpe .gt. 1) then 
	   call commu ( global, ilwork, n, 'out')
	endif
c
c     before doing AP product P must be made periodic
c     on processor slaves did not get updated with the 
c     commu (out) so do it here
c
        do i=1,n
           global(:,i) = global(iper(:),i)  ! iper(i)=i if non-slave so no danger
        enddo
c
c       slave has masters value, for abc we need to rotate it
c        (if this is a vector only no SCALARS)
        if((iabc==1) .and. (n.gt.1)) !are there any axisym bc's
     &     call rotabc(global, iBC,  'out')


c$$$        do j = 1,nshg
c$$$           if (btest(iBC(j),10)) then
c$$$              i = iper(j)
c$$$              res(j,:) = res(i,:) 
c$$$           endif
c$$$        enddo
	
	return 
	end

C============================================================================
	subroutine commIn(  global,  ilwork,  n, 
     &                      iper,    iBC, BC )
	
	include "common.h"
	
	real*8  global(nshg,n), BC(nshg,ndofBC)
	integer ilwork(nlwork), iper(nshg), iBC(nshg)
c
        if((iabc==1) .and. (n.gt.1)) !are there any axisym bc's
     &       call rotabc(global, iBC, 'in ')
c

	if ( numpe .gt. 1 ) then
	   call commu ( global, ilwork, n, 'in ')
	endif
		
	call bc3per ( iBC, global, iper, ilwork, n)
	
	return 
	end

