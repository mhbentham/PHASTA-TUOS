       subroutine getgradphi (x, y,  shp,  shgl, 
     &                        shpb, shglb,  gradphi)
c
c----------------------------------------------------------------------
c
c 
c----------------------------------------------------------------------
c
        use pvsQbi  ! brings in NABI
        use stats
        use pointer_data  ! brings in the pointers for the blocked arrays
c
        include "common.h"
        include "mpif.h"
c
        dimension x(numnp,nsd)               
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &            shpb(MAXTOP,maxsh,MAXQPT),  
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)  
c
        dimension y(nshg,ndof)     
     
	real*8 gradphi(nshg,3)
	real*8 qres_tmp(nshg,idflx), rmass_tmp(nshg)
c
c initialize array
c
        gradphi = zero
	qres_tmp = zero
	rmass_tmp = zero
        idflow = (nflow-1)*nsd
c
c loop over element blocks to gradient of phi
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
c compute local gradient
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c              write(*,*) "getgradphi: calling Asiq for proc ", myrank,
c     &                   " block ", iblk
              call AsIq (y,                x,                       
     &                   shp(lcsyst,1:nshl,:), 
     &                   shgl(lcsyst,:,1:nshl,:),
     &                   mien(iblk)%p,     mxmudmi(iblk)%p,  
     &                   qres_tmp,          rmass_tmp )

           enddo
c
c.... invert the diagonal mass matrix and find q
c
        rmass_tmp = one/rmass_tmp
       
       do i=1, idflx
          qres_tmp(:,i) = rmass_tmp*qres_tmp(:,i)
       enddo
c
c retrieve phi gradient from qres array
c
c	write(*,*) "getgradphi: retrieving grad phi on proc", myrank
        gradphi(:,1) = qres_tmp(:,idflow+1)
	gradphi(:,2) = qres_tmp(:,idflow+2)
        gradphi(:,3) = qres_tmp(:,idflow+3)
c
c Output results to file
c name convention: gradphi.[timestep].[proc]
c

      return
      end
