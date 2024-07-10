        subroutine getelsize (x,  shp,  shgl,  elem_size,
     &                           shpb, shglb,  elemb_size,
     &                           elemvol_glob)
c
c----------------------------------------------------------------------
c
c 
c----------------------------------------------------------------------
c
        use pvsQbi  ! brings in NABI
        use stats   !  
        use pointer_data  ! brings in the pointers for the blocked arrays
        use local_mass
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
	real*8 elem_size(numel), elemb_size(numelb)
        integer numel_array(numpe), numel_tot, index(numpe)
	real*8, allocatable :: elem_size_tot(:), elem_vol_tot(:)
	real*8 elem_vol(numel), elemb_vol(numelb)
        real*8,intent(out):: elemvol_glob(numel)

c
c initialize array
c
        elem_size(:) = zero
        elem_vol(:) = zero
        numel_tot = 0
        elemvol_glob(:) = zero
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

!              if(myrank.eq.master)write(*,*)'iel+npro-1-iel ',iel+npro-1-iel,' npro ',npro
c
c compute local element volume
c
               call e3elsize (shp(lcsyst,1:nshl,:),
     &                        shgl(lcsyst,:,1:nshl,:),
     &         x,mien(iblk)%p,elem_size(iel:iel+npro-1),
     &         elem_vol(iel:iel+npro-1) )

               elemvol_glob(iel:iel+npro-1) = 
     &                  elem_vol(iel:iel+npro-1)

           enddo

c
c... Set Characteristic Size
c    i_spat_var_eps_flag = 0 - no spatially varying eps.  size = 1.0
c                          1 - based on vol**1/3
c                          2 - based on max edge length
c

        if (i_spat_var_eps_flag .eq. 0) then
          elem_size(:) = 1.0
        elseif (i_spat_var_eps_flag .eq. 1) then
          elem_size(:) = elem_vol(:)**(1.0/3.0)
        else
          elem_size(:) = elem_size(:)
        endif
cc
cc Open element volume file
cc
c        if (myrank .eq. master) then
c           open (unit=ivol,  file=fvol,  status='unknown')
c           write(ivol,*) "Processor   element   volume  size"
c       endif
cc
cc Gather number of elements for each processor --> store in
cc numel_array(numpe)
cc
c        if (numpe .gt. 1) then
c           call MPI_ALLGATHER(numel,1,MPI_INTEGER,
c     &         numel_array,1,MPI_INTEGER,
c     &         MPI_COMM_WORLD,ierr)
c          do j=1,numpe
c             numel_tot = numel_tot + numel_array(j)
c          enddo
c         else                   ! single processor
c            numel_array(1) = numel
c            numel_tot = numel
c         endif                  ! if-else f
c         write(*,*) "numel_tot = ", numel_tot
cc
cc Copy element volumes to global array
cc
c        allocate( elem_vol_tot(numel_tot) )
c        allocate( elem_size_tot(numel_tot) )
c        elem_vol_tot = zero
c        index(:) = 0
c        if (numpe .gt. 1) then
c          write(*,*) "proc, numel", myrank, numel
c          do iproc = 2,numpe
c            index(iproc) = index(iproc-1) + numel_array(iproc-1)
c          enddo
c          call MPI_ALLGATHERV(elem_vol,numel,MPI_DOUBLE_PRECISION,
c     &                        elem_vol_tot,numel_array,index,
c     &                        MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
c          call MPI_ALLGATHERV(elem_size,numel,MPI_DOUBLE_PRECISION,
c     &                        elem_size_tot,numel_array,index,
c     &                        MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
c        else
c           elem_vol_tot = elem_vol
c           elem_size_tot = elem_size
c        endif
cc
cc Write element data to file
cc
c        if (myrank .eq. master) then
c         do iproc = 1,numpe
c          do iel = 1, numel_array(iproc)
c           write(ivol,1000) iproc, iel, elem_vol_tot(iel+index(iproc)),
c     &                      elem_size_tot(iel+index(iproc))
c          enddo
c         enddo
c       endif
c 1000  format (i8,2x,i12,2x,e12.6,2x,e12.6)
c
c
c Now compute element volumes for boundary elements
c
c
c initialize array
c
        elemb_vol(:) = zero
        elemb_size(:) = zero	
c
c loop over element blocks to compute volume
c        
        do iblk = 1, nelblb
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel 
          ngauss = nint(lcsyst)
c
c compute local element volume
c
           call e3elsize (shp(lcsyst,1:nshl,:),
     &                    shgl(lcsyst,:,1:nshl,:),
     &         x,mienb(iblk)%p,elemb_size(iel:iel+npro-1),
     &         elemb_vol(iel:iel+npro-1) )

        enddo
c
c... Set Characteristic Size
c    i_spat_var_eps_flag = 0 - no spatially varying eps.  size = 1.0
c                          1 - based on vol**1/3
c                          2 - based on max edge length
c
        if (i_spat_var_eps_flag .eq. 0) then
          elemb_size(:) = 1.0
        elseif (i_spat_var_eps_flag .eq. 1) then
          elemb_size(:) = elemb_vol(:)**(1.0/3.0)
        else
          elemb_size(:) = elemb_size(:)
        endif
c
      return
      end
      
      
c..***********************************************************
c....this small routine just for calculate volume of the element
c*************************************************************
        subroutine e3elsize (shp,shgl,x,ien,loc_el_size,loc_el_vol)

c                                                                      
c----------------------------------------------------------------------
c This routine calculates characteristic size of each element 
c
c input: 
c  shp    (nen,ngauss) 		: element shape-functions
c  shgl   (nsd,nen,ngauss)	: element local-grad-shape-functions
c  x     (numnp,nsd)		: nodal coordinates at current step
c  ien (npro,nshl)
c output:
c  loc_el_vol(npro)		: volume of the element
c  loc_el_size(npro)		: maximum node-to-node distance for element
c  
c
c----------------------------------------------------------------------
c
        include "common.h"
c....Passed arrays
        dimension shp(nshl,ngauss), shgl(nsd,nshl,ngauss),
     &            x(numnp,nsd),ien(npro,nshl)

c
c local arrays
c
        real*8 loc_el_vol(npro)
	real*8 loc_el_size(npro)

        dimension dxidx(npro,nsd,nsd),WdetJ(npro)

        real*8 tmp(npro), disttmp(npro)
        dimension dxdxi(npro,nsd,nsd)

        dimension sgn(npro,nshl),  shape(npro,nshl),
     &  shdrv(npro,nsd,nshl),   xl(npro,nenl,nsd)

c
c Initialize array
c
        loc_el_vol(:) = zero
c
c Set sign
c
        do i=1,nshl
          where ( ien(:,i) < 0 )
            sgn(:,i) = -one
          elsewhere
            sgn(:,i) = one
          endwhere
        enddo
c
c*************Localizing coordinates*************
c
        call localx (x,xl,ien,nsd,'gather  ')
 
c
c.... loop through the integration points
c
        
        do intp = 1, ngauss   

          call getshp(shp,          shgl,      sgn, 
     &                shape,        shdrv)

c
c*************get WdetJ*********************************
c
c
c.... compute the deformation gradient
c
          dxdxi = zero
c
          do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shdrv(:,1,n)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shdrv(:,2,n)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shdrv(:,3,n)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shdrv(:,1,n)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shdrv(:,2,n)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shdrv(:,3,n)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shdrv(:,1,n)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shdrv(:,2,n)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shdrv(:,3,n)
          enddo
c
c.... compute the inverse of deformation gradient
c
          dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3) 
     &                   - dxdxi(:,3,2) * dxdxi(:,2,3)
          dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3) 
     &                   - dxdxi(:,1,2) * dxdxi(:,3,3)
          dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3) 
     &                   - dxdxi(:,1,3) * dxdxi(:,2,2)
          tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1) 
     &                         + dxidx(:,1,2) * dxdxi(:,2,1)  
     &                         + dxidx(:,1,3) * dxdxi(:,3,1) )
          dxidx(:,1,1) = dxidx(:,1,1) * tmp
          dxidx(:,1,2) = dxidx(:,1,2) * tmp
          dxidx(:,1,3) = dxidx(:,1,3) * tmp
          dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1) 
     &                  - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
          dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3) 
     &                  - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
          dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3) 
     &                  - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
          dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2) 
     &                  - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
          dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2) 
     &                  - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
          dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2) 
     &                  - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c
          WdetJ = Qwt(lcsyst,intp)/ tmp
c
          do i=1,nshl
            loc_el_vol(:) = loc_el_vol(:) + abs(shape(:,i)*WdetJ)
          enddo
c
c.... end of the loop over integration points
c
        enddo
c
c... Find maximum distance across nodes for element
c
        loc_el_size(:) = zero
        do n1 = 1, nenl-1
          do n2 = 2, nenl
            disttmp(:) = ( (xl(:,n1,1) - xl(:,n2,1))**2 +
     &                     (xl(:,n1,2) - xl(:,n2,2))**2 +
     &                     (xl(:,n1,3) - xl(:,n2,3))**2 )**0.5
            where (disttmp(:) .gt. loc_el_size(:))
              loc_el_size(:) = disttmp(:)
            endwhere
          enddo
        enddo
c
c      write(*,*) "e3elsize:elem vol value=",loc_el_vol(:)  

       return
       end

