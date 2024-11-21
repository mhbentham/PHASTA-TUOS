        subroutine ElmGMR (u,         y,         ac,
     &                     x,         banma,
     &                     shp,       shgl,      iBC,
     &                     BC,        shpb,      shglb,
     &                     res,       iper,      ilwork,
     &                     rowp,      colm,     lhsK,      
     &                     lhsP,      rerr,     GradV,
     &                     elemvol_global,
     &                     avgxcoordf,  avgycoordf,  avgzcoordf)

c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c Zdenek Johan, Winter 1991.      (Fortran 90)
c Chris Whiting, Winter 1998.     (Matrix EBE-GMRES)
c Alberto Figueroa, Winter 2004.  CMM-FSI
c Irene Vignon, Spring 2004.
c----------------------------------------------------------------------
c
        use pvsQbi      ! brings in NABI
        use stats       !  
        use pointer_data! brings in the pointers for the blocked arrays
        use local_mass
        use spat_var_eps
        use timedata    ! for iblkts usage
        use cf_arrays    ! Magnus, bubble controller access lift control arrays
        use bub_track   ! access to bubble information array
c
        include "common.h"
        include "mpif.h"
c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            u(nshg,nsd),
     &            x(numnp,nsd),               
     &            iBC(nshg),           
     &            BC(nshg,ndofBC),  
     &            res(nshg,nflow),
     &            iper(nshg)
c
        dimension banma(nshg,1)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qres(nshg,idflx),     rmass(nshg)
        dimension GradV(nshg,nsdsq)
c
        dimension ilwork(nlwork)

        integer rowp(nshg*nnz),         colm(nshg+1)

	real*8	lhsK(9,nnz_tot),	lhsP(4,nnz_tot)

        real*8, allocatable, dimension(:,:,:,:) :: xKebe, xGoC

        real*8  rerr(nshg,numerr)

        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

        real*8 spmasstot(20),  ebres(nshg)
        real*8 cfl(nshg), CFLfl_maxtmp
        integer icflhits(nshg)

! Matt T. Bubble Coalescence Control
        real*8 xarray(ibksiz), yarray(ibksiz), zarray(ibksiz)
        integer coordtag(ibksiz) !Passed arrays from e3ivar

        real*8 bubradius, bubradius2

        real*8 avgxcoordf(coalest), avgycoordf(coalest), avgzcoordf(coalest)
!----------------------------------------------------------------------
        real*8 elemvol_global(numel)
        real*8, allocatable::   elemvol_local(:)
c
c!.... set up the timer
c

CAD        call timer ('Elm_Form')
c
c!.... -------------------->   diffusive flux   <--------------------
c
c!.... set up parameters
c
        ires   = 1

        if (idiff==1 .or. idiff==3 .or. isurf==1) then ! global reconstruction
                                                       ! of qdiff
c
c! loop over element blocks for the global reconstruction
c! of the diffusive flux vector, q, and lumped mass matrix, rmass
c
           qres = zero
!MR CHANGE
           if(iter == nitr .and. icomputevort == 1 ) then
             !write(*,*) 'iter:',iter,' - nitr:',nitr
             !write(*,*) 'Setting GradV to zero'
             GradV = zero
           endif
!MR CHANGE END
           rmass = zero
        
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
c!.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass
!MR CHANGE
            if(iter == nitr .and. icomputevort == 1 ) then
               !write(*,*) 'Calling AsIqGradV'
               call AsIqGradV (y,                x,
     &                   shp(lcsyst,1:nshl,:), 
     &                   shgl(lcsyst,:,1:nshl,:),
     &                   mien(iblk)%p,
     &                   GradV)
             endif
!MR CHANGE END
              call AsIq (y,                x,                       
     &                   shp(lcsyst,1:nshl,:), 
     &                   shgl(lcsyst,:,1:nshl,:),
     &                   mien(iblk)%p,     mxmudmi(iblk)%p,  
     &                   qres,             rmass )
           enddo
c           if (myrank .eq. master) write(*,*) 'AsIq is done'
       
c
c!.... form the diffusive flux approximation
c
           call qpbc( rmass, qres, iBC, iper, ilwork )       
!MR CHANGE END
           if(iter == nitr .and. icomputevort == 1 ) then
             !write(*,*) 'Calling solveGradV'
             call solveGradV( rmass, GradV, iBC, iper, ilwork )
           endif
!MR CHANGE END
c
        endif 
c
c!.... -------------------->   interior elements   <--------------------
c
        res    = zero
        if (stsResFlg .ne. 1) then
           flxID = zero
           flxIDsclr = zero
        endif

        if (lhs .eq. 1) then
           lhsp   = zero
           lhsk   = zero
        endif
c
c! initailize array
c
        cfl = zero
        icflhits = 0
c
c!.... loop over the element-blocks
c
c!... arrays used in bubble tracking advanced analysis
        if(iBT .eq. 1) then

        allocate ( procs_dataset(i_num_bubbles, 22) )
        allocate ( procs_coordDn(i_num_bubbles, 3),
     &             procs_coordUp(i_num_bubbles, 3) )
        allocate ( Shear_NodeMin(i_num_bubbles, 2),
     &             Shear_NodeMax(i_num_bubbles, 2) )

        procs_dataset           = zero
        procs_coordDn(:,:)      =  1.0e10
        procs_coordUp(:,:)      = -1.0e10
        Shear_NodeMin(:,:)      =  1.0e10
        Shear_NodeMax(:,:)      = -1.0e10
        endif

!#################################################################
!MB, bubble controller
!       the following subroutine will calculate the control forces in x y, z
!       directions
        if(iCForz.eq.1)then
                CALL CFcalculator()
        end if !iCForz
!#################################################################


!Matt T. Initialize bubble coalecence control variables
        if (coalcon.eq.1) then
!           if (update_coalcon.eq.1) then
            xarray(:) = zero
            yarray(:) = zero
            zarray(:) = zero
            coordtag(:) = 0

            bubradius = 0.0d0

            avgxcoordf(:) = -1.0d3
            avgycoordf(:) = -1.0d3
            avgzcoordf(:) = -1.0d3

!           endif
        endif

        do iblk = 1, nelblk
          iblock = iblk         ! used in local mass inverse (p>2)
          iblkts = iblk         ! used in timeseries
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
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)

c!... assemble local element volume matrix and bubble information array
          if(iBT .eq. 1)  allocate( bub_info(npro,22) )
          if(iCForz.eq.1) allocate( cf_var(npro,15, i_num_bubbles)) ! MB, bubble controller
          allocate(elemvol_local(npro))
          elemvol_local(:)=elemvol_global(iel:iel+npro-1)
c
c.... allocate the element matrices
c
          allocate ( xKebe(npro,9,nshl,nshl) )
          allocate ( xGoC (npro,4,nshl,nshl) )
c
c.... compute and assemble the residual and tangent matrix
c
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
        
          !if (myrank.eq.master) write (*,*) 'calling AsIGMR'
          call AsIGMR (y,         ac,       banma,
     &                 x,         mxmudmi(iblk)%p,
     &                 tmpshp, tmpshgl, mien(iblk)%p,
     &                 res,    qres,
     &                 xKebe,  xGoC,        rerr,
     &                 cfl,    icflhits,    elemvol_local,
     &                 xarray, yarray,  zarray,
     &                 bubradius,  bubradius2,
     &                 coordtag)
        !if (myrank.eq.master) write (*,*) 'called AsIGMR'
!----------------------------------------------------------------------
!       collect the information from the blocks for andvanced analysis
        IF(iBT .eq. 1) then
           call BubASSY()
           deallocate ( bub_info )
        ENDIF !iBT
!----------------------------------------------------------------------

!######################################################################
! MB, required for the bubble controller
        IF(iCForz.eq.1)then
           call CFASSY()
           deallocate ( cf_var )
        ENDIF !iCForz
!######################################################################

c           if (myrank .eq. master) write(*,*) 'AsIGMR is done, iblk = ', iblk
c
c.... satisfy the BC''s on the implicit LHS
c     
          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
             if(ipord.eq.1) 
     &         call bc3lhs (iBC, BC,mien(iblk)%p, xKebe)  
             call fillsparseI (mien(iblk)%p, 
     &                 xKebe,            lhsK,
     &                 xGoC,             lhsP,
     &                 rowp,                      colm)
          endif
c           if (myrank .eq. master) write(*,*) 'fillsparse is done'

          deallocate ( elemvol_local )
          deallocate ( xKebe )
          deallocate ( xGoC  )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c
c!.... end of interior element loop
c
       enddo !iblk
!----------------------------------------------------------------------
!       The following MPI operations are used to assemble the dataset from
!       different processors together and then calculate the average bubble
!       information.                                    Jun, 2015
!----------------------------------------------------------------------
       IF (iBT .eq. 1) THEN
           call BubMPIprocess()
           deallocate ( procs_dataset )
           deallocate ( procs_coordDn, procs_coordUp )
           deallocate ( Shear_NodeMin, Shear_NodeMax )
        ENDIF !iBT
c!.... Matt Talley's Bubble Coalescence Control Subroutine
        if (coalcon.eq.1) then
!          if (update_coalcon.eq.1) then

           call CoalescCon (xarray, yarray, zarray, coordtag,
     &                     bubradius, bubradius2, avgxcoordf,
     &                     avgycoordf, avgzcoordf)
!          endif
        endif

! ##################################################################
!Magnus, bubble controller
c!... The MPI communications for control forces Jun, 2014 Fall
        if(iCForz.eq.1)then
                CALL CFMPIprocess()
        end if !iCForz
! ##################################################################

c$$$       if(ibksiz.eq.20 .and. iwrote.ne.789) then
c$$$          do i=1,nshg
c$$$             write(789,*) 'eqn block ',i 
c$$$             do j=colm(i),colm(i+1)-1
c$$$                write(789,*) 'var block',rowp(j)
c$$$
c$$$                do ii=1,3
c$$$                   write(789,111) (lhsK((ii-1)*3+jj,j),jj=1,3)
c$$$                enddo
c$$$             enddo
c$$$          enddo
c$$$          close(789)
c$$$          iwrote=789
c$$$       endif
c$$$ 111   format(3(e14.7,2x))
c$$$c
c.... add in lumped mass contributions if needed
c
       if((flmpr.ne.0).or.(flmpl.ne.0)) then
          call lmassadd(ac,res,rowp,colm,lhsK,gmass)
       endif

       have_local_mass = 1

c  Divide CFL number by the number of contributors to get
c  the average CFL number.
       cfl = cfl / icflhits

c
c Find element with worst (largest) CFL number
c
        CFLfl_max = cfl(1)
        iCFLfl_maxelem = 1
        do i = 1, nshg
          if (cfl(i) .gt. CFLfl_max) then
            CFLfl_max = cfl(i)
            iCFLfl_maxelem = i
          endif
        enddo
c	write(*,*) 'flow, myrank, cfl = ', myrank, CFLfl_max
c
c Communicate with other processors to get maximum CFL number over
c all processors
c
        if(numpe. gt. 1) then
c            if (myrank .eq. master) write(*,*) 'MPI_reduce is to start'
c     &   , MPI_MAX, MPI_SUM
           call MPI_ALLREDUCE (CFLfl_max, CFLfl_maxtmp, 1,
     &       MPI_DOUBLE_PRECISION,MPI_MAX, MPI_COMM_WORLD,ierr)  
c	  write(*,*) 'myrank, ierr = ', myrank, ierr 
c           if (myrank .eq. master) write(*,*) 'MPI reduce is done'
        else
           CFLfl_maxtmp = CFLfl_max
        endif
        CFLfl_max = CFLfl_maxtmp 
c           if (myrank .eq. master) write(*,*) 'CFLfl_max is done', CFLfl_max
c
c.... time average statistics
c       
       if ( stsResFlg .eq. 1 ) then

          if (numpe > 1) then
             call commu (stsVec, ilwork, nResDims  , 'in ')
          endif
          do j = 1,nshg
             if (btest(iBC(j),10)) then
                i = iper(j)
                stsVec(i,:) = stsVec(i,:) + stsVec(j,:)
             endif
          enddo
c     
          do i = 1,nshg
             stsVec(i,:) = stsVec(iper(i),:)
          enddo
c           if (myrank .eq. master) write(*,*) 'stsVec is done'

          if (numpe > 1) then
             call commu (stsVec, ilwork, nResDims  , 'out')
          endif
          return
          
       endif
c
c.... -------------------->   boundary elements   <--------------------
c
c           if (myrank .eq. master) write(*,*) 'starting b el in elmgmr' 

c
c.... Compute forces over boundary planes for convective
c     boundary conditions and normalization is on
c
        if ((ibcb_conv_p .eq. 1) .and. (ibcb_conv_p_norm .eq. 1)) then
         do iblk = 1, nelblb
c
c.... set up the parameters
c
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


          if(lcsyst.eq.3) lcsyst=nenbl
c
          if(lcsyst.eq.3 .or. lcsyst.eq.4) then
             ngaussb = nintb(lcsyst)
          else
             ngaussb = nintb(lcsyst)
          endif
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c
          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)

          call AsBForce (y,                       x,
     &                   tmpshpb,
     &                   tmpshglb,
     &                   mienb(iblk)%p,           mmatb(iblk)%p,
     &                   miBCB(iblk)%p,           mBCB(iblk)%p)

          deallocate (tmpshpb)
          deallocate (tmpshglb)
c
c.... end of boundary element loop
c
         enddo
c
c Compute average pressure over natural pressure boundary
c
         if (numpe > 1) then
            call MPI_ALLREDUCE (flxID(6,1), Atot, 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (flxID(7,1), Ftot, 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         else
            Ftot=flxID(7,1)
            Atot=flxID(6,1)
         endif
         presavg = Ftot / Atot             

         if (stsResFlg .ne. 1) then
           flxID = zero
           flxIDsclr = zero
         endif
c
        endif
c
c...................................................................
c...................................................................
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
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


          if(lcsyst.eq.3) lcsyst=nenbl
c
          if(lcsyst.eq.3 .or. lcsyst.eq.4) then
             ngaussb = nintb(lcsyst)
          else
             ngaussb = nintb(lcsyst)
          endif
c
c.... allocate the element matrices
c
          allocate ( xKebe(npro,9,nshl,nshl) )
          allocate ( xGoC (npro,4,nshl,nshl) )
          
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c
          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)

          call AsBMFG (u,                       y,
     &                 ac,                      x,
     &                 tmpshpb,
     &                 tmpshglb,
     &                 mienb(iblk)%p,           mmatb(iblk)%p,
     &                 miBCB(iblk)%p,           mBCB(iblk)%p,
     &                 res,                     xKebe)
c
c.... satisfy (again, for the vessel wall contributions) the BC''s on the implicit LHS
c
c.... first, we need to make xGoC zero, since it doesn''t have contributions from the 
c.... vessel wall elements

          xGoC = zero

          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
             if(ipord.eq.1)
     &         call bc3lhs (iBC, BC,mienb(iblk)%p, xKebe)
             call fillsparseI (mienb(iblk)%p,
     &                 xKebe,           lhsK,
     &                 xGoC,             lhsP,
     &                 rowp,                      colm)
          endif

          deallocate ( xKebe )
          deallocate ( xGoC )
          deallocate (tmpshpb)
          deallocate (tmpshglb)
c
c.... end of boundary element loop
c
       enddo

       if(ipvsq.ge.1) then
c
c....  pressure vs. resistance boundary condition sets pressure at
c      outflow to linearly increase as flow through that face increases
c      (routine is at bottom of this file)
c
          call ElmpvsQ (res,y,-1.0d0)     
       endif
           
c
c before the commu we need to rotate the residual vector for axisymmetric
c boundary conditions (so that off processor periodicity is a dof add instead
c of a dof combination).  Take care of all nodes now so periodicity, like
c commu is a simple dof add.
c
       if(iabc==1)              !are there any axisym bc's
     &       call rotabc(res, iBC,  'in ')
c
c
c.... -------------------->   communications <-------------------------
c

       if (numpe > 1) then
          call commu (res  , ilwork, nflow  , 'in ')
       endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3Res (iBC,  BC,  res,  iper, ilwork)
c
c.... return
c
c      call timer ('Back    ')
c           if (myrank .eq. master) write(*,*) 'bc3Res is done'

      return
      end


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!********************************************************************
!--------------------------------------------------------------------

      subroutine ElmGMRSclr (y,         ac,        x,     
     &                       shp,       shgl,      iBC,
     &                       BC,        shpb,      shglb,
     &                       res,       iper,      ilwork,
     &                       rowp,      colm,      lhsS,
     &                       cfl )
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c----------------------------------------------------------------------
c
        use pointer_data
        use local_mass
c
        include "common.h"
        include "mpif.h"
c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            x(numnp,nsd),         iBC(nshg),           
     &            BC(nshg,ndofBC),      res(nshg),
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qres(nshg,nsd),     rmass(nshg)
c
        integer ilwork(nlwork), rowp(nshg*nnz),   colm(nshg+1)

	real*8  lhsS(nnz_tot)

        real*8, allocatable, dimension(:,:,:) :: xSebe
        real*8  cfl(nshg), CFLls_maxtmp, cflold(nshg)
        integer icflhits(nshg)
c
c.... set up the timer
c

CAD        call timer ('Elm_Form')
c
c.... -------------------->   diffusive flux   <--------------------
c
        ires   = 1

        if (idiff==1 .or. idiff==3) then ! global reconstruction of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
           qres = zero
           rmass = zero
           icflhits = 0
           cfl = zero

           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lcsyst = lcblk(3,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              
              ngauss = nint(lcsyst)
c     
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

              call AsIqSclr (y,                   x,                       
     &                       shp(lcsyst,1:nshl,:), 
     &                       shgl(lcsyst,:,1:nshl,:),
     &                       mien(iblk)%p,     qres,                   
     &                       rmass, cfl, icflhits )
       
           enddo

c  Divide CFL number by the number of contributors to get
c  the average CFL number.
           cfl = cfl / icflhits

c
c.... form the diffusive flux approximation
c
           call qpbcSclr ( rmass, qres, iBC, iper, ilwork )       
c
        endif 
c
c.... -------------------->   interior elements   <--------------------
c
        res    = zero
        spmass = zero
        phvol  = zero

        if (lhs .eq. 1) then
           lhsS   = zero
        endif

        if ((impl(1)/10) .eq. 0) then   ! no flow solve so flxID was not zeroed
           flxID = zero
        endif

        icflhits = 0
        cflold = cfl
        cfl = zero
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
          iblock = iblk         ! used in local mass inverse (p>2)
          iel    = lcblk(1,iblk)
          lcsyst = lcblk(3,iblk)
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          ndofl  = lcblk(8,iblk)
          npro   = lcblk(1,iblk+1) - iel

          ngauss = nint(lcsyst)
c
c.... allocate the element matrices
c
          allocate ( xSebe(npro,nshl,nshl) )
c
c.... compute and assemble the residual and tangent matrix
c
          call AsIGMRSclr(y,                   ac,
     &                 x,
     &                 shp(lcsyst,1:nshl,:), 
     &                 shgl(lcsyst,:,1:nshl,:),
     &                 mien(iblk)%p,        res,
     &                 qres,                xSebe, mxmudmi(iblk)%p,
     &                 cfl,  icflhits, cflold )
c
c.... satisfy the BC''s on the implicit LHS
c     
          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
             call fillsparseSclr (mien(iblk)%p, 
     &                 xSebe,             lhsS,
     &                 rowp,              colm)
          endif

          deallocate ( xSebe )
c
c.... end of interior element loop
c
       enddo

c
c.... add in lumped mass contributions if needed
c
       if((flmpr.ne.0).or.(flmpl.ne.0)) then
          call lmassaddSclr(ac(:,isclr), res,rowp,colm,lhsS,gmass)
       endif

       have_local_mass = 1

c  Divide CFL number by the number of contributors to get
c  the average CFL number.
       cfl = cfl / icflhits

c
c Find element with worst (largest) CFL number
c
        CFLls_max = cfl(1)
        iCFLls_maxelem = 1
        do i = 1, nshg
          if (cfl(i) .gt. CFLls_max) then
            CFLls_max = cfl(i)
            iCFLls_maxelem = i
          endif
        enddo
c
c Communicate with other processors to get maximum CFL number over
c all processors
c
        if(numpe. gt. 1) then
           call MPI_ALLREDUCE (CFLls_max, CFLls_maxtmp, 1,
     &       MPI_DOUBLE_PRECISION,MPI_MAX, MPI_COMM_WORLD,ierr) 
c         write(*,*) 'myrank, ierr = ', myrank, ierr
        else
           CFLls_maxtmp = CFLls_max
        endif
        CFLls_max = CFLls_maxtmp
c
c
c  call DtN routine which updates the flux to be consistent with the
c  current solution values.  We will put the result in the last slot of
c  BC (we added a space in input.f).  That way we can localize this
c  value to the boundary elements.  This is important to keep from calling
c  the DtN evaluator more than once per node (it can be very expensive).
c
         if(idtn.eq.1)  call DtN(iBC,BC,y)

c
c Assemble the phasic volumes over all processors
c
c
c Compute average pressure over natural pressure boundary
c
       if (numpe > 1) then
         call MPI_ALLREDUCE (phvol(1), phvol1, 1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE (phvol(2), phvol2, 1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
       else
          phvol1=phvol(1)
          phvol2=phvol(2)
       endif
       phvol(1) = phvol1
       phvol(2) = phvol2

c
c.... -------------------->   boundary elements   <--------------------
c
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lcsyst = lcblkb(3,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel

          if(lcsyst.eq.3) lcsyst=nenbl
          if(lcsyst.eq.3 .or. lcsyst.eq.4) then
             ngaussb = nintb(lcsyst)
          else
             ngaussb = nintb(lcsyst)
          endif
c
c localize the dtn boundary condition
c

          if(idtn.eq.1)   call dtnl(   iBC, BC, mienb(iblk)%p,
     &              miBCB(iblk)%p,  mBCB(iblk)%p)

c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c
          call AsBSclr (y,                       x,
     &                  shpb(lcsyst,1:nshl,:),
     &                  shglb(lcsyst,:,1:nshl,:),
     &                  mienb(iblk)%p,           mmatb(iblk)%p,
     &                  miBCB(iblk)%p,           mBCB(iblk)%p,
     &                  res)
c
c.... end of boundary element loop
c
        enddo
c
c
c.... -------------------->   communications <-------------------------
c

      if (numpe > 1) then
        call commu (res  , ilwork, 1  , 'in ')
      endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3ResSclr (iBC,  res,  iper, ilwork)
c
c.... return
c
CAD      call timer ('Back    ')
      return
      end

        
c
c....routine to compute and return the flow rates for coupled surfaces of a given type
c        
      subroutine GetFlowQ (qsurf,y,srfIdList,numSrfs)
        
      use pvsQbi  ! brings in NABI
c
      include "common.h"
      include "mpif.h"

      real*8  y(nshg,3)
      real*8  qsurf(0:MAXSURF),qsurfProc(0:MAXSURF)
      integer numSrfs, irankCoupled, srfIdList(0:MAXSURF)

c note we only need the first three entries (u) from y

      qsurfProc=zero
      
      do i = 1,nshg
      
        if(numSrfs.gt.zero) then
          do k = 1,numSrfs
            irankCoupled = 0
            if (srfIdList(k).eq.ndsurf(i)) then
              irankCoupled=k
              do j = 1,3              
                 qsurfProc(irankCoupled) = qsurfProc(irankCoupled)
     &                            + NABI(i,j)*y(i,j)
              enddo
              exit
            endif      
          enddo       
        endif
      
      enddo
c      
c     at this point, each qsurf has its "nodes" contributions to Q
c     accumulated into qsurf. Note, because NABI is on processor this
c     will NOT be Q for the surface yet
c
c.... reduce integrated Q for each surface, push on qsurf
c
       npars=MAXSURF+1
       if(impistat.eq.1) iAllR = iAllR+1
       if(impistat2.eq.1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
       if(impistat.eq.1) rmpitmr = TMRC()
       call MPI_ALLREDUCE (qsurfProc, qsurf(:), npars,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
       if(impistat.eq.1) rAllR = rAllR+TMRC()-rmpitmr
  
c
c.... return
c
      return
      end


        
c
c... routine to couple pressure with flow rate for each coupled surface
c
      subroutine ElmpvsQ (res,y,sign)     

      use pvsQbi  ! brings in NABI
      use convolImpFlow !brings in the current part of convol coef for imp BC
c
      include "common.h"
      include "mpif.h"

      real*8 res(nshg,ndof),y(nshg,3)
      real*8 p(0:MAXSURF)
      integer irankCoupled

c
c... get p for the resistance BC
c           
      if(numResistSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistResist,numResistSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        p(:)=sign*p(:)*ValueListResist(:) ! p=QR  now we have the true pressure on each
                                        ! outflow surface.  Note sign is -1
                                        ! for RHS, +1 for LHS
c
c....  multiply it by integral NA n_i
c     
       do i = 1,nshg
          do k = 1,numResistSrfs
              irankCoupled = 0
              if (nsrflistResist(k).eq.ndsurf(i)) then 
                  irankCoupled=k
                  res(i,1:3)=res(i,1:3) + p(irankCoupled)*NABI(i,1:3)     
                  exit 
              endif
          enddo   
       enddo
       
      endif !end of coupling for Resistance BC

      
c
c... get p for the impedance BC
c     
      if(numImpSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistImp,numImpSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        do j = 1,numImpSrfs
            if(sign.lt.zero) then ! RHS so -1
                p(j)= sign*(pold(j) + p(j)*ImpConvCoef(ntimeptpT+2,j)) !pressure p=pold+ Qbeta
            elseif(sign.gt.zero) then ! LHS so sign is positive
                p(j)= sign*p(j)*ImpConvCoef(ntimeptpT+2,j)
            endif
        enddo
             
c
c....  multiply it by integral NA n_i
c     
       do i = 1,nshg
          do k = 1,numImpSrfs
              irankCoupled = 0
              if (nsrflistImp(k).eq.ndsurf(i)) then 
                  irankCoupled=k
                  res(i,1:3)=res(i,1:3) + p(irankCoupled)*NABI(i,1:3)      
                  exit
              endif
          enddo   
       enddo
       
      endif !end of coupling for Impedance BC
c
c.... return
c
      return
      end







