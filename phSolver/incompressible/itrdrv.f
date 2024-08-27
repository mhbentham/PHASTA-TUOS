      subroutine itrdrv (y,         ac,         banma,
     &                   uold,      x,         
     &                   iBC,       BC,         
     &                   iper,      ilwork,     shp,       
     &                   shgl,      shpb,       shglb,
     &                   ifath,     velbar,     nsons ) 
c
c!----------------------------------------------------------------------
c
c! This iterative driver is the semi-discrete, predictor multi-corrector 
c! algorithm. It contains the Hulbert Generalized Alpha method which
c! is 2nd order accurate for Rho_inf from 0 to 1.  The method can be
c! made  first-order accurate by setting Rho_inf=-1. It uses CGP and
c! GMRES iterative solvers.
c
c! working arrays:
c!  y      (nshg,ndof)           : Y variables
c!  x      (nshg,nsd)            : node coordinates
c!  banma  (nshg,1)              : marker field
c!  iBC    (nshg)                : BC codes
c!  BC     (nshg,ndofBC)         : BC constraint parameters
c!  iper   (nshg)                : periodicity table
c
c
c! Zdenek Johan,  Winter 1991.  (Fortran 90)
c! Alberto Figueroa, Winter 2004.  CMM-FSI
c! Irene Vignon, Fall 2004. Impedance BC
c! Jun Fang, Summer 2014. Bubble Tracking
c!----------------------------------------------------------------------
c
      use pvsQbi        !gives us splag (the spmass at the end of this run 
      use specialBC     !gives us itvn
      use timedata      !allows collection of time series
      use convolImpFlow !uses flow history and impedance for convolution
      use spat_var_eps  !use spatial varying eps_ls
      use bub_track     !access to bubble information array
      use turbsa        !access to d2wal info
      use redist_freeze !access to BC arrays if freezing value of primary LS vertices
      
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
        include "svLS.h"   ! MB, include svLS lib
c
        real*8 NewQImp(0:MAXSURF) !temporary unknown for the flow
                        !rate that needs to be added to the flow history
        
        real*8    y(nshg,ndof),              ac(nshg,ndof),           
     &	          yold(nshg,ndof),           acold(nshg,ndof),
     &            u(nshg,nsd),               uold(nshg,nsd),
     &            x(numnp,nsd),              solinc(nshg,ndof),
     &            BC(nshg,ndofBC),           tf(nshg,ndof),
     &            GradV(nshg,nsdsq)

c
!----------------------------------------------------------------------
        dimension banma(nshg,1)
        real*8, allocatable ::  breakupSeederTMP(:,:)
        real*8    elemvol_global(numel)   !Element volume array
        integer   ibreakupFlag, nBubbleLTS
        character*100 fbub
!----------------------------------------------------------------------
c
        real*8    res(nshg,ndof)
c     
        real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        integer   rowp(nshg,nnz),         colm(nshg+1),
     &            iBC(nshg),
     &            ilwork(nlwork),
     &            iper(nshg),            ifuncs(6)

        integer stopjob, iphase, recnum, numvar_rec
        character*10 cname2
        character*5  cname
        integer i_redist_counter
        real*8 redist_toler_previous
        logical iloop
c
c  stuff for dynamic model s.w.avg and wall model
c
        dimension ifath(numnp),    velbar(nfath,ndof),  nsons(nfath)

        dimension wallubar(2),walltot(2)
c     
c.... For Farzin''s Library
c
        integer eqnType, prjFlag, presPrjFlag, verbose
c
        real*8, allocatable, dimension(:,:) :: aperm,  atemp, atempS
        real*8, allocatable, dimension(:,:,:) :: apermS

        real*8, allocatable, dimension(:,:) :: lhsP, lhsK, lhsS
        real*8   almit, alfit, gamit
c
        character*1024    servername
        character*20    fname1,fmt1
        character*20    fname2,fmt2,fnamer2
        character*60    fnamepold, fvarts, fvartsb
        character*4     fieldybar
        integer         iarray(50) ! integers for headers
        
!MR CHANGD
        real*8, allocatable, dimension(:,:) :: strain, vorticity ! for vorticity
        real*8, allocatable, dimension(:,:) :: vorticitybar      ! for vorticity
        real*8, allocatable, dimension(:,:) :: wallssVec, wallssVecbar
        real*8 tcorecp(2)
!MR CHANGE END

        real*8 rerr(nshg,numerr),ybar(nshg,5) , uhess(nshg,27),
     &         gradu( nshg, 9 )
        integer, allocatable, dimension(:) :: ivarts
        integer, allocatable, dimension(:) :: ivartsg
        real*8, allocatable, dimension(:) :: vartssoln
        real*8, allocatable, dimension(:) :: vartssolng
        real*8 elem_size(numel), elemb_size(numelb)
        real*8 elem_size_min, elem_size_mintmp
      
        real*8 xi2
c
        real*8 gradphi(nshg,3), gradphimag(nshg), maxgradphi
	integer igradphi
c
c Redistancing option of fixing phi of primary vertices
c
        real*8  primvertval(nshg,2)
        integer primvert(nshg)
        integer i_primvert,numpv,numpvset
        integer  iredist_flag, numrun

        REAL*8,  allocatable :: BCredist(:)
        integer, allocatable :: iBCredist(:)

        real*8 CFLls(nshg)


c! Mass conservarion Epsilon_lsd adjustment:
	real*8 vf, vf0, eps_new, mass_error

c!....Matt Talley's Bubble Coal Control
        real*8 avgxcoordf(coalest), avgycoordf(coalest),
     &         avgzcoordf(coalest), avgxcoordold2(coalest),
     &         avgycoordold2(coalest), avgzcoordold2(coalest),
     &         app_time(coalest,2)

        real*8 itrtimestp

! ---------------------------------------------------------------------------
! MB block
!  Setting up the svLS
      integer svLS_nFaces_sc, svLS_nFaces, gnNo, nNo, faIn, facenNo
      integer svLS_nFacesT, gnNoT, nNoT, faInT, facenNoT
      integer, allocatable :: ltg(:), gNodes(:), gNodesT(:)
      real*8, allocatable :: sV(:,:), svT(:,:)
      
      character*128 fileName
      TYPE(svLS_commuType) communicator
      TYPE(svLS_lhsType) svLS_lhs, svLS_lhs_sc, svLS_lhsT
      TYPE(svLS_lsType) svLS_ls, svLS_sc, svLST 
      character*10 cname2nd
      real*8 sumtime
! ---------------------------------------------------------------------------
         if (myrank.eq.master) write(*,*) "starting the itrdrv subroutine..."
!----------------------------------------------------------------------
!       Initialize the current void fraction and constant used in
!       interface adjustment
!                                               Jun, March, 2014
!----------------------------------------------------------------------
        C_int_adjust    = 0.0d0
        vf_now          = 0.0d0
        if(myrank .eq. master) write(*,*) 'id2w is', id2w
!----------------------------------------------------------------------
!       Open files for bubble analysis
!----------------------------------------------------------------------
        if(iBT.eq.1 .and. iLSet .eq. 2) then
           if(myrank.eq.master) write(*,*)
           if(myrank.eq.master) write(*,*) 'Bubble tracking: Enabled!'
           if(iClrLiq.eq.1) then        !Specify the tracked liquid region
              phi_inner = phi_inner * epsilonBT
              phi_outer = phi_outer * epsilonBT
           else
              phi_inner = 0.0d0
              phi_outer = 0.0d0
           endif
!       get the total number of bubbles in the entire domain (i_num_bubbles)
           call CountIniBub(banma)
           allocate( bub_cent(i_num_bubbles,5) )
           ibreakupFlag         = 0
           if(iBK.eq.1) then
              allocate( breakupSeeder(i_num_bubbles,6) )
              breakupSeeder     = zero
              do i = 1, i_num_bubbles
                 breakupSeeder(i,1) = real(i)
                 breakupSeeder(i,6) = real(i)
              enddo
           endif
           nBubbleLTS           = i_num_bubbles
           ts_hold              = lstep
           GhostRatio           = 0.2
           NbrhdRatio           = 0.2

           if(myrank.eq.master) then
              call OpenBubFiles(nBubbleLTS, i_num_bubbles,
     &                          C_int_adjust)
              fbub = '../bubStat/bubble'
              fbub = trim(fbub)//trim(cname2(lstep))
              fbub = trim(fbub)//'.dat'
              fbub = trim(fbub)
              open(unit=20202, file=fbub, status="unknown",
     &        form="formatted", recl=2*8+13*14,
     &        access='direct')
              write(*,*) 'OpenBubFiles is done!'
           endif
           if(icoalCtrl.eq.1) then
              allocate(coalCenter(1,3))
              ncoalEvent = 0
           endif
        endif
!----------------------------------------------------------------------

!MR CHANGE - Scaling statistics
        impistat = 0
        impistat2 = 0
        iISend = 0
        iIRecv = 0
        iWaitAll = 0
        iAllR = 0
        rISend = zero
        rIRecv = zero
        rWaitAll = zero
        rAllR = zero
        rCommu = zero
!MR CHANGE END
!--------------------------------------------------------------------------------------
! MB block -- block added to initialise the svLS solver
      if(myrank.eq.master) write(*,*) "svLSFlag is set to ", svLSFlag
      IF (svLSFlag.EQ. 1) THEN
        if(myrank.eq.master) write(*,*) "calling svLS_LS_CREATE"
        call svLS_LS_CREATE(svLS_ls, LS_TYPE_GMRES, dimKry=Kspace,
     2   relTol=epstol(8), relTolIn=(/epstol(1), epstol(7)/),
     3   maxItr=nPrjs, maxItrIn=(/nGMRES, maxIters/))
         if(myrank.eq.master) write(*,*) "called svLS_LS_CREATE"

         nsolt = mod(impl(1), 2)
         nsclrsol=nsolt+nsclr

         if (nsclrsol.gt.0) then
            if(myrank.eq.master) write(*,*) "calling svLS_LS_CREATE"
            call svLS_LS_CREATE(svLS_ls, LS_TYPE_GMRES, dimKry=Kspace,
     2      relTol=epstol(8), relTolIn=(/epstol(1), epstol(7)/),
     3      maxItr=nPrjs, maxItrIn=(/nGMRES, maxIters/))
            if(myrank.eq.master) write(*,*) "called svLS_LS_CREATE"
         end if

         call svLS_COMMU_CREATE(communicator, MPI_COMM_WORLD)  ! MB, added to prevent MPI_ALLREDUCE problem

! Assuming the protocal to read the ltg files and set gnNo, nNO and ltg is not required
! we can simply comment all of this section out. This is what we are currently trying

         if (numpe.gt.1) then
            if (myrank.eq.master) write(*,*) "starting to write the ltg files"
            write(fileName,*) myrank+1
            fileName = "ltg.dat." //ADJUSTL(TRIM(fileName))
            if (numpe.gt.idirtrigger) then
               fileName = trim(cname2nd(int(myrank/dirstep)*idirstep))
     1         //"-set/"//trim(fileName)
            end if
            open(1, FILE=fileName)
            read(1,*) gnNo
            read(1,*) nNo
            allocate(ltg(nNo))
            read(1,*) ltg
            close(1)
            if (myrank.eq.master) write(*,*)"finished writing and reading the ltg files"
         else
! MB, I think this part of the code causes the problem with the bounds in svLS_LHS_CREATE
            ! memLS syntax
            gnNo = nshg
            nNo = nshg
            allocate(ltg(nNo))
            do i=1, nNo
               ltg(i) = i
            end do
         end if
            ! COLORADO SYNTAX
!            nNo = nshg
!            gnNo = nshgt
            


      ELSE 
!--------------------------------------------------------------------------------------
        if (myrank.eq.master) write(*,*) "svLSFlag is 0, so Acusim solver is used"
        call SolverLicenseServer(servername)

      END IF   ! MB
c
c only master should be verbose
c

        if(numpe.gt.0 .and. myrank.ne.master)iverbose=0  
c
c        call MPI_Barrier(MPI_COMM_WORLD)
c           if (myrank .eq. master) write(*,*) 'Starting itrdrv...' 

        inquire(file='xyzts.dat',exist=exts)
        lskeep = lstep 
        if(exts) then
           tssearch = 1
c           if (myrank .eq. master) write(*,*) 'Open the xyzts.dat...' 
           
           open(unit=626,file='xyzts.dat',status='old')
           read(626,*) ntspts, freq, tolpt, iterat, numvar, numrun 
           if (myrank .eq. master) write(*,*) 'numvar = ', numvar
           call sTD             ! sets data structures
           if (myrank .eq. master) write(*,*) 'sTD is done'
           

           do jj=1,ntspts       ! read coordinate data where solution desired
              read(626,*) ptts(jj,1),ptts(jj,2),ptts(jj,3)
           enddo
           close(626)
c           if (myrank .eq. master) write(*,*) 'ptts is read'

           statptts(:,:) = 0
           parptts(:,:) = zero
           varts(:,:) = zero

           allocate (ivarts(ntspts*numvar))
           allocate (ivartsg(ntspts*numvar))
           allocate (vartssoln(ntspts*numvar))
           allocate (vartssolng(ntspts*numvar))

c          if (myrank .eq. master) write(*,*) 'allocation is done'


           if (myrank .eq. master) then
                 fvartsb='varts/varts'
                 fvartsb=trim(fvartsb)//trim(cname2(lstep))
                 fvartsb=trim(fvartsb)//'.run'
                 fvartsb=trim(fvartsb)//trim(cname2(numrun))
                 fvartsb=trim(fvartsb)//'.dat'
                 fvartsb=trim(fvartsb)
              numvar_rec = max(15, numvar)  ! Create 19 variables for LS gradient case
                 open(unit=10101, file=fvartsb, status='unknown',
     1   form='formatted', recl=2*8+3+15*numvar_rec, access='direct')   ! Structured output file
           endif ! myrank

        endif   ! exts

c        call MPI_Barrier(MPI_COMM_WORLD)
c          if (myrank .eq. master) write(*,*) 'l.176, itvn: ', itvn

        if (itvn.gt.0.and.myrank.eq.master) then
          open(unit=1101, file='../bctinflow.dat', status='unknown')
          close(1101)
        end if 

c
c.... open history and aerodynamic forces files
c
        if (myrank .eq. master) then
           open (unit=ihist,  file=fhist, status='unknown')
           open (unit=iforce, file=fforce, status='unknown')
           open (unit=76, file="fort.76", status='unknown')
          write(*,*) ' l. 192; iLSet = ', iLSet
           if (iLSet .eq. 2) then
             open (unit=ivhist,  file=fvhist, status='unknown')
           endif
        endif
c
c.... initialize
c     
c         if (myrank .eq. master) write(*,*) 'l.201'
        ifuncs(:)  = 0              ! func. evaluation counter
        istep  = 0
        yold   = y
        acold  = ac

        rerr = zero
        ybar = zero
c        call MPI_Barrier(MPI_COMM_WORLD)
c          if (myrank .eq. master) write(*,*) 'l.210'



!MR CHANGE
        if(ivort == 1) then
          allocate(strain(nshg,6))
          allocate(vorticity(nshg,5))
          if(ioybar .eq. 1) then
            allocate(vorticitybar(nshg,4))
            vorticitybar=zero ! initialization for avg vorticity
          endif
        endif

        if(abs(itwmod).ne.1 .and. iowflux.eq.1) then
          allocate(wallssVec(nshg,3)) 
          if(ioybar .eq. 1) then
            allocate(wallssVecbar(nshg,3))
            wallssVecbar = zero ! Initialization important if mean wss computed
          endif
        endif
!MR CHANGE END

c
c.... ---------> initialize Farzin''s Library <---------------
c
c.... assign parameter values
c     
        do i = 1, 100
           numeqns(i) = i
        enddo
        nKvecs       = Kspace
        prjFlag      = iprjFlag
        presPrjFlag  = ipresPrjFlag
        verbose      = iverbose
c
c.... determine how many scalar equations we are going to need to solve
c
      nsolt=mod(impl(1),2)      ! 1 if solving temperature
      nsclrsol=nsolt+nsclr      ! total number of scalars solved At
                                ! some point we probably want to create
                                ! a map, considering stepseq(), to find
                                ! what is actually solved and only
                                ! dimension lhs to the appropriate
                                ! size. (see 1.6.1 and earlier for a
                                ! "failed" attempt at this).


      nsolflow=mod(impl(1),100)/10  ! 1 if solving flow
c        call MPI_Barrier(MPI_COMM_WORLD)
c          if (myrank .eq. master) write(*,*) 'l. 221 itrdrv.f'

c
c.... Now, call Farzin''s lesNew routine to initialize
c     memory space
c
c        call MPI_Barrier(MPI_COMM_WORLD)
c           if (myrank .eq. master) write(*,*) 'call genadj...', icnt
      call genadj(colm, rowp, icnt )  ! preprocess the adjacency list
c          if (myrank .eq. master) write(*,*) 'l. 230 itrdrv.f'

      nnz_tot=icnt ! this is exactly the number of non-zero blocks on
                   ! this proc

      if (nsolflow.eq.1) then
         lesId   = numeqns(1)
         eqnType = 1
         nDofs   = 4
!#########################################################################
! MB in this section we will initialse the linear solver parameters required for svLS
         ! method

         IF (svLSFlag .EQ. 1) THEN
            svLS_nFaces = 1
            write(*,*) 'myrank, gnNo =', myrank, gnNo
            !if (myrank.eq.master) write(*,*) 'calling svLS_BC_CREATE'
            !if (myrank.eq.master) write(*,*) 'gNodes is ', gNodes
            
            !if (myrank.eq.master) write(*,*) 'nNo is ', nNo
            !if (myrank .eq. master) write(*,*) 'calling svLS_LHS_CREATE'
            call svLS_LHS_CREATE(svLS_lhs, communicator, gnNo, nNo,
     2         nnz_tot, ltg, colm, rowp, svLS_nFaces)
            !if (myrank .eq. master) write(*,*) 'called svLS_LHS_CREATE'       

            faIn = 1
            facenNo = 0

            DO i=1, nshg
               IF (IBITS(iBC(i),3,3) .NE. 0) facenNo = facenNo + 1
            END DO

            allocate(gNodes(facenNo), sV(nsd, facenNo))
            sV = 0D0
            j = 0

            DO i=1, nshg
               IF (IBITS(iBC(i),3,3) .NE. 0) THEN
                  j = j + 1
                  gNodes(j) = i
                  IF (.NOT.BTEST(iBC(i),3)) sV(1,j) = 1D0
                  IF (.NOT.BTEST(iBC(i),4)) sV(2,j) = 1D0
                  IF (.NOT.BTEST(iBC(i),5)) sV(3,j) = 1D0
               END IF
            END DO
            

            call svLS_BC_CREATE(svLS_lhs, faIn, facenNo,
     2         nsd, BC_TYPE_Dir, gNodes, sV)
            if (myrank.eq.master) write(*,*) 'called svLS_BC_CREATE'

         ELSE  ! if svLS not available

!#########################################################################
            write(*,*) 'myrank, gnNo = ', myrank, gnNo
            if (myrank .eq. master) write(*,*) 'calling myfLesNew'
            call myfLesNew( lesId,   41994,
     &                 eqnType,
     &                 nDofs,          minIters,       maxIters,
     &                 nKvecs,         prjFlag,        nPrjs,
     &                 presPrjFlag,    nPresPrjs,      epstol(1),
     &                 prestol,        verbose,        statsflow,
     &                 nPermDims,      nTmpDims,      servername  )

            if (myrank .eq. master) write(*,*) 'called myfLesNew'

            allocate (aperm(nshg,nPermDims))
            allocate (atemp(nshg,nTmpDims))
         
         END IF      ! end of leslib vs svLS coondition ! MB

         allocate (lhsP(4,nnz_tot))
         allocate (lhsK(9,nnz_tot))

         if (myrank .eq. master) write(*,*) 'calling readLesRestart'
         
         call readLesRestart( lesId,  aperm, nshg, myrank, lstep,
     &                        nPermDims )
         if (myrank .eq. master) write(*,*) 'called readLesRestart'
      else
         nPermDims = 0
         nTempDims = 0
      endif

c          if (myrank .eq. master) write(*,*) 'l. 262 itrdrv.f'
!########################################################################
      if (myrank.eq.master) write(*,*) 'beginning solving for scalars block'
      if (nsclrsol.gt.0) then    ! solving for scalars
         IF (svLSFlag .EQ. 1) THEN
            ! Possible add separate memLS_lhs definition for each scalar ? Important for
            ! b.c.s
            
         if (nsolt.gt.0) then
               if (myrank.eq.master)
     &                  write(*,*) 'Setting BC for temperature'

         svLS_nFacesT = 1

         CALL svLS_LHS_CREATE(svLS_lhsT, communicator, gnNo, nNo,
     2         nnz_tot, ltg, colm, rowp, svLS_nFacesT)

         faInT = 1
         facenNoT = 0
         DO i=1, nshg
            IF (btest(iBC(i),1))  facenNoT = facenNoT + 1  ! Temperature check
         END DO
         ALLOCATE(gNodesT(facenNoT), sVT(1,facenNoT))
         sVT = 0D0
         j = 0
         DO i=1, nshg
            IF (btest(iBC(i),1)) THEN
               j = j + 1
               gNodesT(j) = i
               IF (.NOT.BTEST(iBC(i),1)) sV(1,j) = 1D0
            END IF
         END DO
         CALL svLS_BC_CREATE(svLS_lhsT, faInT, facenNoT,
     2         1, BC_TYPE_Dir, gNodesT, sVT)

      end if

      if (iLSet.gt.0) then
         svLS_nFaces_sc = 0

      CALL svLS_LHS_CREATE(svLS_lhs_sc, communicator, gnNo, nNo,
     2         nnz_tot, ltg, colm, rowp, svLS_nFaces_sc)

      if (myrank.eq.master) 
     &  write(*,*) 'Setting no-BC for level set and ls distance'

      end if

         ELSE     ! acusim

!########################################################################

      !if(nsclrsol.gt.0) then
       do isolsc=1,nsclrsol
         lesId       = numeqns(isolsc+1)
         eqnType     = 2
         nDofs       = 1
         presPrjflag = 0        
         nPresPrjs   = 0       
         prjFlag     = 1
         indx=isolsc+2-nsolt ! complicated to keep epstol(2) for
                             ! temperature followed by scalars
         call myfLesNew( lesId,            41994,
     &                 eqnType,
     &                 nDofs,          minIters,       maxIters,
     &                 nKvecs,         prjFlag,        nPrjs,
     &                 presPrjFlag,    nPresPrjs,      epstol(indx),
     &                 prestol,        verbose,        statssclr,
     &                 nPermDimsS,     nTmpDimsS,   servername )
       enddo
c
c  Assume all scalars have the same size needs
c
c          if (myrank .eq. master) write(*,*) 'l. 285 itrdrv.f'

       allocate (apermS(nshg,nPermDimsS,nsclrsol))
       allocate (atempS(nshg,nTmpDimsS))  !they can all share this

      END IF !end solver choice
      
      allocate (lhsS(nnz_tot,nsclrsol))
c
c actually they could even share with atemp but leave that for later
c
 
      else
         nPermDimsS = 0
         nTmpDimsS  = 0
      endif




c
c...  prepare lumped mass if needed
c
c      if((flmpr.ne.0).or.(flmpl.ne.0)) call genlmass(x, shp,shgl)
      call genlmass(x, shp, shgl, iBC, iper, ilwork)
c... compute element volumes
c
      allocate(elem_local_size(numel))
      if (numelb .gt. 0) then
        allocate(elemb_local_size(numelb))
      else
        allocate(elemb_local_size(1))
      endif
      call getelsize(x,  shp,  shgl,  elem_local_size,
     &               shpb, shglb,  elemb_local_size,
     &               elemvol_global)  


c
c! Can determine psuedo time step for redistancing now that 
c! element size is known
c
c!      if (i_dtlset_cfl .eq. 1) then
c!        elem_size_min = minval(elem_local_size)
c!        if(numpe. gt. 1) then
c!           call MPI_ALLREDUCE (elem_size_min, elem_size_mintmp, 1,
c!     &       MPI_DOUBLE_PRECISION,MPI_MIN, MPI_COMM_WORLD,ierr)
c!        else
c!           elem_size_mintmp = elem_size_min
c!        endif
c!        elem_size_min = elem_size_mintmp
c!        dtlset = dtlset_cfl * elem_size_min
c!      endif 
c
c! Initialize Level Set CFL array
c
       CFLls = zero
c
c! set flag for freezing LS Scalar 2
c
            if (iSolvLSSclr2.eq.2) then
              allocate (iBCredist(nshg))
              allocate (BCredist(nshg))
              BCredist(:) = zero
              iBCredist(:) = 0
            endif
c
c! Redistancing option of fixing phi of primary vertices
c
        i_primvert = 0
        if (i_primvert .eq. 1) then
          do inode = 1, nshg
           primvert(inode) = 0
           primvertval(inode,1) = zero 
           primvertval(inode,2) = zero
c           write(*,*) "primvertval, inode = ", primvertval(inode,1),
c     &                                  primvertval(inode,2), inode
          enddo
        endif
c           if (myrank .eq. master) write(*,*) 'initialization is done'

c
c!.... -----------------> End of initialization <-----------------
c

c!.... Initialize for Matt Talley's Coalescence Control
      if (coalcon.eq.1) then
         avgxcoordold(:) = -1.0d3
         avgycoordold(:) = -1.0d3
         avgzcoordold(:) = -1.0d3

         avgxcoordold2(:) = -1.0d3
         avgycoordold2(:) = -1.0d3
         avgzcoordold2(:) = -1.0d3

         app_time(:,2) = 0.0d0
         coalcon_rem(:) = 1
         itrtimestp = 0.0d0
      endif

c!.....open the necessary files to gather time series
c
c          if (myrank .eq. master) write(*,*) 'l. 360 itrdrv.f'

      lstep0 = lstep+1
c
c.... loop through the time sequences
c

      do 3000 itsq = 1, ntseq
         itseq = itsq

CAD         tcorecp1 = second(0)
CAD         tcorewc1 = second(-1)
c
c.... set up the time integration parameters
c         
         nstp   = nstep(itseq)
         nitr   = niter(itseq)
         LCtime = loctim(itseq)
         dtol(:)= deltol(itseq,:)

         call itrSetup ( y, acold )

c          if (myrank .eq. master) write(*,*) 'l. 382 itrdrv.f'



c IRENE: THIS SHOULD GO IN A SUBROUTINE
c...initialize the coefficients for the impedance convolution,
c   which are functions of alphaf so need to do it after itrSetup
         if(numImpSrfs.gt.zero) then
            allocate (ConvCoef(ntimeptpT+2,2)) !same time discret. for all imp. BC
            do j=1,ntimeptpT+2
                ConvCoef(j,:)=Delt(1)/2.0
            enddo
            ConvCoef(1,1)=ConvCoef(1,1)*(1.0-alfi)*(1.0-alfi)
            ConvCoef(1,2)=zero
            ConvCoef(2,2)=-ConvCoef(2,2)*(1.0-alfi)*(1.0-alfi)
            ConvCoef(ntimeptpT+1,1)=ConvCoef(ntimeptpT+1,1)*
     &                              alfi*(2.0-alfi)
            ConvCoef(ntimeptpT+2,2)=ConvCoef(ntimeptpT+2,2)*alfi*alfi
            ConvCoef(ntimeptpT+2,1)=zero  
            ConvCoef=ConvCoef/(ntimeptpT*Delt(1)) !divide by period T=N*dt
c
c...calculate the coefficients for the impedance convolution
c 
            allocate (ImpConvCoef(ntimeptpT+2,numImpSrfs))
            do j=2,ntimeptpT+1
                ImpConvCoef(j,:) = ValueListImp(j-1,:)*ConvCoef(j,2)
     &                             + ValueListImp(j,:)*ConvCoef(j,1)  
            enddo
            ImpConvCoef(1,:) = ValueListImp(1,:)*ConvCoef(1,1)
            ImpConvCoef(ntimeptpT+2,:) = 
     &           ValueListImp(ntimeptpT+1,:)*ConvCoef(ntimeptpT+2,2)
         endif
c
c  find the last solve of the flow in the step sequence so that we will
c         know when we are at/near end of step
c
c         ilast=0
c          if (myrank .eq. master) write(*,*) 'l. 417 itrdrv.f'

         nitr=0  ! count number of flow solves in a step (# of iterations)
         do i=1,seqsize
            if(stepseq(i).eq.0) nitr=nitr+1
         enddo

         if (numpe > 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!MR CHANGE
         tcorecp(:) = zero ! used in solfar.f (solflow)
!MR CHANGE END
         if(myrank.eq.0)  then
            tcorecp1 = TMRC()
         endif




c
c.... loop through the time steps
c
         istop=0
         rmub=datmat(1,2,1)
         if(rmutarget.gt.0) then
            rmue=rmutarget
         else
            rmue=datmat(1,2,1) ! keep constant
         endif
!         do 2000 istp = 1, nstp

        istp = 1
        do while (istp .le. nstp)
           if(iBT.eq.1 .and. myrank.eq.master)
     &        allocate( avg_info(i_num_bubbles,22) )


           if(iBT.eq.1 .and. iBK.eq.1 .and. ibreakupFlag.eq.1) then
              nBubbleLTS = i_num_bubbles
              call breakupConfirmer(y, banma)
              if(nBubbleLTS.lt.i_num_bubbles.and.myrank.eq.master)then
                 call OpenBubFiles(nBubbleLTS, i_num_bubbles,
     &                             C_int_adjust)
                 deallocate( avg_info )
                 allocate  ( avg_info(i_num_bubbles,22) )
              endif
           endif


           call rerun_check(stopjob)
           if(stopjob.ne.0) goto 2001
c
c Modify time step based on CFL number
c
           call calc_delt(istp)
c 
c            xi=istp*1.0/nstp
c            datmat(1,2,1)=rmub*(1.0-xi)+xi*rmue
c  ************ Time varying Viscosity and Density for second phase:
c            if (iramp.gt.0) then
c              xi2 = (real(lstep)-real(nrts0))/(real(nrts1)-real(nrts0))
c               if (lstep.lt.nrts0) xi2 = 0.0
c               if (lstep.gt.nrts1) xi2 = 1.0
            if (iramp.gt.0) then
              xi2 = (real(time)-real(qrts0))/(real(qrts1)-real(qrts0))
               if (time.lt.qrts0) xi2 = 0.0
               if (time.gt.qrts1) xi2 = 1.0
              datmat(1,2,2) = omu*(1.0-xi2)+xi2*tmu   ! 2nd phase viscosity 
              datmat(1,1,2) = orho*(1.0-xi2)+xi2*trho   ! 2nd phase density 
c              beta_rho = (trho/orho)**((time-qrts0)/(qrts1-qrts0)) 
c               beta_mu  =   (tmu/omu)**((time-qrts0)/(qrts1-qrts0))
c               if (time.gt.qrts0.and.time.lt.qrts1) then
c                  datmat(1,2,2) = omu*beta_mu
c                  datmat(1,1,2) = orho*beta_rho	
c               end if
c               if (time.ge.qrts1) then
c                   datmat(1,2,2) = tmu
c                   datmat(1,1,2) = trho
c               end if

                if(myrank.eq.master) then
                 if (xi2.eq.0.0) then
                   write(*,*) "Ramping properties will start at time = ", qrts0
                 else if (xi2.eq.1.0) then
c                  write(*,*) "Ramping properties is done at time = ", qrts1
                 else
                   write(*,*) "Second rho/mu",
     &    " is changed to ", datmat(1,1:2,2), " xi = ", xi2 
                 end if
                endif
            end if     
c **************************** above lines are added on 02/11/2009 by Igor Bolotnov; updated on 03/26/2009 - geometric ramp is introduced
c.... if we have time varying boundary conditions update the values of BC.
c     these will be for time step n+1 so use lstep+1
c     
c           if(itvn.gt.0) call BCint((lstep)*Delt(1), shp, shgl, 
c    &                               shpb, shglb, x, BC, iBC)
            if(itvn.gt.0) call BCint(time,shp, shgl, 
     &                               shpb, shglb, x, BC, iBC)
            if(itvbc.gt.0) call setBC(time, x, iBC, BC)
c
c ... calculate the pressure contribution that depends on the history 
c     for the impendance BC
c
c          if (myrank .eq. master) write(*,*) 'l. 489 itrdrv.f'


            if(numImpSrfs.gt.0) call pHist(pold,QHistImp,ImpConvCoef,
     &                                          ntimeptpT,numImpSrfs)
c
c Decay of scalars
c
           if(nsclr.gt.0 .and. tdecay.ne.1) then
              yold(:,6:ndof)=y(:,6:ndof)*tdecay
              BC(:,7:6+nsclr)= BC(:,7:6+nsclr)*tdecay
           endif




           if(nosource.eq.1) BC(:,7:6+nsclr)= BC(:,7:6+nsclr)*0.8


            if(iLES.gt.0) then  !complicated stuff has moved to
                                        !routine below
               call lesmodels(yold,  acold,     shgl,      shp, 
     &                        iper,  ilwork,    rowp,      colm,
     &                        nsons, ifath,     x,   
     &                        iBC,   BC)

            
            endif

c.... set traction BCs for modeled walls
c
            if (itwmod.ne.0) then
               call asbwmod(yold,   acold,   x,      BC,     iBC,
     &                      iper,   ilwork,  ifath,  velbar)
            endif
c          if (myrank .eq. master) write(*,*) 'l. 521 itrdrv.f'
!MR CHANGE
c
c.... Determine whether the vorticity field needs to be computed for this time step or not
c
            icomputevort = 0
            if (ivort == 1) then ! Print vorticity = True in solver.inp
              ! We then compute the vorticity only if we 
              ! 1) we write an intermediate checkpoint
              ! 2) we reach the last time step and write the last checkpoint
              ! 3) we accumulate statistics in ybar for every time step
              ! BEWARE: we need here lstep+1 and istep+1 because the lstep and 
              ! istep gets incremened after the flowsolve, further below
              if (((irs .ge. 1) .and. (mod(lstep+1, ntout) .eq. 0)) .or.
     &                   istep+1.eq.nstep(itseq) .or. ioybar == 1) then
                icomputevort = 1
              endif
            endif

!            write(*,*) 'icomputevort: ',icomputevort, ' - istep: ',
!     &                istep,' - nstep(itseq):',nstep(itseq),'- lstep:',
!     &                lstep, '- ntout:', ntout
!MR CHANGE END


c
c.... -----------------------> predictor phase <-----------------------
c
c           if (myrank .eq. master) write(*,*) 'predictor phase'
            call itrPredict(yold, y,   acold,  ac ,  uold,  u)


            call itrBC (y,  ac, banma, iBC,  BC,  iper,ilwork)
c           if (myrank .eq. master) write(*,*) 'l. 529 itrdrv.f'
            if(nsolt.eq.1) then
               isclr=0
               call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
            endif
            do isclr=1,nsclr
               call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
            enddo
            iter=0
            ilss=0  ! this is a switch thrown on first solve of LS redistance
             

                
c
c set the initial tolerance for the redistance loop
c
            if (i_redist_loop_flag.eq.1) then
              redist_toler_previous = 100.0
            endif
c
c LOOP OVER SEQUENCES
c
            istepc = 1
            iloop = .true.
            i_redist_counter=0



c            do istepc=1,seqsize
             do while (iloop) 
               icode=stepseq(istepc)
c               write(*,'(i5)') icode
c               
               if(mod(icode,5).eq.0) then ! this is a solve
c                write(*,'(a)') 'icode is a multiplier of 5'
                  isolve=icode/10
                  if(icode.eq.0) then ! flow solve (encoded as 0)
c
                     iter   = iter+1
                     ifuncs(1)  = ifuncs(1) + 1
c     
                     Force(1) = zero
                     Force(2) = zero
                     Force(3) = zero
                     HFlux    = zero
                     entrop   = zero
                     lhs = 1 - min(1,mod(ifuncs(1)-1,LHSupd(1))) 
c             if (myrank .eq. master) write(*,*) 'SolFlow is about to be called'

                     avgxcoordf(:) = -1.0d3
                     avgycoordf(:) = -1.0d3
                     avgzcoordf(:) = -1.0d3

      if (myrank.eq.master) write(*,*) 'ndof=',ndof
      if (myrank.eq.master) write(*,*) 'ndof2=',ndof2
                     call SolFlow(y,          ac,        u,
     &                         banma,
     &                         yold,          acold,     uold,
     &                         x,             iBC,
     &                         BC,            res,
     &                         nPermDims,     nTmpDims,  aperm,
     &                         atemp,         iper,          
     &                         ilwork,        shp,       shgl,
     &                         shpb,          shglb,     rowp,     
     &                         colm,          lhsK,      lhsP,
     &                         solinc,        rerr,      tcorecp,
     &                         GradV,         elemvol_global,
     &                         avgxcoordf, avgycoordf, avgzcoordf,
     &                         svLS_lhs, svLS_ls, svLS_nFaces)


c!....Matt Talley's Coalescence Contorl
                      if (coalcon.eq.1) then
                         if (coaltimtrak.eq.1) then

                            do k = 1, coalest
                               avgxcoordold(k) = avgxcoordf(k)
                               avgycoordold(k) = avgycoordf(k)
                               avgzcoordold(k) = avgzcoordf(k)

!                               if (avgxcoordold(k).gt.-1.0d3) then
!                                  if (myrank.eq.master) write(*,*) 'Coalescence',
!     &                            ' Event #: ', k
!                                  if (myrank.eq.master) write(*,*) 'x average',
!     &                            ' position:', avgxcoordold(k)
!                                  if (myrank.eq.master) write(*,*) 'y average',
!     &                            ' position:', avgycoordold(k)
!                                  if (myrank.eq.master) write(*,*) 'z average',
!     &                            ' position:', avgzcoordold(k)
!                               endif

                            enddo ! k

                         else

                            itrtimestp = Delt(1)

                            call CoalescAppTime (avgxcoordf, avgycoordf,
     &                                          avgzcoordf, avgxcoordold2,
     &                                          avgycoordold2, avgzcoordold2,
     &                                          app_time, itrtimestp)
                         endif ! coaltimtrak
                      endif ! coalcon

!MR CHANGE END
                  else          ! scalar type solve
                     if (icode.eq.5) then ! Solve for Temperature
                                ! (encoded as (nsclr+1)*10)

                        isclr=0
                        ifuncs(2)  = ifuncs(2) + 1
                        j=1

                     else       ! solve a scalar  (encoded at isclr*10)
                        isclr=isolve  
                        ifuncs(isclr+2)  = ifuncs(isclr+2) + 1
                        j=isclr+nsolt
c  Modify psuedo time step based on CFL number for redistancing
                        if((iLSet.eq.2).and.(ilss.ge.1).and.
     &                     (i_dtlset_cfl.eq.1).and.
     &                     (isclr.eq.2)) then
                           call calc_deltau()
                           Delt(1) = dtlset ! psuedo time step for level set
                           Dtgl = one / Delt(1)
                           ilss = ilss+1
                        endif
c
                        if((iLSet.eq.2).and.(ilss.eq.0)
     &                       .and.(isclr.eq.2)) then 
                           ilss=1 ! throw switch (once per step)
                           y(:,7)=y(:,6) ! redistance field initialized
                           ac(:,7)   = zero
!23456789012345678901234567890123456789012345678901234567890123456789012
                           if (iSolvLSSclr2.eq.2)  then
                             call get_bcredist(x,y,iBCredist,BCredist,
     &                                      primvert, primvertval(:,1))
                             primvertval(:,1) = BCredist(:)
                             ib=5+isclr
                             ibb=ib+1
                             do inode = 1, nshg
                              if (iBCredist(inode).eq.1) then
                               if (btest(iBC(inode),ib)) then
                              write(*,*) "WARNING -- Bit 7 already set"
                               endif
                              endif
                             enddo
c
                             where (iBCredist.eq.1) 
                              iBC(:) = iBC(:) + 128   ! set scalar 2 (bit 7)    
                              BC(:,ibb) = BCredist(:)
                             endwhere
                             numpv = 0
                             numpvset = 0
                             do inode = 1, nshg
                               if (primvert(inode) .gt. 0) then
                               numpv = numpv + 1
                                 if (primvert(inode).eq.2) then
                                   numpvset = numpvset + 1
                                 endif
                               endif
                             enddo
                             write(*,*) lstep+1,
     &                                  " Primary Verts: set/exist = ",
     &                                  numpvset, numpv
                           endif
                           call itrBCSclr (  y,  ac,  iBC,  BC, iper,
     &                          ilwork)


c     
c....store the flow alpha, gamma parameter values and assigm them the 
c....Backward Euler parameters to solve the second levelset scalar
c     
                           alfit=alfi
                           gamit=gami
                           almit=almi
                           Deltt=Delt(1)
                           Dtglt=Dtgl
                           alfi = 1
                           gami = 1
                           almi = 1
c     Delt(1)= Deltt ! Give a pseudo time step
                           Delt(1) = dtlset ! psuedo time step for level set
                           Dtgl = one / Delt(1)
                        endif  ! level set eq. 2
                     endif ! deciding between temperature and scalar
                        


                     lhs = 1 - min(1,mod(ifuncs(isclr+2)-1,
     &                                   LHSupd(isclr+2))) 

                     if((isclr.eq.1.and.iSolvLSSclr1.eq.1) .or. 
     &                  (isclr.eq.2.and.iSolvLSSclr2.eq.1) .or.
     &                  (isclr.eq.2.and.iSolvLSSclr2.eq.2)) then
                      
                      write(*,'(a)') 'This big block is executed'
                      lhs=0
                      call SolSclrExp(y,          ac,        yold,
     &                         acold,         x,         iBC,
     &                         BC,            nPermDimsS,nTmpDimsS,  
     &                         apermS(1,1,j), atempS,    iper,          
     &                         ilwork,        shp,       shgl,
     &                         shpb,          shglb,     rowp,     
     &                         colm,          lhsS(1,j), 
     &                         solinc(1,isclr+5), CFLls)
                     else
!                     write(*,'(a)') 'This second big block is executed'
!                      write(*,*)istepc

                     call SolSclr(y,          ac,        u,
     &                         yold,          acold,     uold,
     &                         x,             iBC,
     &                         BC,            nPermDimsS,nTmpDimsS,  
     &                         apermS(1,1,j), atempS,    iper,          
     &                         ilwork,        shp,       shgl,
     &                         shpb,          shglb,     rowp,     
     &                         colm,          lhsS(1,j), 
     &                         solinc(1,isclr+5), CFLls)
                        
                                
                     endif

                        
                  endif         ! end of scalar type solve

               else ! this is an update  (mod did not equal zero)
                  iupdate=icode/10  ! what to update
!                  write(*,*)'update',istepc

                  if(icode.eq.1) then !update flow  
                     call itrCorrect ( y,    ac,    u,   solinc)
                     call itrBC (y,  ac, banma, iBC,  BC, iper, ilwork)
                  else  ! update scalar
                     isclr=iupdate  !unless
!                     write(*,*)'update',istepc
                     if(icode.eq.6) isclr=0
                     if(iRANS.lt.0)then  ! RANS
                        call itrCorrectSclrPos(y,ac,solinc(1,isclr+5))
                     else
                        call itrCorrectSclr (y, ac, solinc(1,isclr+5))
                     endif
                     if (ilset.eq.2 .and. isclr.eq.2)  then
                        if (ivconstraint .eq. 1) then
                           call itrBCSclr (  y,  ac,  iBC,  BC, iper,
     &                          ilwork)
c                    
c ... applying the volume constraint on second level set scalar
c
                           call solvecon (y,    x,      iBC,  BC, 
     &                          iper, ilwork, shp,  shgl)
c
                        endif   ! end of volume constraint calculations
                     endif      ! end of redistance calculations
c                           write(*,*)istepc
                          call itrBCSclr (  y,  ac,  iBC,  BC, iper,
     &                       ilwork)

c
c ... update the old value for second level set scalar
c
                     if (ilset.eq.2 .and. isclr.eq.2)  then
                         call itrUpdateDist( yold, acold, y, ac)
                     endif   

                     endif      ! end of flow or scalar update
                  endif         ! end of switch between solve or update
c
c** Conditions for Redistancing Loop **
c Here we test to see if the following conditions are met:
c	no. of redistance iterations < i_redist_max_iter
c	residual (redist_toler_curr) > redist_toler
c If these are true then we continue in the redistance loop
c
                 if(i_redist_loop_flag.eq.1) then
                   if (icode .eq. 21) then ! only check after a redistance update
                     if((ilset.eq.2).and.(isclr.eq.2)) then !redistance condition
       if (redist_toler_curr.gt.redist_toler.or.
     &  i_redist_counter.eq.0) then !condition 1
              if (i_redist_counter.lt.i_redist_max_iter) then ! condition 2
                        i_redist_counter = i_redist_counter + 1
                        istepc = istepc - 2  ! repeat the 20 21 step
                        if(redist_toler_curr.gt.redist_toler_previous)
     &                  then
c                         if(myrank.eq.master) then
c                          write(*,*) "Warning: diverging!"
c                         endif
                        endif
                       else
                        iloop = .false. 
                        if(myrank.eq.master) then  
                         write(*,*) "Exceeded Max # of the iterations: "
     &                              , i_redist_max_iter
                        endif
                       endif
                       redist_toler_previous=redist_toler_curr
                      else
                       if(myrank.eq.master) then
                        write(*,*) "Redistance loop converged in ",
     &                       i_redist_counter," iterations"
                       endif
                       iloop = .false. 
                      endif
                     endif
                   endif !end of the redistance condition
                 endif !end of the condition for the redistance loop
c
                 if (istepc .eq. seqsize) then
                   iloop = .false.
                 endif
                 istepc = istepc + 1
             
c
c**End of loop condition for Redistancing equation**
c		 		  
               end do      ! end while loop over sequence steps
c
c Check if interface has moved into region of larger interface
c
             if ((iLSet.eq.2).and.(i_check_prox.eq.1)) then 
               call check_proximity(y, stopjob)
	       if(stopjob.ne.0) then
                   lstep = lstep + 1
                   goto 2001
               endif
             endif          
c
c print out phasic volume for Level Set
c
!        open(unit=737,file='vftarget.dat',status='old')
!        read(737,*) vf0
!        close(737)

        if (iLSet.eq.2) then
           vf =  phvol(2)/(phvol(1)+phvol(2))
           if(istp.eq.1)then
              vf_obj = vf_target 
              if(myrank.eq.master)then
                 write(*,*)'Target/Initial vf is ',vf_obj
              end if
           end if
           C_int_adjust = C_int_adjust - (vf - vf_obj)*vfcontrcoeff

           if(iuse_vfcont_cap .eq. 1)then
              if(C_int_adjust .gt. C_int_cap)then
                 C_int_adjust = C_int_cap
              elseif(C_int_adjust.lt.-C_int_cap)then
                 C_int_adjust = -C_int_cap
              endif
           end if

           if (myrank.eq.master) then 
               write (ivhist,800) lstep+1, phvol(1), 
     &  phvol(2), vf, epsilon_lsd
               write (*,801) lstep+1, vf*1.0D+02, 
     &  (vf/vf_target-1.0D+00)*1.00D+02
           endif 

        endif      !iLSet.eq.2
 800    format(1x,i6,4e15.6)
 801    format(i6,' Void Fraction: ',F15.10, 
     1 '%; Mass Error: ', F15.10,'%')
     
c
c.... obtain the time average statistics
c
c	if (myrank.eq.master) write(*,*) 'GetStats is about to run' 
            if (ioform .eq. 2) then

               call stsGetStats( y,      yold,     ac,     acold,
     &                           u,      uold,     x,
     &                           shp,    shgl,     shpb,   shglb,
     &                           iBC,    BC,       iper,   ilwork,
     &                           rowp,   colm,     lhsK,   lhsP )
            endif

c     
c  Find the solution at the end of the timestep and move it to old
c
c  
c ...First to reassign the parameters for the original time integrator scheme
c
            if((iLSet.eq.2).and.(ilss.gt.0)) then 
               alfi =alfit
               gami =gamit
               almi =almit 
               Delt(1)=Deltt
               Dtgl =Dtglt
            endif          
            call itrUpdate( yold,  acold,   uold,  y,    ac,   u)
            call itrBC (yold, acold, banma, iBC,  BC,  iper,ilwork)

c        if (myrank.eq.master) write(*,*) 'itrBC is done'

!----------------------------------------------------------------------
!       print out bubble information and some other processing at the end of
!       each time iteration
!----------------------------------------------------------------------
        if(iBT.eq.1 .and. i_num_bubbles.ne.0) then
           Nbubtot    = 0          !total number of bubbles
           Nghost     = 0          !number of ghost bubbles
           deallocate(bub_cent)
           if(myrank.eq.master) call reCountBub()
           if(numpe > 1) then
              call MPI_Bcast(Nbubtot,1,MPI_INTEGER,
     &                    master,MPI_COMM_WORLD,ierr)
              call MPI_Bcast(Nghost, 1,MPI_INTEGER,
     &                    master,MPI_COMM_WORLD,ierr)
              call MPI_Bcast(i_num_bubbles, 1,MPI_INTEGER,
     &                    master,MPI_COMM_WORLD,ierr)
           endif
           allocate(bub_cent(i_num_bubbles+Nghost,5))
           bub_cent = zero
!       breakupSeeder will be reshaped if total number of bubbles changes
           if(iBK.eq.1 .and. nBubbleLTS .lt. i_num_bubbles) then
              allocate(breakupSeederTMP(nBubbleLTS,6))
              breakupSeederTMP(1:nBubbleLTS,:) =
     &           breakupSeeder(1:nBubbleLTS,:)
              deallocate(breakupSeeder)
              allocate(breakupSeeder(i_num_bubbles,6))
              breakupSeeder(:,:) = 0
                 breakupSeeder(1:nBubbleLTS,:) =
     &        breakupSeederTMP(1:nBubbleLTS,:)
              deallocate(breakupSeederTMP)
              nBubbleLTS = i_num_bubbles
           endif

           if(myrank.eq.master) then
!       save the bubble curvature for the detection of bubble breakup in the
!       next time step
              if(iBK.eq.1) call breakupDetector(ibreakupFlag)
!       write out bubble information, and prepare the bubble center information
!       (including real bubbles and ghost bubbles at periodic faces)
              write(*,'(1x,A,I5)') 'Nbubtot =', Nbubtot
              write(*,'(1x,A,I5)') 'Nghost  =', Nghost
              call BubPrintOut(vf)
              do i=1,i_num_bubbles
                 if(avg_info(i,4).gt.0.0d0)
     &           write(20202,'(2I8, 22ES14.4)')
     &           lstep+1, i, avg_info(i,1:22)
              enddo

              write(*,*) 'BubPrintOut is done!'
              deallocate(avg_info)
           endif
!       broadcast the bubbles centers and radii for the marker field updating
!       algorithm
           if(numpe > 1) then
              call MPI_Bcast(bub_cent,(i_num_bubbles+Nghost)*5,
     &             MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
              if(iBK.eq.1) call MPI_Bcast(ibreakupFlag,1,
     &             MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
              if(iBK.eq.1) call MPI_Bcast(breakupSeeder,i_num_bubbles*6,
     &             MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
           endif
        endif     !iBT
!----------------------------------------------------------------------
!       Jun Fang Coalescence control based on bubble tracking information
!----------------------------------------------------------------------
        if(icoalCtrl.eq.1) then
           deallocate( coalCenter )
!       Obtain the centers of coalescence events
           if(myrank.eq.master.and.i_num_bubbles.gt.1)
     &     call coalescenceDetection()
           if(numpe > 1) then
              call MPI_Bcast(ncoalEvent,1,MPI_INTEGER,
     &                    master,MPI_COMM_WORLD,ierr)
              if(myrank.ne.master) allocate(coalCenter(ncoalEvent,3))
              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
              call MPI_Bcast(coalCenter,ncoalEvent*3,
     &             MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
           endif
        endif !end of coalescence part


            istep = istep + 1
            lstep = lstep + 1
            time  = time + delt(itseq)
c
c... Post process data and write out
c
            if ((i_gradphi.eq.1) .and. (mod(lstep, ntout).eq.0)) then
	      idflx = 0 
              if(idiff >= 1 )  idflx= (nflow-1) * nsd
              if (isurf == 1) idflx=nflow*nsd
	      call getgradphi(x, y, shp,       shgl,
     &                         shpb,          shglb, gradphi)
              gradphimag(:) = ( gradphi(:,1)**2 + 
     &                          gradphi(:,2)**2 + 
     &                          gradphi(:,3)**2 )**0.5
              maxgradphi = maxval(gradphimag)
	      write(*,1000) maxgradphi, myrank
              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	      call write_gradphi(myrank, lstep, nshg, 3, 
     &                           x, y, gradphi, gradphimag)
 1000 format ("Maximum LS gradient = ",f12.6," on proc",i6)
	    endif
c
c...Reset BC codes of primary vertices
c (need to remove dirchlet bc on primary vertices
c  by removing the 8th bit (val=128)
c
            if((iLSet.eq.2).and.(ilss.eq.1).and.(iSolvLSSclr2.eq.2))
     &      then
               where (iBCredist.eq.1)
                 iBC(:) = iBC(:) - 128   ! remove prescription on scalar 2
               endwhere
            endif

c
c... Write out Redistancing primary vertice information
c
            if ((i_primvert.eq.1) .and. (mod(lstep, ntout).eq.0)) then
                primvertval(:,2) = primvert(:) * 1.0
                call write_primvert(myrank, lstep, nshg, 2,
     &                              primvertval)


            endif

!MR CHANGE
            if ( icomputevort == 1) then
 
              ! vorticity components and magnitude
              vorticity(:,1) = GradV(:,8)-GradV(:,6) !omega_x
              vorticity(:,2) = GradV(:,3)-GradV(:,7) !omega_y
              vorticity(:,3) = GradV(:,4)-GradV(:,2) !omega_z
              vorticity(:,4) = sqrt(   vorticity(:,1)*vorticity(:,1)
     &                               + vorticity(:,2)*vorticity(:,2)
     &                               + vorticity(:,3)*vorticity(:,3) )
              ! Q
              strain(:,1) = GradV(:,1)                  !S11
              strain(:,2) = 0.5*(GradV(:,2)+GradV(:,4)) !S12
              strain(:,3) = 0.5*(GradV(:,3)+GradV(:,7)) !S13
              strain(:,4) = GradV(:,5)                  !S22
              strain(:,5) = 0.5*(GradV(:,6)+GradV(:,8)) !S23
              strain(:,6) = GradV(:,9)                  !S33
 
              vorticity(:,5) = 0.25*( vorticity(:,4)*vorticity(:,4)  !Q
     &                            - 2.0*(      strain(:,1)*strain(:,1)
     &                                    + 2* strain(:,2)*strain(:,2)
     &                                    + 2* strain(:,3)*strain(:,3)
     &                                    +    strain(:,4)*strain(:,4)
     &                                    + 2* strain(:,5)*strain(:,5)
     &                                    +    strain(:,6)*strain(:,6)))

            endif
!MR CHANGE END

c
c .. write out the solution
c
            if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
                 IF(iBT.eq.1) call banmaCorrector(banma, ibreakupFlag, yold)
                 call restar ('out ',  yold  ,ac)
                 IF(iBT.eq.1) yold(:,7) = yold(:,6)
c               if(iofieldv.ne.0) then
c                 call pstdrv (y,         ac,         x,
c     &                        iBC,       BC,
c     &                        iper,      ilwork,     shp,
c     &                        shgl,      shpb,       shglb,
c     &                        ifath,     velbar,     nsons )
c
c               endif
               if(ideformwall.eq.1) 
     &              call write_displ(myrank, lstep, nshg, 3, uold ) 
               if (ioform .eq. 2) 
     &              call stsWriteStats()
               if(ivort == 1) then 
                 call write_field(myrank,'a','vorticity',9,vorticity,
     &                       'd',nshg,5,lstep)
               endif


            endif
c SUBROUTINE
c ... update the flow rate history for the impedance convolution
c    
            if(numImpSrfs.gt.zero) then
                call GetFlowQ(NewQImp,y,nsrflistImp,numImpSrfs) !flow Q for imp BC
                do j=1, ntimeptpT
                    QHistImp(j,:)=QHistImp(j+1,:)
                enddo
                QHistImp(ntimeptpT+1,1:numImpSrfs) = 
     &          NewQImp(1:numImpSrfs)

c
c.... write out the new history of flow rates to Qhistor.dat
c      
                if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and.
     &               (myrank .eq. zero)) then
                   open(unit=816, file='Qhistor.dat',status='replace')
                   write(816,*) ntimeptpT
                   do j=1,ntimeptpT+1
                      write(816,*) (QHistImp(j,n),n=1, numImpSrfs)
                   enddo
                   close(816)
                endif
             endif
c
c.... compute the consistent boundary flux
c
            if(abs(itwmod).ne.1 .and. iowflux.eq.1) then
               call Bflux ( yold,      acold,      uold,    x,
     &                      shp,       shgl,       shpb,   
     &                      shglb,     ilwork,     iBC,
     &                      BC,        iper,       wallssVec)
            endif

c...  dump TIME SERIES

            if (exts) then
               if (mod(lstep-1,freq).eq.0) then

                  if (numpe > 1) then
                     do jj = 1, ntspts
             vartssoln((jj-1)*numvar+1:jj*numvar)=varts(jj,1:numvar)
                        ivarts=zero
                     enddo

                     do k=1,numvar*ntspts
                        if(vartssoln(k).ne.zero) ivarts(k)=1
                     enddo
           call MPI_REDUCE(vartssoln, vartssolng, numvar*ntspts,
     &                    MPI_DOUBLE_PRECISION, MPI_SUM, master,
     &                    MPI_COMM_WORLD, ierr)

              call MPI_REDUCE(ivarts, ivartsg, numvar*ntspts,
     &                    MPI_INTEGER, MPI_SUM, master,
     &                    MPI_COMM_WORLD, ierr)

                     if (myrank.eq.zero) then
                        do jj = 1, ntspts

                           indxvarts = (jj-1)*numvar
                           do k=1, numvar
                              if(ivartsg(indxvarts+k).ne.0) then ! none of the vartssoln(parts) were non zero
                                 varts(jj,k)=vartssolng(indxvarts+k)/
     &                                ivartsg(indxvarts+k)
                              endif
                           enddo ! do k
                        enddo ! do jj
                     endif !only on master
                  endif !only if numpe > 1

c Do not search anymore:
	tssearch = 0

                  if (myrank.eq.zero) then
                     do jj = 1, ntspts
			if (numvar.ge.15) then
			   if (varts(jj,15).le.0) then
			     iphase = 1
			   else
			     iphase = 0
			   end if
			else
			 iphase = 0
			end if
			varts(jj,5) = delt(itseq)*real(freq)   ! delta t has to be recorded; multiply it by the varts stepping
			recnum = (istep-1)/freq*ntspts+jj
                        if (numvar.eq.19) then
			varts(jj,19) = one  ! Put the Heps value here (depends upon varts(jj,15) = distance field value)
             if  (abs(varts(jj,15)) .le. epsilon_lsd_tmp) then
          varts(jj,19) = 0.5*(one + varts(jj,15) / epsilon_lsd_tmp +
     &                    (sin(pi*varts(jj,15)/epsilon_lsd_tmp))/pi)
              elseif (varts(jj,15) .lt. - epsilon_lsd_tmp) then
                 varts(jj,19) = zero
              endif
                         write(10101,'(2I8, I3, 19E15.7)',REC=recnum)
     &            lstep-1,jj,iphase,(varts(jj,k), k=1, 19)   !  including the distance field (k = 15), gradient (16..18) and Heps value (19)
			else if (numvar.eq.15) then
			 write(10101,'(2I8, I3, 15E15.7)',REC=recnum) 
     &            lstep-1,jj,iphase,(varts(jj,k), k=1, 15)   !  including the distance field (k = 15)
			else
c			if (jj.eq.1) write(*,*) 'numvar = ', numvar, '; point is there'
                         write(10101,'(2I8, I3, 15E15.7)',REC=recnum)
     &            lstep-1,jj,iphase,(varts(jj,k), k=1, 14), 0.0E0  
			end if
                        ifile = 1000+jj
                        if (ntspts.gt.80) then                  ! no more than 96 files can be opened simulatneously
                           fvarts='varts/varts'
                           fvarts=trim(fvarts)//trim(cname2(jj))
                           fvarts=trim(fvarts)//trim(cname2(lskeep))
                           fvarts=trim(fvarts)//'.dat'
                           fvarts=trim(fvarts)
c                           open(unit=ifile, file=fvarts,
c     &                          position='append')
                        end if
c            write(ifile,555) lstep-1, (varts(jj,k), k=1, numvar) 
c                        call flush(ifile)
                        if (ntspts.gt.80) then
c                           close(ifile)
                        end if

                       if (((irs .ge. 1) .and. (ntspts.le.80) .and.
     &                       (mod(lstep, ntout) .eq. 0))) then
c                           close(ifile)
                           fvarts='varts/varts'
                           fvarts=trim(fvarts)//trim(cname2(jj))
                           fvarts=trim(fvarts)//trim(cname2(lskeep))
                           fvarts=trim(fvarts)//'.dat'
                           fvarts=trim(fvarts)
c                           open(unit=ifile, file=fvarts,
c     &                          position='append')
                        endif !only when dumping restart
                     enddo
                  endif !only on master

                  varts(:,:) = zero ! reset the array for next step

 555              format(i6,20(2x,E12.5e2))

               endif
            endif

c
c.... update and the aerodynamic forces
c
!	if (myrank.eq.master) write(*,*) 'calling forces, l. 974 itrdrv'
            call forces ( yold,  ilwork )
            
            if((irscale.ge.0).or.(itwmod.gt.0).or.(iDNS.ne.0)) 
     &           call getvel (yold,     ilwork, iBC,
     &                        nsons,    ifath, velbar)

            if((irscale.ge.0).and.(myrank.eq.master)) then
               call genscale(yold,       x,    iBC)   ! Corrected by Igor A. Bolotnov, 07/06/2009
c               call genscale(yold,       x,       iper, 
c     &                       iBC,     ifath,   velbar,
c     &                       nsons)
            endif
c
c....  print out results.
c
            ntoutv=max(ntout,100)   ! velb is not needed so often
            if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
               if( (mod(lstep, ntoutv) .eq. 0) .and.
     &              ((irscale.ge.0).or.(itwmod.gt.0) .or. 
     &              ((nsonmax.eq.1).and.(iLES.gt.0))))
     &              call rwvelb  ('out ',  velbar  ,ifail)
            endif
c
c.... end of the NSTEP and NTSEQ loops
c
c
c.... -------------------> error calculation  <-----------------
c 
            if(ierrcalc.eq.1 .or. ioybar.eq.1) then
c$$$c
c$$$c compute average
c$$$c
c$$$               tfact=one/istep
c$$$               ybar =tfact*yold + (one-tfact)*ybar

c compute average
c ybar(:,1) - ybar(:,3) is average velocity components
c ybar(:,4) is average pressure
c ybar(:,5) is average speed
c averaging procedure justified only for identical time step sizes
c istep is number of time step
c
               tfact=one/istep

c ybar to contain the averaged ((u,v,w),p)-field
c and speed average, i.e., sqrt(u^2+v^2+w^2)

               ybar(:,1) = tfact*yold(:,1) + (one-tfact)*ybar(:,1)
               ybar(:,2) = tfact*yold(:,2) + (one-tfact)*ybar(:,2)
               ybar(:,3) = tfact*yold(:,3) + (one-tfact)*ybar(:,3)
               ybar(:,4) = tfact*yold(:,4) + (one-tfact)*ybar(:,4)
               ybar(:,5) = tfact*sqrt(yold(:,1)**2+yold(:,2)**2+
     &                     yold(:,3)**2) + (one-tfact)*ybar(:,5)

!MR CHANGD
               if(ivort == 1) then 
                 vorticitybar(:,1) = tfact*vorticity(:,1) + 
     &                           (one-tfact)*vorticitybar(:,1)
                 vorticitybar(:,2) = tfact*vorticity(:,2) + 
     &                           (one-tfact)*vorticitybar(:,2)
                 vorticitybar(:,3) = tfact*vorticity(:,3) + 
     &                           (one-tfact)*vorticitybar(:,3)
                 vorticitybar(:,4) = tfact*vorticity(:,4) + 
     &                           (one-tfact)*vorticitybar(:,4)
               endif

               if(abs(itwmod).ne.1 .and. iowflux.eq.1) then 
                 wallssVecBar(:,1) = tfact*wallssVec(:,1)
     &                               +(one-tfact)*wallssVecBar(:,1)
                 wallssVecBar(:,2) = tfact*wallssVec(:,2)
     &                               +(one-tfact)*wallssVecBar(:,2)
                 wallssVecBar(:,3) = tfact*wallssVec(:,3)
     &                               +(one-tfact)*wallssVecBar(:,3)
               endif
!MR CHANGE END

c
c compute rms
c
               rerr(:, 7)=rerr(:, 7)+(yold(:,1)-ybar(:,1))**2
               rerr(:, 8)=rerr(:, 8)+(yold(:,2)-ybar(:,2))**2
               rerr(:, 9)=rerr(:, 9)+(yold(:,3)-ybar(:,3))**2
               rerr(:,10)=rerr(:,10)+(yold(:,4)-ybar(:,4))**2

            endif
            
            if(istop.eq.1000) exit ! stop when delta small (see rstatic)
! 2000    continue
        istp = istp + 1
        enddo
 2001    continue

CAD         tcorecp2 = second(0)
CAD         tcorewc2 = second(-1)
         
CAD         write(6,*) 'T(core) cpu-wallclock = ',tcorecp2-tcorecp1,
CAD     &                                        tcorewc2-tcorewc1

         if (numpe > 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         if(myrank.eq.0)  then
            tcorecp2 = TMRC()
            write(6,*) 'T(core) cpu = ',tcorecp2-tcorecp1
            write(6,*) '(Elm. form.',tcorecp(1),',Lin. alg. sol.',
     &                    tcorecp(2),')'
         endif

         call print_system_stats(tcorecp)
         call print_mesh_stats()
         call print_mpi_stats()
         if (numpe > 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)

 3000 continue
c      open(unit=87,file="spmass.dat",status="unknown")
c      write(87,197)(splag(j),j=1,10)
c      close(87)
 197  format(10(2x,e14.7))
c
c.... ---------------------->  Post Processing  <----------------------

	if (myrank.eq.master.and.exts)
     1   close(10101)		! Close the binary output file

c
c.... print out the last step
c
      if ((irs .ge. 1) .and. ((mod(lstep, ntout) .ne. 0) .or.
     &     (nstp .eq. 0) .or. (stopjob.ne.0) )) then
         if(
     &              ((irscale.ge.0).or.(itwmod.gt.0) .or. 
     &              ((nsonmax.eq.1).and.iLES.gt.0)))
     &              call rwvelb  ('out ',  velbar  ,ifail)
         IF(iBT.eq.1) call banmaCorrector(banma, ibreakupFlag, yold)
         call restar ('out ',  yold  ,ac)
         IF(iBT.eq.1) yold(:,7) = yold(:,6)

c         if(iofieldv.ne.0)
c     &           call pstdrv (yold,      ac,         x,
c     &                        iBC,       BC,
c     &                        iper,      ilwork,     shp,
c     &                        shgl,      shpb,       shglb,
c     &                        ifath,     velbar,     nsons )
         if (i_gradphi.eq.1) then
	    idflx = 0 
            if(idiff >= 1 )  idflx= (nflow-1) * nsd
            if (isurf == 1) idflx=nflow*nsd
	    call getgradphi(x, y, shp,       shgl,
     &                      shpb,          shglb, gradphi)
            gradphimag(:) = ( gradphi(:,1)**2 + 
     &                        gradphi(:,2)**2 + 
     &                        gradphi(:,3)**2 )**0.5
            maxgradphi = maxval(gradphimag)
	    write(*,1000) maxgradphi, myrank
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	    call write_gradphi(myrank, lstep, nshg, 3, 
     &                           x, y, gradphi, gradphimag)
	 endif
c
c... Write out Redistancing primary vertice information
c
         if (i_primvert.eq.1) then
             primvertval(:,2) = primvert(:) * 1.0
             call write_primvert(myrank, lstep, nshg, 2,
     &                           primvertval)
         endif
         if(ideformwall.eq.1) 
     &        call write_displ(myrank, lstep, nshg, 3, u ) 
      endif


c
c Free memory for freezing LS Scalar 2
c
            if (iSolvLSSclr2.eq.2) then
              deallocate (iBCredist)
              deallocate (BCredist)
            endif

         lesId   = numeqns(1)
         call saveLesRestart( lesId,  aperm , nshg, myrank, lstep,
     &                        nPermDims )


      if(ierrcalc.eq.1) then
c
c.....smooth the error indicators
c
        do i=1,ierrsmooth
            call errsmooth( rerr, x, iper, ilwork, shp, shgl, iBC )
        end do
c
c.... open the output file
c
         call write_field(myrank, 'a', 'errors', 6,
     &                    rerr, 'd', nshg,numerr,lstep)
      endif

      if(ioybar.eq.1) then

         call write_field(myrank,'a','ybar',4,
     &                    ybar,'d',nshg,ndof,lstep)

!MR CHANGE
         if(ivort == 1) then
           call write_field(myrank,'a','vorticitybar',12,
     &                    vorticitybar,'d',nshg,4,lstep)
           deallocate(vorticitybar)
         endif
!MR CHANGE END



c         write (fmt2,"('(''restart.'',i',i1,',1x)')") itmp
c         write (fname2,fmt2) lstep
c
c         fname2 = trim(fname2) // cname(myrank+1)
c
c.... open  files
c
c         call openfile( trim(fname2)//char(0), 'append?'//char(0), 
c     &                   irstin )
c
c         fnamer2 = 'ybar'
c         isize = nshg*5
c         nitems = 3
c         iarray(1) = nshg
c         iarray(2) = 5
c         iarray(3) = lstep
c         call writeheader(irstin, trim(fnamer2)//char(0),iarray,
c      &                    nitems, isize, 'double'//char(0), iotype )
c
c         nitems = nshg*5
c         call writedatablock(irstin, trim(fnamer2)//char(0),ybar, 
c      &         nitems, 'double'//char(0), iotype)
c
c         call closefile( irstin, 'append'//char(0) )

      endif

!MR CHANGE
      if(ivort == 1) then
         deallocate(strain,vorticity)
      endif
!MR CHANGE END


      if ( ( ihessian .eq. 1 ) .and. ( numpe < 2 )  )then
          uhess = zero
          gradu = zero
          tf = zero

          do ku=1,nshg
c           tf(ku,1) = x(ku,1)**2+2*x(ku,1)*x(ku,2)
            tf(ku,1) = x(ku,1)**3
          end do

          call hessian( yold, x,     shp,  shgl,   iBC, 
     &                  shpb, shglb, iper, ilwork, uhess, gradu )

          call write_hessian( uhess, gradu, nshg )
      endif

!... save the d2wal information into the last restart files
      if (iRANS .lt. 0 .or. (iDNS .ne. 0 .and. abs(itwmod) .eq. 1)
     &                 .or. id2w .eq. 1) then
        call write_field(myrank, 'a', 'dwal', 4, d2wal, 'd', nshg,
     &                       4, lstep)
        deallocate(d2wal,d2wall,x2wall,y2wall,z2wall)
        if(myrank.eq.master) write(*,*)'Save dwal to restarts!'
      endif

c
c write phasta status file
c
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (myrank .eq. master) then
         open (unit=istat,  file=fstat,  status='unknown')
         write(istat,*) "done"
      endif
      close(istat)
c
c.... close history and aerodynamic forces files
c
      if (myrank .eq. master) then
         close (ihist)
         close (iforce)
         if (iLSet .eq. 2) then
           close (ivhist)
         endif
         if(exts) then
            do jj=1,ntspts
               close(1000+jj)
            enddo
         endif
      endif
      do isrf = 0,MAXSURF
         if ( nsrflist(isrf).ne.0 ) then
            iunit=60+isrf
            close(iunit)
         endif
      enddo
 5    format(1X,F15.10,3X,F15.10,3X,F15.10,3X,F15.10)
 444  format(6(2x,e14.7))
c
c.... end
c
      if(nsolflow.eq.1) then
         deallocate (lhsK)
         deallocate (lhsP)
         deallocate (aperm)
         deallocate (atemp)
      endif
      if(nsclrsol.gt.0) then
         deallocate (lhsS)
         deallocate (apermS)
         deallocate (atempS)
      endif
      
      if(iabc==1) deallocate(acs)

      return
      end
      
      subroutine lesmodels(y,     ac,        shgl,      shp, 
     &                     iper,  ilwork,    rowp,      colm,    
     &                     nsons, ifath,     x,   
     &                     iBC,   BC)
      
      include "common.h"

      real*8    y(nshg,ndof),              ac(nshg,ndof),           
     &            x(numnp,nsd),
     &            BC(nshg,ndofBC)
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT)

c
      integer   rowp(nshg,nnz),         colm(nshg+1),
     &            iBC(nshg),
     &            ilwork(nlwork),
     &            iper(nshg)
      dimension ifath(numnp),    nsons(nfath)

      real*8, allocatable, dimension(:) :: fwr2,fwr3,fwr4
      real*8, allocatable, dimension(:) :: stabdis,cdelsq1
      real*8, allocatable, dimension(:,:) :: xavegt, xavegt2,xavegt3

      if( (iLES.gt.1) )   then ! Allocate Stuff for advanced LES models
         allocate (fwr2(nshg))
         allocate (fwr3(nshg))
         allocate (fwr4(nshg))
         allocate (xavegt(nfath,12))
         allocate (xavegt2(nfath,12))
         allocate (xavegt3(nfath,12))
         allocate (stabdis(nfath))
      endif

c.... get dynamic model coefficient
c
      ilesmod=iLES/10  
c
c digit bit set filter rule, 10 bit set model
c
      if (ilesmod.eq.0) then    ! 0 < iLES< 10 => dyn. model calculated
                                ! at nodes based on discrete filtering


         if(isubmod.eq.2) then
            call SUPGdis(y,      ac,        shgl,      shp, 
     &                   iper,   ilwork,    
     &                   nsons,  ifath,     x,   
     &                   iBC,    BC, stabdis, xavegt3)
         endif

         if( ((isubmod.eq.0).or.(isubmod.eq.2)))then ! If no
                                                     ! sub-model
                                                     ! or SUPG
                                                     ! model wanted

            if(i2filt.eq.0)then ! If simple filter
              
               if(modlstats .eq. 0) then ! If no model stats wanted
                  call getdmc (y,       shgl,      shp, 
     &                         iper,       ilwork,    nsons,
     &                         ifath,      x)
               else             ! else get model stats 
                  call stdfdmc (y,       shgl,      shp, 
     &                          iper,       ilwork,    nsons,
     &                          ifath,      x)
               endif            ! end of stats if statement  

            else                ! else if twice filtering

               call widefdmc(y,       shgl,      shp, 
     &                       iper,       ilwork,    nsons,
     &                       ifath,      x)

               
            endif               ! end of simple filter if statement

         endif                  ! end of SUPG or no sub-model if statement


         if( (isubmod.eq.1) ) then ! If DFWR sub-model wanted
            call cdelBHsq (y,       shgl,      shp, 
     &                     iper,       ilwork,    nsons,
     &                     ifath,      x,         cdelsq1)
            call FiltRat (y,       shgl,      shp, 
     &                    iper,       ilwork,    nsons,
     &                    ifath,      x,         cdelsq1,
     &                    fwr4,       fwr3)

            
            if (i2filt.eq.0) then ! If simple filter wanted
               call DFWRsfdmc(y,       shgl,      shp, 
     &                        iper,       ilwork,    nsons,
     &                        ifath,      x,         fwr2, fwr3) 
            else                ! else if twice filtering wanted 
               call DFWRwfdmc(y,       shgl,      shp, 
     &                        iper,       ilwork,    nsons,
     &                        ifath,      x,         fwr4, fwr4) 
            endif               ! end of simple filter if statement
             
         endif                  ! end of DFWR sub-model if statement

         if( (isubmod.eq.2) )then ! If SUPG sub-model wanted
            call dmcSUPG (y,           ac,         shgl,      
     &                    shp,         iper,       ilwork,    
     &                    nsons,       ifath,      x,
     &                    iBC,    BC,  rowp,       colm,
     &                    xavegt2,    stabdis)
         endif

         if(idis.eq.1)then      ! If SUPG/Model dissipation wanted
            call ediss (y,        ac,      shgl,      
     &                  shp,      iper,       ilwork,    
     &                  nsons,    ifath,      x,
     &                  iBC,      BC,  xavegt)
         endif

      endif                     ! end of ilesmod
      
      if (ilesmod .eq. 1) then  ! 10 < iLES < 20 => dynamic-mixed
                                ! at nodes based on discrete filtering
         call bardmc (y,       shgl,      shp, 
     &                iper,    ilwork,    
     &                nsons,   ifath,     x) 
      endif
      
      if (ilesmod .eq. 2) then  ! 20 < iLES < 30 => dynamic at quad
                                ! pts based on lumped projection filt. 

         if(isubmod.eq.0)then
            call projdmc (y,       shgl,      shp, 
     &                    iper,       ilwork,    x) 
         else
            call cpjdmcnoi (y,      shgl,      shp, 
     &                      iper,   ilwork,       x,
     &                      rowp,   colm, 
     &                      iBC,    BC)
         endif

      endif

      if( (iLES.gt.1) )   then ! Deallocate Stuff for advanced LES models
         deallocate (fwr2)
         deallocate (fwr3)
         deallocate (fwr4)
         deallocate (xavegt)
         deallocate (xavegt2)
         deallocate (xavegt3)
         deallocate (stabdis)
      endif
      return
      end
