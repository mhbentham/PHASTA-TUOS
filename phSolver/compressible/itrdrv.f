              subroutine itrdrv (y,         ac,         x,         
     &                   iBC,       BC,         
     &                   iper,      ilwork,     shp,       
     &                   shgl,      shpb,       shglb,
     &                   ifath,     velbar,     nsons ) 
c
c----------------------------------------------------------------------
c
c This iterative driver is the semi-discrete, predictor multi-corrector 
c algorithm. It contains the Hulbert Generalized Alpha method which
c is 2nd order accurate for Rho_inf from 0 to 1.  The method can be
c made  first-order accurate by setting Rho_inf=-1. It uses a
c GMRES iterative solver.
c
c working arrays:
c  y      (nshg,ndof)           : Y variables
c  x      (nshg,nsd)            : node coordinates
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  iper   (nshg)                : periodicity table
c
c shape functions:
c  shp    (nshape,ngauss)        : interior element shape functions
c  shgl   (nsd,nshape,ngauss)    : local shape function gradients
c  shpb   (nshapeb,ngaussb)      : boundary element shape functions
c  shglb  (nsd,nshapeb,ngaussb)  : bdry. elt. shape gradients
c
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use pvsQbi     !gives us splag (the spmass at the end of this run 
      use specialBC  !gives us itvn
      use timedata   !allows collection of time series

        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
      
c
        dimension y(nshg,ndof),            ac(nshg,ndof),           
     &		  yold(nshg,ndof),         acold(nshg,ndof),           
     &            x(numnp,nsd),            iBC(nshg),
     &            BC(nshg,ndofBC),         ilwork(nlwork),
     &            iper(nshg)
c
        dimension res(nshg,nflow),         BDiag(nshg,nflow,nflow),
     &            rest(nshg),              solinc(nshg,ndof)
c     
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
        real*8   almit, alfit, gamit, soln, asoln, asolng
        dimension ifath(numnp),    velbar(nfath,ndof),  nsons(nfath)
        real*8 rerr(nshg,10),ybar(nshg,nflow-1)
        character*20    fname1,fmt1
        character*5  cname
        integer ifuncs(6)
c
c  Here are the data structures for sparse matrix GMRES
c
       integer, allocatable, dimension(:,:) :: rowp
       integer, allocatable, dimension(:) :: colm
       real*8, allocatable, dimension(:,:) :: lhsK
       real*8, allocatable, dimension(:,:) :: EGmass
       real*8, allocatable, dimension(:,:) :: EGmasst

       inquire(file='xyzts.dat',exist=exts)
       
       if(exts) then
         
          open(unit=626,file='xyzts.dat',status='old')
          read(626,*) ntspts, freq, tolpt, iterat, varcod
          call sTD              ! sets data structures
          
          do jj=1,ntspts        ! read coordinate data where solution desired
             read(626,*) ptts(jj,1),ptts(jj,2),ptts(jj,3)
          enddo
          
          varts = zero
          
       endif
 
c
c
c.... open history and aerodynamic forces files
c
        if (myrank .eq. master) then
          open (unit=ihist,  file=fhist,  status='unknown')
          open (unit=iforce, file=fforce, status='unknown')
        endif
c
c Ridiculous hack to get stable outflow
c
c  to enable this you must ask for a sponge AND then set its growth rate
c        to zero.  Hopefully this cryptic combination does not cause
c        other users grief (until we do things a bit better)

       
clast-working G6        if(matflg(5,1).eq.4.and.grthOSponge.eq.zero) then
        if(matflg(5,1).eq.4.and.grthISponge.eq.zero) then
           write(*,*) "WARNING: modification of outflow BC at"
        do i =1,nshg
           if(x(i,2).gt.radsponge.and.btest(ibc(i),2)) then
cONLY IF IC GOOD           if(btest(ibc(i),0)) then
! this is a pressure outflow and above the cutoff for the edge of the
 ! shear layer
              if(btest(ibc(i),1)) then !temperature already set
              else
                 write(*,407) i,x(i,1),x(i,2),x(i,3), ibc(i)
                 ibc(i)=ibc(i)+2 ! sets temperature
                 BC(i,2)=.00688
cONLY IF IC GOOD                  BC(i,2)=y(i,5)
              endif
           endif
        enddo
 407    format(i6,3(2x,e14.7),2x,i4)
        endif

c
c.... initialize
c
        ifuncs  = 0                      ! func. evaluation counter
        istep  = 0
        ntotGM = 0                      ! number of GMRES iterations
        time   = 0
        yold   = y
        acold  = ac
        if (mod(impl(1),100)/10 .eq. 1) then
c
c     generate the sparse data fill vectors
c
           allocate  (rowp(nshg,nnz))
           allocate  (colm(nshg+1))
           call genadj(colm, rowp, icnt ) ! preprocess the adjacency list

           nnz_tot=icnt         ! this is exactly the number of non-zero 
                                ! blocks on this proc
           allocate (lhsK(nflow*nflow,nnz_tot))
        endif
        if (mod(impl(1),100)/10 .eq. 3) then
c
c     generate the ebe data fill vectors
c
           nedof=nflow*nshape
           allocate  (EGmass(numel,nedof*nedof))
        endif
cc
        rerr = zero
        ybar = y
c
c.... loop through the time sequences
c
        do 3000 itsq = 1, ntseq
        itseq = itsq
c
c.... set up the current parameters
c
        nstp   = nstep(itseq)
        nitr   = niter(itseq)
        LCtime = loctim(itseq)
c
        call itrSetup ( y,  acold)
        isclr=0

        niter(itseq)=0          ! count number of flow solves in a step
                                !  (# of iterations)
        do i=1,seqsize
           if(stepseq(i).eq.0) niter(itseq)=niter(itseq)+1
        enddo
        nitr = niter(itseq)
c
c.... determine how many scalar equations we are going to need to solve
c
        nsclrsol=nsclr          ! total number of scalars solved. At
                                ! some point we probably want to create
                                ! a map, considering stepseq(), to find
                                ! what is actually solved and only
                                ! dimension EGmasst to the appropriate
                                ! size.

        if(nsclrsol.gt.0)allocate  (EGmasst(numel*nshape*nshape
     &                              ,nsclrsol))
 
c
c.... loop through the time steps
c
        ttim(1) = REAL(secs(0.0)) / 100.
        ttim(2) = secs(0.0)

        tcorecp1 = REAL(secs(0.0)) / 100.
        tcorewc1 = secs(0.0)

        rmub=datmat(1,2,1)
        if(rmutarget.gt.0) then
           rmue=rmutarget
           xmulfact=(rmue/rmub)**(1.0/nstp)
           if(myrank.eq.0) then
              write(*,*) 'viscosity will by multiplied by ', xmulfact
              write(*,*) 'to bring it from ', rmub,' down to ', rmue
           endif
           datmat(1,2,1)=datmat(1,2,1)/xmulfact ! make first step right
        else
           rmue=datmat(1,2,1)   ! keep constant
           xmulfact=one
        endif
c
c keep sponge from conflicting with BC's
c
c$$$        if(matflg(5,1).eq.4) then
c$$$        do i=1,numnp
c$$$           r=sqrt(x(i,1)**2+x(i,2)**2)
c$$$           if(r.gt.radst.and.btest(iBC(i),1))  iBC(i)=iBC(i)-2
c$$$           if(r.gt.radst.and.btest(iBC(i),2))  iBC(i)=iBC(i)-4
c$$$           if(z.gt.zstart.and.btest(iBC(i),1)) iBC(i)=iBC(i)-2
c$$$           if(r.gt.zstart.and.btest(iBC(i),2)) iBC(i)=iBC(i)-4
c$$$        enddo
c$$$        endif

        do 2000 istp = 1, nstp

c           call rerun_check(stopjob)
cc          if(stopjob.ne.0) goto 2001
c
c Decay of scalars
c
           if(nsclr.gt.0 .and. tdecay.ne.1) then
              yold(:,6:ndof)=y(:,6:ndof)*tdecay
              BC(:,7:6+nsclr)= BC(:,7:6+nsclr)*tdecay
           endif

           if(nosource.eq.1) BC(:,7:6+nsclr)= BC(:,7:6+nsclr)*0.8


c           xi=istp*one/nstp
c           datmat(1,2,1)=rmub*(1.0-xi)+xi*rmue
           datmat(1,2,1)=xmulfact*datmat(1,2,1)


            if(iLES.gt.0) then
c
c.... get dynamic model coefficient
c
            ilesmod=iLES/10  
c
c digit bit set filter rule, 10 bit set model
c
            if (ilesmod.eq.0) then ! 0 < iLES < 10 => dyn. model calculated
                                   ! at nodes based on discrete filtering
               call getdmc (yold,       shgl,      shp, 
     &                      iper,       ilwork,    nsons,
     &                      ifath,      x)
            endif
            if (ilesmod .eq. 1) then ! 10 < iLES < 20 => dynamic-mixed
                                     ! at nodes based on discrete filtering
               call bardmc (yold,       shgl,      shp, 
     &                      iper,       ilwork,    
     &                      nsons,      ifath,     x) 
            endif
            if (ilesmod .eq. 2) then ! 20 < iLES < 30 => dynamic at quad
                                     ! pts based on lumped projection filt. 
               call projdmc (yold,       shgl,      shp, 
     &                       iper,       ilwork,    x) 
            endif
c
            endif


c
c.... set traction BCs for modeled walls
c
            if (itwmod.ne.0) then   !wallfn check
               call asbwmod(yold,   acold,   x,      BC,     iBC,
     &                      iper,   ilwork,  ifath,  velbar)
            endif
c
c.... -----------------------> predictor phase <-----------------------
c
            call itrPredict(   yold,    acold,    y,   ac )
            call itrBC (y,  ac,  iBC,  BC,  iper, ilwork)
            isclr = zero
            if (nsclr.gt.zero) then
            do isclr=1,nsclr
               call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
            enddo
            endif
c
c.... --------------------> multi-corrector phase <--------------------
c
           iter=0
            ilss=0  ! this is a switch thrown on first solve of LS redistance
            do istepc=1,seqsize
               icode=stepseq(istepc)
               if(mod(icode,10).eq.0) then ! this is a solve
                  isolve=icode/10
                  if(isolve.eq.0) then   ! flow solve (encoded as 0)
c
                     etol=epstol(1)
                     iter   = iter+1
                     ifuncs(1)  = ifuncs(1) + 1
c     
c.... reset the aerodynamic forces
c     
                     Force(1) = zero
                     Force(2) = zero
                     Force(3) = zero
                     HFlux    = zero
c     
c.... form the element data and solve the matrix problem
c     
c.... explicit solver
c     
                     if (impl(itseq) .eq. 0) then
                        if (myrank .eq. master)
     &                       call error('itrdrv  ','impl ',impl(itseq))
                     endif
                     if (mod(impl(1),100)/10 .eq. 1) then  ! sparse solve
c     
c.... preconditioned sparse matrix GMRES solver
c     
                        lhs = 1 - min(1,mod(ifuncs(1)-1,LHSupd(1))) 
                        iprec=lhs
                        nedof = nflow*nshape
c                        write(*,*) 'lhs=',lhs
                      call SolGMRs (y,             ac,            yold,
     &                       acold,         x,
     &                       iBC,           BC,
     &                       colm,          rowp,          lhsk,
     &                       res,
     &                       BDiag,         a(mHBrg),      a(meBrg),
     &                       a(myBrg),      a(mRcos),      a(mRsin),
     &                       iper,          ilwork,
     &                       shp,           shgl,
     &                       shpb,          shglb,         solinc,
     &                       rerr)
                      else if (mod(impl(1),100)/10 .eq. 2) then ! mfg solve
c     
c.... preconditioned matrix-free GMRES solver
c     
                        lhs=0
                        iprec = 1 - min(1,mod(ifuncs(1)-1,LHSupd(1))) 
                        nedof = 0
                        call SolMFG (y,             ac,            yold,
     &                       acold,         x,
     &                       iBC,           BC,
     &                       res,           
     &                       BDiag,         a(mHBrg),      a(meBrg),
     &                       a(myBrg),      a(mRcos),      a(mRsin),
     &                       iper,          ilwork,
     &                       shp,           shgl,
     &                       shpb,          shglb,         solinc, 
     &                       rerr)
c     
                     else if (mod(impl(1),100)/10 .eq. 3) then ! ebe solve
c.... preconditioned ebe matrix GMRES solver
c     

                        lhs = 1 - min(1,mod(ifuncs(1)-1,LHSupd(1))) 
                        iprec = lhs
                        nedof = nflow*nshape
c                        write(*,*) 'lhs=',lhs
                      call SolGMRe (y,             ac,            yold,
     &                       acold,         x,
     &                       iBC,           BC,
     &                       EGmass,        res,
     &                       BDiag,         a(mHBrg),      a(meBrg),
     &                       a(myBrg),      a(mRcos),      a(mRsin),
     &                       iper,          ilwork,
     &                       shp,           shgl,
     &                       shpb,          shglb,         solinc,
     &                       rerr)
                     endif
c     
                else          ! solve a scalar  (encoded at isclr*10)
                     ifuncs(isclr+2)  = ifuncs(isclr+2) + 1
                     etol=epstol(isclr+1)
                     isclr=isolve
                     if((iLSet.eq.2).and.(ilss.eq.0)
     &                    .and.(isclr.eq.2)) then 
                        ilss=1  ! throw switch (once per step)
c     
c... copy the first scalar at t=n+1 into the second scalar of the 
c... level sets
c     
                     y(:,6)    = yold(:,6)  + (y(:,6)-yold(:,6))/alfi
                     y(:,7)    =  y(:,6)
                     yold(:,7) = y(:,7)
                     ac(:,7)   = zero
c     
                     call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
c     
c....store the flow alpha, gamma parameter values and assigm them the 
c....Backward Euler parameters to solve the second levelset scalar
c     
                        alfit=alfi
                        gamit=gami
                        almit=almi
                        alfi = 1
                        gami = 1
                        almi = 1
                     endif
c     
                     lhs = 1 - min(1,mod(ifuncs(isclr+2)-1,
     &                                       LHSupd(isclr+2)))
                     iprec = lhs
                     call SolGMRSclr(y,             ac,         yold,
     &		          acold,         EGmasst(1,isclr),
     &                    x,
     &                    iBC,           BC,          
     &                    rest,           
     &                    a(mHBrg),      a(meBrg),
     &                    a(myBrg),      a(mRcos),    a(mRsin),
     &                    iper,          ilwork,
     &                    shp,           shgl,
     &                    shpb,          shglb, solinc(1,isclr+5))
c     
                  endif         ! end of scalar type solve
c     
c     
c.... end of the multi-corrector loop
c     
 1000             continue      !check this

               else             ! this is an update  (mod did not equal zero)
                  iupdate=icode/10 ! what to update
                  if(iupdate.eq.0) then !update flow  
                     call itrCorrect ( y, ac, yold, acold, solinc)
                     call itrBC (y,  ac,  iBC,  BC, iper, ilwork)
c Elaine-SPEBC
                     if((irscale.ge.0).and.(myrank.eq.master)) then
                        call genscale(y, x, iBC)
c                       call itrBC (y,  ac,  iBC,  BC, iper, ilwork)
                     endif
                  else          ! update scalar
                     isclr=iupdate !unless
                     if(iupdate.eq.nsclr+1) isclr=0
                     call itrCorrectSclr ( y, ac, yold, acold,
     &                                     solinc(1,isclr+5))
                     if (ilset.eq.2 .and. isclr.eq.2)  then
                        fct2=one/almi
                        fct3=one/alfi
                        acold(:,7) = acold(:,7) 
     &                             + (ac(:,7)-acold(:,7))*fct2
                        yold(:,7)  = yold(:,7)  
     &                             + (y(:,7)-yold(:,7))*fct3  
                        call itrBCSclr (  yold,  acold,  iBC,  BC, 
     &                                    iper,  ilwork)
                        ac(:,7) = acold(:,7)*(one-almi/gami)
                        y(:,7)  = yold(:,7)
                        ac(:,7) = zero 
                        if (ivconstraint .eq. 1) then
     &                       
c ... applying the volume constraint
c
                           call solvecon (y,    x,      iBC,  BC, 
     &                                    iper, ilwork, shp,  shgl)
c
                        endif   ! end of volume constraint calculations
                     endif
                     call itrBCSclr (  y,  ac,  iBC,  BC, iper, ilwork)
                  endif
               endif            !end of switch between solve or update
            enddo               ! loop over sequence in step
c     
c     Find the solution at the end of the timestep and move it to old
c  
c.... First to reassign the parameters for the original time integrator scheme
c
            if((iLSet.eq.2).and.(ilss.eq.1)) then 
               alfi =alfit
               gami =gamit
               almi =almit  
            endif          
            call itrUpdate( yold,  acold,   y,    ac)
            call itrBC (yold, acold,  iBC,  BC, iper,ilwork)  
c Elaine-SPEBC      
            if((irscale.ge.0).and.(myrank.eq.master)) then
                call genscale(yold, x, iBC)
c               call itrBC (y,  ac,  iBC,  BC, iper, ilwork)
            endif           
            do isclr=1,nsclr
               call itrBCSclr (yold, acold,  iBC, BC, iper, ilwork)
            enddo
c     
            istep = istep + 1
c     
c.... compute boundary fluxes and print out
c     
            lstep = lstep + 1
            if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
c
c here is where we save our averaged field.  In some cases we want to
c write it less frequently
               ntoutv=max(ntout,100) 
               if( (mod(lstep, ntoutv) .eq. 0) .and.
     &              ((irscale.ge.0).or.(itwmod.gt.0) .or. 
     &              ((nsonmax.eq.1).and.(iLES.gt.0))))
     &              call rwvelb  ('out ',  velbar  ,ifail)

               call Bflux  (yold,          acold,     x,
     &              shp,           shgl,      shpb,
     &              shglb,         nodflx,    ilwork)
            endif
c     
c.... end of the NSTEP and NTSEQ loops
c     

c...  dump TIME SERIES
            
            if (exts) then
               
               do jj = 1, ntspts
                  
                  if (numpe > 1) then
                     
                     soln = varts(jj, k)
                     asoln = abs(soln)
                     
c                     if(jj.eq.24) then
c                        write(*,*) soln
c                        write(*,*)"and..."
c                     endif

                     if (asoln.ne.zero) then
                        sgn = soln/asoln
                     else
                        sgn = 1
                     endif
                     
                     call MPI_ALLREDUCE ( asoln, asolng, 1, 
     &                    MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD,
     &                    ierr)
                     varts(jj, 1) = sgn * asolng
                     
                  endif
                  
                  if (myrank.eq.zero) then
                     ifile = 1000+jj
                     write(ifile,555) varts(jj, k)
                     call flush(ifile)
                  endif
                  
               enddo
               

               varts = zero  ! reset the array for next step
               
c$$$  do jjj=1,ntspts
c$$$  ndts=nodetimeseries(jjj)
c$$$  if(ndts.ne.0) then ! this processor has min
c$$$  ifile=jjj+1000+1000*myrank
c$$$  write(ifile,555) yold(ndts,4)
c$$$  call flush(ifile)
c$$$  endif
c$$$  enddo
 555           format(e18.11)
               
            endif
c$$$            
c$$$c...  dump DKE
c$$$
c$$$c
c$$$c  START OF DISTURBANCE BLOCK
c$$$c
c$$$
c$$$           
c$$$            if (numpe > 1) then
c$$$               call MPI_REDUCE (dke,dkesum,1,MPI_DOUBLE_PRECISION,
c$$$     &              MPI_SUM, master, MPI_COMM_WORLD,ierr)
c$$$            endif
c$$$            if (numpe.eq.1)
c$$$     &           write (776,1001) lstep+1, dke,dke
c$$$            
c$$$            if ((myrank .eq. master).and.(numpe > 1)) then 
c$$$               write (776,1001) lstep+1, dkesum,dkesum
c$$$c     
c$$$               call flush(776)
c$$$c     
c$$$            endif
c$$$            dke=0.0             ! we must zero it back out for next step
c$$$            
c$$$ 1001       format(1p,i6,8e13.5)
c$$$            
c$$$            
c$$$C
c$$$C     END OF DISTURBANCE BLOCK
c$$$            
            
 2000    continue
 2001    continue
c
c.... -------------------> error calculation  <-----------------
c 
            if(ierrcalc.eq.1) then
c
c compute average
c
               tfact=one/istep
               ybar =tfact*yold + (one-tfact)*ybar
c
c compute rms
c
               rerr(:, 7)=rerr(:, 7)+(yold(:,1)-ybar(:,1))**2
               rerr(:, 8)=rerr(:, 8)+(yold(:,2)-ybar(:,2))**2
               rerr(:, 9)=rerr(:, 9)+(yold(:,3)-ybar(:,3))**2
               rerr(:,10)=rerr(:,10)+(yold(:,4)-ybar(:,4))**2
c               rerr= 10.0*rerr
            endif
c            call restar('out ',rerr(:,1:5),ac) ! for debugging

         ttim(1) = REAL(secs(0.0)) / 100. - ttim(1)
         ttim(2) = secs(0.0)                 - ttim(2)

         tcorecp2 = REAL(secs(0.0)) / 100.
         tcorewc2 = secs(0.0)
         
         if (myrank .eq. master) then
            write(6,*) 'T(core) cpu-wc = ',tcorecp2-tcorecp1,
     &           tcorewc2-tcorewc1
         endif

c     call wtime

 3000 continue
c     
c.... ---------------------->  Post Processing  <----------------------
c     
c.... print out the last step
c     
      if ((irs .ge. 1) .and. ((mod(lstep, ntout) .ne. 0) .or.
     &     (nstp .eq. 0))) then
         if( (mod(lstep, ntoutv) .eq. 0) .and.
     &        ((irscale.ge.0).or.(itwmod.gt.0) .or. 
     &        ((nsonmax.eq.1).and.(iLES.gt.0))))
     &        call rwvelb  ('out ',  velbar  ,ifail)

         call Bflux  (yold,  acold,     x,
     &        shp,           shgl,      shpb,
     &        shglb,         nodflx,    ilwork)
      endif
c     
      if(ierrcalc.eq.1) then
c
c.... open the output file
c
         itmp = 1
         if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
         write (fmt1,"('(''error.'',i',i1,',1x)')") itmp
         write (fname1,fmt1) lstep
         
         fname1 = trim(fname1) // cname(myrank+1)
         open ( unit = 87, file = fname1, status = 'unknown',
     &        form = 'unformatted')
         
         write (87) machin, nshg, lstep
         write (87) rerr
         close (87)
      endif
c$$$      lstep = lstep +1
c$$$      call restar('out ',rerr(:,1:5),ac)  ! for debugging
c
c.... close history and aerodynamic forces files
c     
      if (myrank .eq. master) then
         close (ihist)
         close (iforce)
      endif
      close (iecho)
      if(iabc==1) deallocate(acs)
c
c.... end
c
        return
        end
