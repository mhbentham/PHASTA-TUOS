        subroutine settauw (y,  x,  BC, 
     &                      ifath,   velbar)
c
c----------------------------------------------------------------------
c
c This routine computes the time varying viscous flux for turbulence wall
c  boundary elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use pointer_data
      use turbSA
      include "common.h"
c
      dimension y(nshg,ndofl),            x(numnp,nsd), 
     &          BC(nshg,ndofBC),
     &          ifath(numnp),             velbar(nfath,nflow),
     &          ull(nsd),                 trx(numnp,nsd),
     &          ullb(nsd),                dull(nsd),
     &          evisc(numnp)
      real*8 outvec(nshg,nsd+nsd+nsd)
      real*8 trx_tmp(1,nsd), trx_diff  !JF: array for trx bounding

c JMR Debug
      character*20 fname
      character*5 cname
c End JMR Debug
c
c calculate the traction magnitude for all nodes on the wall
c
      trx       = zero
      trx_tmp   = zero
      evisc     = zero
      effvisc   = zero
      rm        = datmat(1,2,1)

c	if (myrank.eq.master) write(*,*) 'nflow, nfath = ', nflow, nfath
c JMR Debug
c      fname = "node_coord_" // cname(myrank)
c      open (unit=700,file=fname,status='unknown');
c      write(700,1000)
c 1000 format ("  Node        x         y          z          otwn")
c 1001 format(i6,3(1x,e12.6),1x,i6)
c End JMR Debug
      do nodw = 1, numnp ! loop over nodes
c JMR Debug
c      write(700,1001) nodw, x(nodw,1), x(nodw,2), x(nodw,3), otwn(nodw)
c End JMR Debug
      if ( otwn(nodw).ne.nodw ) then ! wall node check
         if (iLES.gt.0.or.(iDNS.ne.0.and.abs(itwmod).eq.1)) then ! iLES + slip velocity DNS
            ull(:)=velbar(ifath(otwn(nodw)),2:4) ! recall old #ing
c       write(*,*) 'ull, otwn, ifath, y(1) = ', ull(1), 
c     1  otwn(nodw), ifath(otwn(nodw)), y(otwn(nodw),1)
         endif
         if (iRANS.lt.0) then
            ull(:)=y(otwn(nodw),1:3)
         endif
c	if (myrank.eq.master) write(*,*) 'iDNS = ', iDNS
         if (iDNS.gt.0.and.abs(itwmod).eq.2) then  ! wall (eff. viscosity) model w/o turb model
c for now will use the local off-the-wall velocity but will need to
c make this an averaged velocity in the future.
            ull(:)=y(otwn(nodw),1:3)
         end if
         ub=ull(1)*wnrm(nodw,1) !
     &     +ull(2)*wnrm(nodw,2) ! store u.n here for now
     &     +ull(3)*wnrm(nodw,3) !
c     u_parallel_to_boundary=u-(u.n)n  :
         ull(:)=ull(:)-ub*wnrm(nodw,:) ! ull in flow
c     u_b is || u_ll ||
         ub=sqrt(ull(1)**2+ull(2)**2+ull(3)**2)
c     perp d2wall is (x2-x1).n   :
         dw=abs(
     &         (x(otwn(nodw),1)-x(nodw,1))*wnrm(nodw,1)
     &        +(x(otwn(nodw),2)-x(nodw,2))*wnrm(nodw,2)
     &        +(x(otwn(nodw),3)-x(nodw,3))*wnrm(nodw,3)
     &         )
         if(ub.eq.0) then
            ut=zero
            twoub=zero
         else
            ut=utau(ub,dw,rm,nodw,x(nodw,:))   ! find utau
            twoub=-ut*ut/ub     ! find -tau_w/ub
         endif 
c
c for LES ull has the mean-parallel velocity vector.  We want the 
c instantaneous-parallel velocity vector
c
         if (iLES.gt.0.or.(iDNS.gt.0.and.abs(itwmod).eq.1)) then
            ullb=ull    ! saved the ubar
            ull(:)=y(otwn(nodw),1:3)
            ubn=ull(1)*wnrm(nodw,1) !
     &        +ull(2)*wnrm(nodw,2) ! store u.n here for now
     &        +ull(3)*wnrm(nodw,3) !
c
c     u_parallel_to_boundary=u-(u.n)n  :
c
            ull(:)=ull(:)-ubn*wnrm(nodw,:) ! ull in flow
c
c hack a limiter into this fluctuating vector. Early transients have
c huge differences from mean values
c
            dull=ull-ullb ! the current vector difference
            dullm=sqrt(dull(1)*dull(1)+dull(2)*dull(2)+dull(3)*dull(3))
     &           + 1.0e-9
            ullbm=ub ! mag of ullb still there
c
c limit the magnitude of the difference to a 40% change from the mean.
c if less than that already we will take the whole difference, otherwise
c only take a 40% change.
c

            dullmod=min(one,0.4*ullbm/dullm)
            ull=ullb+dullmod*dull
         endif

!         trx(nodw,:)=twoub*ull(:)
!c JF CHANGE
        trx_tmp(1,:) = twoub*ull(:)
        trx_diff = (trx_tmp(1,1)-trx(nodw,1))**2
     &           + (trx_tmp(1,2)-trx(nodw,2))**2
     &           + (trx_tmp(1,3)-trx(nodw,3))**2
        trx_diff = sqrt(trx_diff)
        trx_diff = trx_diff / sqrt(trx(nodw,1)**2
     &           + trx(nodw,2)**2 + trx(nodw,3)**2)

        if(trx_diff .lt. 0.1) then
           trx(nodw,:) = trx_tmp(1,:)
        else
           trx(nodw,:) = 0.90*trx(nodw,:) + 
     &                   0.10*trx_tmp(1,:)      !10% adjustment 

!           if(nodw.ne.1 .and. nodw.ne.numnp) then
!              trx(nodw,:) = 0.5*(trx(nodw-1,:) +
!     &                           trx(nodw+1,:) ) 
!           endif
        endif !trx_diff

!        if(myrank.eq.master.and.nodw.eq.18) write(*,*) trx(nodw,:)
!c JF CHANGE ENDS
c         trx(nodw,:)=twoub*ullb(:)    ! attempt to use a mean velocity (not instantenous)
c	if (myrank.eq.master.and.nodw.eq.18) write(*,*) 'nodw, otwn(), twoub = ', nodw, otwn(otwn), twoub

         if (abs(itwmod) .eq. 2) then  ! effective vicosity
           if((iRANS.lt.0).and.(itwmod.eq.-2)) then ! w/ RANS (redundant logic here)
            tauw=ut*ut
            BC(nodw,7)=tauw*dw/ub-rm
           else if ((iDNS.gt.0).and.(itwmod.eq.-2)) then ! DNS (redundant logic here)
              tauw=ut*ut
              effvisc(nodw)=tauw*dw/ub-rm
           else if ((iLES>0).and.(itwmod.eq.2)) then ! LES
c
c  mag of u instantaneous
c
            ullm=sqrt(ull(1)*ull(1)+ull(2)*ull(2)+ull(3)*ull(3))
            tauw=ut*ut*ullm/ub
            evisc(nodw)=tauw*dw/ub-rm
           endif
         endif
         if((itwmod.eq.-1).and.(iDNS.eq.0)) then ! slip-velocity RANS, not DNS!
            up=sqrt(
     &           y(nodw,1)**2   !
     &           +y(nodw,2)**2  ! flow is auto-|| at boundary
     &           +y(nodw,3)**2  !
     &           )/ut
            BC(nodw,7)=savarw(up,rm,saCv1,nodw,x(nodw,:))
         endif
      endif ! wallnode check
      enddo ! loop over nodes
c
c JMR Debug
c      close(700)
c End JMR Debug
c
c Write the traction vectors to a file restart.4077.n
c
c      if (myrank.eq.master) write(*,*) 'Preparing to print ' 
c      ilstep=4077
c      outvec(:,1:3)=trx(:,1:3)
c      outvec(:,4:6)=0
c      do i=1,numnp
c         if(otwn(i).ne. i ) outvec(i,4:6)=y(otwn(i),1:3)
c      enddo
c      outvec(:,7:9)=wnrm(:,1:3)
c      if (myrank.eq.master) write(*,*) 'Calling write_restart... '
c      call write_restart(myrank,ilstep,numnp,ndof,time,outvec,outvec)
c       if (myrank.eq.master) 
c     1 write(*,*) 'Traction dumped to restart.4077.*'
c
c Put traction calculations into BCB
c
      do iblk = 1, nelblb
         iel    = lcblkb(1,iblk)
         nenl   = lcblkb(5,iblk) ! no. of vertices per element
         nenbl  = lcblkb(6,iblk) ! no. of vertices per bdry. face
         ndofl  = lcblkb(8,iblk)
         npro   = lcblkb(1,iblk+1) - iel 
c For all elements in this block that lie on a wall, assign the traction
         do i=1,npro
            if(btest(miBCB(iblk)%p(i,1),4)) then ! wall elt
c         if (myrank.eq.master.and.i.eq.1) write(*,*) 'turb wall is here'
               do j = 1, nenbl
                  if(itwmod.eq.-2) then ! effective-viscosity
                     mBCB(iblk)%p(i,j,3:5)=0.0  ! set traction to zero - won't use it
                  endif
                  if((itwmod.eq.-1).or.(itwmod.eq.1)) then ! slip-velocity
c         if (myrank.eq.master.and.i.eq.1.and.j.eq.1) 
c     1     write(*,*) 'traction is set to ', trx(mienb(iblk)%p(i,j),:)
                     mBCB(iblk)%p(i,j,3:5)=trx(mienb(iblk)%p(i,j),:) ! set traction
                  endif
               enddo
            endif
         enddo
      enddo


      if(itwmod.eq.2) then   !effective viscosity
c
c For the elements which touch a modeled wall,
c modify the eddy viscosity at all quadrature points to be the average
c of the "optimal" nodal values.  This is an element-wise constant.
c
         do iblk = 1,nelblk
            nenl   = lcblk(5,iblk) ! no. of vertices per element
            iel    = lcblk(1,iblk)
            npro   = lcblk(1,iblk+1) - iel
            lcsyst = lcblk(3,iblk)
            ngauss = nint(lcsyst)
            do i=1,npro
               xnave=zero
               avevisc=zero
               do j=1,nenl
                  ev = evisc(mien(iblk)%p(i,j))
                  avevisc=avevisc + ev
                  if(ev.ne.zero) xnave=xnave+1
               enddo
               if(xnave.ne.0) mxmudmi(iblk)%p(i,1:ngauss)=avevisc/xnave
            enddo
         enddo
      endif
c
c.... end
c
        return
        end

        function utaul(u,y,rm)
        real*8 utaul
        utaul=.2
        return
        end

        function utau(u,y,rm,nodw,x)
c---------------------------------------------------------------------------
c
c  Computes the friction velocity, u_tau, that is consistent
c  with the off-the-wall parallel velocity (u) and normal distance (y).
c  The Spalding (1961) composite log-law formula is used within.
c
c---------------------------------------------------------------------------
        implicit none
        real*8 u,err,utau,yrmi,y,rm,yp,up,kup,f,dfds,rat,efac
        real*8 kup2,kup3,pt5,sxth,kappa
        integer iter
        integer nodw
        real*8 x(3)
        pt5=0.5
        sxth=0.1666666666666667
        err=1.0d-6
        utau=0.04
        yrmi=y/rm
        kappa=0.4
c$$$        B=5.5
        efac=0.1108 ! exp(-kappa*B)
c
c This is an iterative solve   
c    1 - guess u_tau
c    2 - compute difference between y+=y*u_tau/visc and that from Spalding (f term below)
c    3 - compute gradient of difference wrt u_tau (dfds term below)
c    4 - compute correction on u_tau: du_tau = df / (dfds)
c    5 - update u_tau = u_tau_old + correction
c    6 - iterate again if not within tolerance
c
        do iter=1,500
           yp=yrmi*utau
           up=u/utau
           kup=kappa*up
           kup2=kup*kup
           kup3=kup*kup2
           f=   up-yp+efac*(exp(kup)-1.0-kup - kup2*pt5 - kup3*sxth)
           dfds=up+yp+efac*(exp(kup)*kup-kup - kup2 - kup3*pt5) !/-utau in rat
           rat=f*utau/dfds ! this is -f/dfds but we did not do 1/(-utau) yet.
           utau=utau+rat
           if(abs(rat).le.err) goto 20
        enddo
c        write(*,*)'utau failed to converge at ',nodw,x(1),x(2),x(3)
c     1  , 'u,          dwallperp,      mu', u, y, rm, 
c     1  'dfds,         rat,         utau', dfds,rat,utau
c        stop
c        utau=0.0
c
c  if the above fails then try a simple no-slip traction
c
        utau=sqrt(u/y*rm)

 20     continue
        return
        end  

c$$$        function utau_log_layer_only(u,y,rm)
c$$$        real*8 u,err,utau,yrmi,y,rm,lnyp,f,dfds,rat
c$$$        err=1.0d-6
c$$$        utau=0.04
c$$$        yrmi=y/rm
c$$$        do iter=1,50
c$$$           lnyp=log(yrmi*utau)
c$$$           f=u-utau*(2.5*lnyp+5.2)
c$$$           dfds=-2.5*(lnyp+6.2)
c$$$           rat=-f/dfds
c$$$           utau=utau+rat
c$$$           if(utau.gt.0.5) then !flow is laminar for now give the parabolic
c$$$              utau=sqrt(2.0*rm*u)/y
c$$$              goto 20
c$$$           endif
c$$$              
c$$$           if(abs(rat).le.err) goto 20
c$$$        enddo
c$$$        write(*,*)'utau failed to converge',r,dfds,rat,utau,y,u,rm
c$$$        stop
c$$$ 20     continue
c$$$        return
c$$$        end  

        function savarw(up,rm,cv1,nodw,x)
        implicit none
        real*8 err,savarw,rm,up,cv1,f,dfds,rat,efac
        real*8 pt5,kappa,B,xmut,chi3,denom,cv1_3
        integer iter
        integer nodw
        real*8 x(3)
        pt5=0.5
        err=1.0d-6
        savarw=rm*cv1*1.2599 ! inflection point chi=cv1*cuberoot(2)
        kappa=0.4
c$$$        B=5.5
        efac=0.1108 ! exp(-kappa*B)
        xmut=rm*kappa*efac*(exp(kappa*up)-1.0-kappa*up
     &       -pt5*(kappa**2)*(up**2))
        do iter=1,500
           chi3=savarw/rm
           chi3=chi3*chi3*chi3
           cv1_3=cv1**3
           denom=chi3+cv1_3

           f=savarw*chi3/denom - xmut
           dfds=chi3*(chi3+4.0*cv1_3)/(denom**2)
           rat=-f/dfds
           savarw=savarw+rat
           if(abs(rat).le.err) goto 20
        enddo
        write(*,*)'savarw failed to converge at ',nodw,x(1),x(2),x(3)
c        write(*,*) 'uplus, mu', up, rm
c        write(*,*) 'dfds,        rat,        savarw', dfds, rat, savarw
        savarw=10000.0*rm
 20     continue
        return
        end  
