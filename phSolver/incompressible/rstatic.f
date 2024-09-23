        subroutine rstatic (res, y, Dy)
c
c----------------------------------------------------------------------
c
c This subroutine calculates the statistics of the residual.
c
c input:
c  res   (nshg,nflow)   : preconditioned residual
c
c output:
c  The time step, cpu-time and entropy-norm of the residual
c     are printed in the file HISTOR.DAT.
c  
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        dimension res(nshg,nflow),    mproc(1)
!SCATTER        dimension  rvec(numpe)
        dimension rtmp(nshg),        nrsmax(1)
!SCATTER        dimension irecvcount(numpe), resvec(numpe)

        real*8    y(nshg,ndof),    Dy(nshg,4)
c        integer tmrc
c
c$$$	ttim(68) = ttim(68) - tmr()
c
c.... compute max delta y
c
        rdy1 = zero
        rdy2 = zero
        rdy4 = zero
        rdy5 = zero
        call sumgatN( abs(gami*Delt(itseq) 
     &                * Dy(1:numnp,1:3)),3,rdy1, numnp)
        call sumgatN( abs( y(1:numnp,1:3)),3,rdy2,numnp)
        call sumgatN( abs(gami*alfi*Delt(itseq)
     &                * Dy(1:numnp,4)),1,rdy4,numnp)
        call sumgatN( abs( y(1:numnp,4)),  1,rdy5,numnp)
        rmaxdyU = rdy1/rdy2
        rmaxdyP = rdy4/rdy5
	
c
c..... Signal to quit if delta is very small. look in itrdrv.f for the
c      completion of the hack.
c
	if( rmaxdyU .lt. dtol(1) .and. rmaxdyP .lt. dtol(2)) then
           istop = 1000
        endif

        if (numpe == 1) nshgt=nshg   ! global = this processor
c
c
c.... ----------------------->  Convergence  <-------------------------
c
c.... compute the maximum residual and the corresponding node number
c
        rtmp = zero
        if (impl(itseq) .ge. 9) then
          do i = 1, nflow
            rtmp = rtmp + res(:,i)**2    ! only add continuity and momentum
          enddo
        endif

        call sumgat (rtmp, 1, resnrm)
        
        resmaxl = maxval(rtmp)

        irecvcount = 1
        resvec = resmaxl
        if (numpe > 1) then
!MR CHANGE
           if(impistat.eq.1) iAllR = iAllR+1
           if(impistat2.eq.1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
           if(impistat.eq.1) rmpitmr = TMRC()
!MR CHANGE END
           call MPI_ALLREDUCE (resvec, resmax, irecvcount,
     &          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
!MR CHANGE
           if(impistat.eq.1) rAllR = rAllR+TMRC()-rmpitmr
!MR CHANGE END
c           call MPI_REDUCE_SCATTER (resvec, resmax, irecvcount,
c     &          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
           if (resmax .eq. resvec ) then
c           if (resmax .eq. resvec(1) ) then
              mproc(1) = myrank
              nrsmax   = maxloc(rtmp)
           else
              mproc(1)  = -1
              nrsmax(1) = -1
           endif
           resvec = nrsmax(1)
!MR CHANGE
           if(impistat.eq.1) iAllR = iAllR+1
           if(impistat2.eq.1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
           if(impistat.eq.1) rmpitmr = TMRC()
!MR CHANGE END
           call MPI_ALLREDUCE (resvec, rvec, irecvcount,
     &          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
!MR CHANG
           if(impistat.eq.1) rAllR = rAllR+TMRC()-rmpitmr
!MR CHANGE END
c           call MPI_REDUCE_SCATTER (resvec, rvec, irecvcount,
c     &          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
           nrsmax = rvec
c           nrsmax = rvec(1)
           resvec = mproc(1)
!MR CHANGE
           if(impistat.eq.1) iAllR = iAllR+1
           if(impistat2.eq.1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
           if(impistat.eq.1) rmpitmr = TMRC()
!MR CHANGE END
           call MPI_ALLREDUCE (resvec, rvec, irecvcount,
     &          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
!MR CHANGE
           if(impistat.eq.1) rAllR = rAllR+TMRC()-rmpitmr
!MR CHANGE END
c           call MPI_REDUCE_SCATTER (resvec, rvec, irecvcount,
c     &          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
           mproc(1) = rvec
c           mproc = rvec(1)
      else
          resmax = resmaxl
          nrsmax = maxloc(rtmp)
          mproc(1) = 0
      endif
c
c.... correct the residuals
c
        if (loctim(itseq) .eq. 0) then
          resnrm = resnrm 
          resmax = resmax
        else
          resnrm = resnrm
          resmax = resmax
        endif
c
c.... approximate the number of entries
c
        !write(*,*) 'nshgt is ', nshgt   !MB
        totres = resnrm / abs(real(nshgt))
        totres = dsqrt(totres)
        resmax = dsqrt(resmax)
        if (resfrt .eq. zero) resfrt = totres
        jtotrs = int  ( 10.d0 * dlog10 ( totres / resfrt ) )
        jresmx = int  ( 10.d0 * dlog10 ( resmax / totres ) )
c     
c.... get the CPU-time
c
CAD        cputme = (second(0) - ttim(100))
        rsec=TMRC()
        cputme = (rsec - ttim(100))
c
c.... output the result
c
        if (numpe > 1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
        
        if (myrank .eq. master) then

!        print *, nshgt, totres, resmax, jtotrs, jresmx, resfrt, resnrm
c
c.... results of continuity and momentum 
c
           
           print 2000, lstep+1, delt(itseq), time+delt(itseq),
     &          cputme, totres, jtotrs, rmaxdyU,
     &          rmaxdyP,nrsmax,
     &          mproc(1)+1, jresmx, int(statsflow(4)),
     &          int(statsflow(1)), CFLfl_max, iCFLfl_maxelem
           write (ihist,2000) lstep+1, delt(itseq), time+delt(itseq), 
     &          cputme, totres, jtotrs, 
     &          rmaxdyU, rmaxdyP, nrsmax,
     &          mproc(1)+1,jresmx,int(statsflow(4)),
     &          int(statsflow(1)), CFLfl_max, iCFLfl_maxelem
           
           call flush(ihist)
        endif
        if(numpe>1) call MPI_BARRIER (MPI_COMM_WORLD,ierr)

c$$$	ttim(68) = ttim(68) + tmr()

c
c.... return
c
        return
c
 1000   format(1p,i6,5e13.5)
c 2000   format(1p,i6,e10.3,e10.3,2x,'(',i4,')',2x,e10.3,2x,e10.3,
c     &       2x,'<',i6,'-',i5,'|',
c     &       i4,'>', ' [', i4,' -',i4,']')
 2000   format(1p,i6,e10.3, e10.3,e10.3,e10.3,2x,'(',i6,')',2x,e10.3,
     &       2x,e10.3,2x,'<',i6,'-',i5,'|',
     &       i6,'>', ' [', i4,' -',i4,']',e10.3,i8)
 3000   format(1p,i6,e10.3,e10.3,3x,'(',i6,')',3x,'<',i6,'-',i5,'|',
     &       i6,'>', ' [', i4,' -',i4,' -',i4,']')

c
        end


        subroutine rstaticSclr (res, y, Dy, icomp)
c
c----------------------------------------------------------------------
c
c This subroutine calculates the statistics of the residual
c
c----------------------------------------------------------------------
c
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        dimension res(nshg)
        dimension rtmp(nshg)
        real*8    y(nshg,ndof),    Dy(nshg), nrm
c
        real*8 CFL_max_tmp
        integer iCFL_maxelem_tmp
c        integer tmrc
c
c.... compute max delta y
c
        rdy1 = zero
        rdy2 = zero
c
c.... normalize turbulence with molecular viscosity
c        
        if ( (icomp .eq. 6).and. (iRANS.eq.-1) ) then
           nrm = datmat(1,2,1)
        else 
           nrm = zero
        endif
        call sumgat( abs(gami*Delt(itseq)*Dy(:)),1,rdy1)
        call sumgat( abs( y(:,icomp)),1,rdy2)
        rmaxdyT = rdy1/(rdy2+nrm)
c
c.... compute the maximum residual and the corresponding node number
c
        rtmp = zero
        rtmp = rtmp + res**2 ! add temperature also
        call sumgat (rtmp, 1, resnrm)

        if (numpe == 1) nshgt=nshg ! global = this processor
        totres = resnrm / float(nshgt)
        totres = sqrt(totres)

c     redistance tolerance
         if (isclr.eq.2) then
              redist_toler_curr = totres
         endif

c	if (myrank.eq.master) then
c      write(*,*) "iSolvLSSclr1,2,ILSet = ", iSolvLSSclr1, iSolvLSSclr2, iLSet
c	end if

c     redistance CFL
         if ((iLSet.eq.2) .and. (isclr.gt.0)) then
           CFL_max_tmp = CFLls_max
           iCFL_maxelem_tmp = iCFLls_maxelem
         else
           CFL_max_tmp = CFLfl_max
           iCFL_maxelem_tmp = iCFLfl_maxelem
         endif
        

c        if (mod(impl(1),100)/10 .eq. 0) then  !not solving flow
           if (myrank .eq. master) then
c     
c.... get the CPU-time
c
              rsec=TMRC()
              cputme = (rsec - ttim(100))

           print 802, lstep+1, delt(itseq), time+delt(itseq),
     &                cputme, totres, rmaxdyT,
     &                int(statssclr(1)), CFL_max_tmp, iCFL_maxelem_tmp
           write (ihist,802) lstep+1, delt(itseq), time+delt(itseq),
     &          cputme, totres, rmaxdyT,
     &          int(statssclr(1)), CFL_max_tmp, iCFL_maxelem_tmp
           
               call flush(ihist)
           endif
c        else 
c           if (myrank .eq. master) then
c              print 803, totres, rmaxdyT, int(statssclr(1))
c              write(ihist,803) totres, rmaxdyT, int(statssclr(1))
c           endif
c        endif

        return
        
 802    format(1p,i6,e10.3,e10.3,e10.3,e10.3,10X,e10.3,
     &         31X,'[',i6,']',e10.3,i8)
 803    format(1p,16x,20x,e10.3,10x,e10.3,31X,'[',i10,']')
    
        end



