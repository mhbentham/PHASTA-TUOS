      subroutine forces (y ,  ilwork)
c
c----------------------------------------------------------------------
c
c This subroutine calculates and updates the aerodynamic forces
c
c----------------------------------------------------------------------
c
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
c     
      integer ilwork(nlwork)

      real*8  y(nshg,ndof)
      real*8  Forin(5), Forout(5),spmasss
      real*8  ftots(3,0:MAXSURF),ftot(3)
      real*8  spmasss1, spmasss2, Atots1, Atots2
c
c  START OF DISTURBANCE BLOCK
c
      dimension diste(1), distesum(1)
      distcalc=0  ! we would want to activate this for T-S instability studies
      if(distcalc.eq.1) then
        if(iter.eq.nitr) then  ! we have completed the last iteration
          if (numpe > 1) then
          call MPI_REDUCE (dke,dkesum,1,MPI_DOUBLE_PRECISION,
     &                     MPI_SUM, master, MPI_COMM_WORLD,ierr)
          endif
          if (numpe.eq.1)
     &      write (76,1000) lstep+1, dke,dke

          if ((myrank .eq. master).and.(numpe > 1)) then 
             write (76,1000) lstep+1, dkesum,dkesum
c
             call flush(76)
c     
          endif
          dke=0.0               ! we must zero it back out for next step

       endif
       endif

C
C END OF DISTURBANCE BLOCK
C

c
c.... -------------------->  Aerodynamic Forces  <----------------------
c
c.... output the forces and the heat flux
c
      if (numpe > 1) then
         Forin = (/ Force(1), Force(2), Force(3), HFlux, 
     &              entrop /)
         call MPI_ALLREDUCE (Forin(1), Forout(1),5,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         Force = Forout(1:3)
         HFlux = Forout(4)
         entrop= Forout(5)
      endif

      if (numpe > 1) then
         call MPI_ALLREDUCE (flxID(2,isrfIM), spmasss,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE (flxID(1,isrfIM), Atots,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE (flxID(3,:), Ftots(1,:),MAXSURF+1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE (flxID(4,:), Ftots(2,:),MAXSURF+1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE (flxID(5,:), Ftots(3,:),MAXSURF+1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         if (iLSet.ne.0) then
          call MPI_ALLREDUCE (flxIDsclr(1,isrfIM), Atots1,1,
     &         MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
          call MPI_ALLREDUCE (flxIDsclr(2,isrfIM), spmasss1,1,
     &         MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
          call MPI_ALLREDUCE (flxIDsclr(3,isrfIM), Atots2,1,
     &         MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
          call MPI_ALLREDUCE (flxIDsclr(4,isrfIM), spmasss2,1,
     &         MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         endif
      else
         Ftots=flxID(3:5,:)
         Atots=flxID(1,isrfIM)
         spmasss=flxID(2,isrfIM)
        if (iLSet.ne.0) then
         Atots1 = flxIDsclr(1,isrfIM)
         Atots2 = flxIDsclr(3,isrfIM)
         spmasss1 = flxIDsclr(2,isrfIM)
         spmasss2 = flxIDsclr(4,isrfIM)
        endif
      endif
      call sforce()

      ftot(1)=sum(Ftots(1,0:MAXSURF))
      ftot(2)=sum(Ftots(2,0:MAXSURF))
      ftot(3)=sum(Ftots(3,0:MAXSURF))
   
      if (isrfIM .ge. 0) then 
      if((spmasss.ne.0) .and. (matflg(5,1).ne.0))then ! we are doing
                                                      ! force adjustment 
         tmp=Dtgl*0.025         ! .025= 1/2/tmp1 when tmp1=20
         
         select case ( matflg(5,1) )
         case ( 1 )             ! standard linear body force
            fwall=Force(1)
            vlngth=xlngth
         case ( 2 )       
         case ( 3 )     
            fwall=Force(3)
            vlngth=zlngth
         end select

         write(222,*) spmasss
c         fnew=(vel-spmasss/(ro*Atots))*tmp + 
c     &        fwall/(vlngth*Atots)
         fold = datmat(1,5,1)
         fnew = fold + (vel-spmasss/Atots)*tmp*0.01
        if (iLSet.ne.0) then
         fnew = fold + (vel-spmasss1/Atots1)*tmp*0.01
        endif
         if(myrank.eq.0) then
            write(880,*) datmat(1,5,1),fnew
            call flush(880)
         endif
         datmat(1,5,1)=fnew     !* (one - tmp2) + datmat(1,5,1) * tmp2
      endif
      endif

      if (myrank .eq. master) then
         write (iforce,1000) lstep, (Force(i), i=1,nsd), 
     &        HFlux,spmasss,(ftot(i),i=1,nsd)
         call flush(iforce)
      endif
c
      if (pzero .eq. 1) then
c

c  also adjust the pressure here so that the mean stays at zero (prevent drift)
c
c  This is only valid if there are NO Dirichlet Boundary conditions on p!!!
c
c  (method below counts partition boundary nodes twice. When we get ntopsh
c  into the method it will be easy to fix this.  Currently it would require
c  some work and it is probably not worth it 
c
c  switched to avg of the on-processor averages to get around hierarchic
c  modes problem
c
c              pave=sum(y(1:numnp,1)) 
c              call MPI_ALLREDUCE (pave, pavet, 1,
c     &             MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
c     pavet=pavet/nshgt
         pave=sum(y(1:numnp,4)) 
         if (numpe .gt. 1) then
            xnpts=numnp
            call MPI_ALLREDUCE (pave, pavet, 1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (xnpts, xnptst, 1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            pavet=pavet/xnptst 
         else
            pavet = pave/numnp
         endif

         y(1:numnp,4)=y(1:numnp,4)-pavet
         pave=sum(y(1:numnp,4))
c         if(myrank.eq.0) then
c            write(90+myrank,*) pavet,pave
c            call flush(900)
c         endif
c
      endif

      return
c
 1000 format(1p,i6,8e13.5)
c
      end

      subroutine sforce()

      include "common.h"
      include "mpif.h"

      real*8  SFlux(5),SFluxg(5)
      character*20 fname1
      character*5  cname

      integer icalled, lskeep 
      data icalled /0/
      save icalled
      save lskeep

      if(icalled.eq.0) then
        lskeep = lstep
        icalled = 1
        do isrf = 0,MAXSURF
          if ( nsrflist(isrf).ne.0 ) then
            iunit=60+isrf
            fname1 = 'forces_s'//trim(cname(isrf))// cname(lstep) 
            open(unit=iunit,file=trim(fname1),status='unknown')
          endif
        enddo
      endif

      do isrf = 0,MAXSURF
        if ( nsrflist(isrf).ne.0 ) then
          SFlux = ( / flxID(1,isrf),                                    ! Area
     .                flxID(2,isrf),                                    ! Mass Flux
     .                flxID(3,isrf),                                    ! Force 1
     .                flxID(4,isrf),                                    ! Force 2
     .                flxID(5,isrf)/)                                   ! Force 3
          if ( numpe > 1 ) then
            call MPI_ALLREDUCE (SFlux, SFluxg,5,
     .        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            SFlux = SFluxg
          endif
          if(myrank.eq.0) then
            iunit=60+isrf
            write(iunit,"(i7,1p5e14.5)")lstep,SFlux
            call flush(iunit)
            if(((icalled .gt. 0) .and.
     &        (mod(lstep, ntout) .eq. 0))) then
              close(iunit)
              fname1 = 'forces_s'//trim(cname(isrf))// cname(lskeep)
              open(iunit,file=fname1, position='append')
             endif
          endif
        endif
      enddo

      return
      end
