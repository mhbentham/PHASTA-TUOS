!----------------------------------------------------------------------
!
!       This file contains several subroutines which are dealing with
!       bubble tracking capability and advanced co-processing. The
!       subroutines are called in many subroutines, such as itrdrv.f,
!       elmgmr.f and so on. 
!
!       Jun Fang,                               Summer, 2015
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
!       Name              Description
!       procs_dataset   : data collected from a single processor
!       unive_dataset   : data collected from all processors
!       bub_info        : bubble information array for each element
!       avg_info        : bubble information array for each bubble
!       one_procs       : dummy of procs_dataset
!       all_procs       : dummy of unive_dataset
!       phi_min         : the min level set value in a bubble
!       phi_tmp         : dummy of phi_min
!       eq_rad          : equivalent radius for deformed bubble
!       i_mrk           : index of colored bubbles
!       elemvol_global  : array contains element volumes
!       elemvol_local   : array contains element volumes for each block
!       Spnt_min        : local shear calculation point with min y
!       coordinate
!       Spnt_max        : local shear calculation point with max y
!       coordinate
!       Spnt            : two points array used in local shear
!                         calculation
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
        subroutine CountIniBub(banma)
!----------------------------------------------------------------------
!       This subroutine is used to count the total number of bubbles in
!       the initial profile. 
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        USE, INTRINSIC :: ISO_C_BINDING !for calling C++ routines
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        dimension banma(nshg,1)
        real*8    local_bmmax, global_bmmax


        if((myrank.eq.master).and.(iClrLiq.eq.1)) then
           write(*,'(1x,A)') 'Tracking local liquid: Yes!'
           write(*,'(1x,A,ES15.7)') 'epsilonBT = ', epsilonBT
           write(*,'(1x,A,ES15.7)') 'phi_inner = ', phi_inner
           write(*,'(1x,A,ES15.7)') 'phi_outer = ', phi_outer
        endif

        local_bmmax     = 0.0d0
        global_bmmax    = 0.0d0
        do i = 1, nshg
           if(local_bmmax .lt. banma(i,1))
     &        local_bmmax = banma(i,1)
        enddo
        call MPI_ALLREDUCE (local_bmmax, global_bmmax, 1,
     &       MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

        i_num_bubbles = INT(global_bmmax)
        if(myrank.eq.master) write(*,*) 'Number of bubbles: ',
     &  i_num_bubbles

        XLEN = DomainSize(2)-DomainSize(1)
        YLEN = DomainSize(4)-DomainSize(3)
        ZLEN = DomainSize(6)-DomainSize(5)
!       The domain center coordinates will be used in bubble center
!       adjustment
!       (for those crossing periodic boundaries)
        Xmid = 0.5d0* ( DomainSize(2) + DomainSize(1))
        Ymid = 0.5d0* ( DomainSize(4) + DomainSize(3))
        Zmid = 0.5d0* ( DomainSize(6) + DomainSize(5))

        if(myrank.eq.master) then
           write(*,*)
           write(*,*) 'Note BT requires fixed epsilon_ls'
           write(*,*) 'Note BT requires d2wall information'
           write(*,*)
           write(*,'(1x,A)') 'Domain dimensions: '
           write(*,'(1x,A,2ES15.7)') 'x range: ', DomainSize(1:2)
           write(*,'(1x,A,2ES15.7)') 'y range: ', DomainSize(3:4)
           write(*,'(1x,A,2ES15.7)') 'z range: ', DomainSize(5:6)

        endif

        end     !CountIniBub ends 
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
        subroutine OpenBubFiles(nBubbleLTS, i_num_bubbles, 
     &                          C_int_adjust)
!----------------------------------------------------------------------
!       Open files for saving the data about bubble properties or
!       behaviors (which is part of bubble tracking algorithm)
!----------------------------------------------------------------------

        integer:: i_num_bubbles
        INTEGER:: file_number, il, irstart, ierror
        integer:: nl_alpha, inistep
        integer,allocatable :: alpha01(:,:)
        real*8  C_int_adjust
        real*8, allocatable :: alpha02(:,:)
        logical:: ObjExt1, ObjExt2

        C_int_adjust = 0.0d0
        write(*,*)
        inquire(file='../bubStat/',exist=ObjExt1)
        if(ObjExt1.eqv..False.) then
           call system('mkdir ../bubStat')
           write(*,*) 'bubStat not found, plz create one manually!'
           write(*,*)
        else
           write(*,*) 'bubStat already exists'
           write(*,*)
        end if
        write(*,*) 'In OpenBubFiles, i_num_bubbles = ', i_num_bubbles
!       open file to record void fraction during simulations and trim
!       the existing files such that new data can fit without any
!       overlap with the previous data
        inquire(file='../bubStat/alpha.dat',exist=ObjExt2)
        if(ObjExt2.eqv..False.) then   !no data file found, create them

        open(unit=740,file='../bubStat/alpha.dat',status="replace",
     &       action="write",iostat=ierror)
        if(ierror /= 0) STOP "Error creating file 740"
        close(740)

!----------------------------------------------------------------------
        else                            !old data files found, trim them
                                        !if necessary
        call CountAlphaLines(nl_alpha,inistep)
        write(*,*) 'Number of timesteps recorded = ', nl_alpha
        write(*,*) 'initial timestep in alpha.dat = ',inistep

        open(unit=72,file='numstart.dat',status='old')
        read(72,*) irstart
        close(72)

        write(*,*) 'Time step when processing bubble files is', 
     &             irstart+1

!----------------------------------
        if(irstart+1 .lt. (inistep+nl_alpha) 
     &                  .and. irstart+1.ne.inistep) then 
          write(*,*) 'Bubble files are being trimed!'

          allocate (alpha01(irstart+1-inistep,3))
          allocate (alpha02(irstart+1-inistep,3))

          open(unit=740,file='../bubStat/alpha.dat',status="old",
     &       action="read",iostat=ierror)
          if(ierror /= 0) STOP "Error reading file 740"
          do il = 1, (irstart+1-inistep)
             read(740,813) alpha01(il,1:3), alpha02(il,1:3)
          enddo
          close(740)
          open(unit=740,file='../bubStat/alpha.dat',status="replace",
     &       action="write",iostat=ierror)
          if(ierror /= 0) STOP "Error creating file 740"
          do il = 1, (irstart+1-inistep)
             write(740,813) alpha01(il,1:3), alpha02(il,1:3)
!       Use the C_int_adjust from previous simulation that can better
!       preserve void fraction 
             if(il.eq.(irstart+1-inistep)) C_int_adjust=alpha02(il,1)
          enddo
          close(740)

          deallocate(alpha01)
          deallocate(alpha02)
!
!----------------------------------
        elseif(irstart+1.eq.inistep) then
!       restart the case from the every beginning so we can totally
!       replace the previous data files        
        open(unit=740,file='../bubStat/alpha.dat',status="replace",
     &       action="write",iostat=ierror)
        if(ierror /= 0) stop "Error creating file 740"
        close(740)

        endif   !irstart .lt. (inistep+nl_alpha)
        endif   !existence of alpha.dat

 813  format(I6, 1x, I4, 1x, I4, 3x, ES14.7, 1x, ES14.7, 1x, ES14.7)

        end     !OpenBubFiles ends
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
        subroutine BubASSY()
!----------------------------------------------------------------------
!       Called in elmgmr.f
!       This subroutine is used to assembly bubbles' information after
!       the loop over blocks on each processor
!
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        integer i_mrk, ib

        do i = 1, npro
           i_mrk = INT(bub_info(i,11))
!       The marker re-arrangement when auxillary index is seeded during
!       bubble breakup recognition
           if(iBK.eq.1 .and. i_mrk.gt.i_num_bubbles) then
              do ib = 1,i_num_bubbles
                 if(i_mrk.eq.int(breakupSeeder(ib,6)))then
                    i_mrk = int(breakupSeeder(ib,1))
!       j=15: phasic volume occupied by auxiliary ID
                    procs_dataset(i_mrk,15) = procs_dataset(i_mrk,15)
     &                                      + bub_info(i,4)
                    exit
                 endif
              enddo
           endif
!       Here starts the extraction part
           if (i_mrk .gt. 0) then
!       j=1-3: x coord, y coord, z coord (weighted by element volume)
              do j = 1, 3
              procs_dataset(i_mrk,j) = procs_dataset(i_mrk,j)
     &                               + bub_info(i,j)*bub_info(i,4)
              enddo
!       j=4-5: bubble element volume, gas mass
              do j = 4, 5
              procs_dataset(i_mrk,j) = procs_dataset(i_mrk,j)
     &                               + bub_info(i,j)
              enddo
!       j=6-8: local liquid velocity in x1, x2, x3 directions
              do j = 6, 8
              procs_dataset(i_mrk,j) = procs_dataset(i_mrk,j)
     &                               + bub_info(i,j)*bub_info(i,17)
              enddo
!       j=9: minimum level-set value in bubble region
              if(bub_info(i,10).lt.procs_dataset(i_mrk,9))
     &        procs_dataset(i_mrk,9) = bub_info(i,10)
!       j=10: the highest process rank the bubble elements exist
              if(bub_info(i,16).gt.procs_dataset(i_mrk,10))
     &           procs_dataset(i_mrk,10) = bub_info(i,16)
!       j=11: summation of element volume in local liquid shell
              procs_dataset(i_mrk,11) = procs_dataset(i_mrk,11)
     &                                + bub_info(i,17)
!       j=12-14: bubble velocity in x1, x2, x3 directions 
              do j = 12, 14
              procs_dataset(i_mrk,j) = procs_dataset(i_mrk,j)
     &                               + bub_info(i,j)
              enddo
!===============================================================================
!-------------------------------------------------------------------------------
!       This part is for temperature collection
!       j=16: element count in bubble region
!       j=18: temperature
!       j=19: temperature gradient
!       j=20: The number of element in the outer shell
!       j=21: The number of element in the inner shell
!-------------------------------------------------------------------------------          
              do j=18, 19
              procs_dataset(i_mrk,j) = procs_dataset(i_mrk,j)
     &                                 + bub_info(i,j)
              enddo
!              if(bub_info(i,18).ne.0.0d0) then
!              write(*,*)'procs_dataset(1,18)', procs_dataset(1,18)
!              endif
              do j = 20, 22
              procs_dataset(i_mrk,j) = procs_dataset(i_mrk,j)
     &                                 + bub_info(i,j)
              enddo
!              if(procs_dataset(i_mrk,21).ne.0.0d0)
!     &        write(*,*) 'procs_dataset(i_mrk,21)= ',
!     &                    procs_dataset(i_mrk,21)
!              if(procs_dataset(i_mrk,20).ne.0.0d0)
!     &        write(*,*) 'procs_dataset(i_mrk,20)= ',
!     &                    procs_dataset(i_mrk,20)
!              write(*,*)'bub_info 20,21',bub_info(i,20),bub_info(i,21)
!-------------------------------------------------------------------------------
!===============================================================================

!       j=16: element count in bubble region
              if(bub_info(i,4).gt.0.0d0)
     &        procs_dataset(i_mrk,16) = procs_dataset(i_mrk,16) + 1.0d0

!       Find out the upper and lower bounds of element coordinates
!       inside each bubble
              IF(bub_info(i,4).ne.0.0d0) THEN   !Inside bubble
              if(bub_info(i,1).lt.procs_coordDn(i_mrk,1)) 
     &           procs_coordDn(i_mrk,1) = bub_info(i,1)
              if(bub_info(i,2).lt.procs_coordDn(i_mrk,2))
     &           procs_coordDn(i_mrk,2) = bub_info(i,2)
              if(bub_info(i,3).lt.procs_coordDn(i_mrk,3))
     &           procs_coordDn(i_mrk,3) = bub_info(i,3)

              if(bub_info(i,1).gt.procs_coordUp(i_mrk,1))
     &           procs_coordUp(i_mrk,1) = bub_info(i,1)
              if(bub_info(i,2).gt.procs_coordUp(i_mrk,2))
     &           procs_coordUp(i_mrk,2) = bub_info(i,2)
              if(bub_info(i,3).gt.procs_coordUp(i_mrk,3))
     &           procs_coordUp(i_mrk,3) = bub_info(i,3)
              ENDIF

!       Find out the points with min and max d2wall and the associated
!       velocities for each bubble
              if(bub_info(i,6).ne.0.0d0) then
                 if(bub_info(i,9).lt.Shear_NodeMin(i_mrk,1)) then
                    Shear_NodeMin(i_mrk,1) = bub_info(i,9)         !min d2w
                    Shear_NodeMin(i_mrk,2) = bub_info(i,6)         !x vel
                 endif
                 if(bub_info(i,9).gt.Shear_NodeMax(i_mrk,1)) then
                    Shear_NodeMax(i_mrk,1) = bub_info(i,9)         !max d2w
                    Shear_NodeMax(i_mrk,2) = bub_info(i,6)         !x vel
                 endif
              endif

           endif !i_mrk
        enddo !npro

        end     !BubASSY = bubble assembly
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
        subroutine BubMPIprocess()
!----------------------------------------------------------------------
!       Called in elmgmr.f
!       This subroutine is dealing with the MPI processing of bubble
!       information from different processors
!
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        real*8 one_procs,       all_procs
        real*8 phi_tmp,         phi_min,        phi_max
        real*8 eq_rad
        real*8 Shear_dummy01(2), Shear_dummy02(2), Shear_dummy03
        real*8 diffinY
        real*8 bubble_vol_temp,  numshell_temp1,   numshell_temp2,
     &         bubble_tempG_temp, bubbleID_temp

        allocate ( unive_dataset(i_num_bubbles, 22) )
        allocate ( unive_coordDn(i_num_bubbles, 3)  )
        allocate ( unive_coordUp(i_num_bubbles, 3)  )
        allocate ( bubbl_coordDf(i_num_bubbles, 3)  )
        allocate ( Shear(i_num_bubbles,4)           )

        unive_dataset   = zero 
        unive_coordDn   = zero
        unive_coordUp   = zero
        bubbl_coordDf   = zero
        Shear           = zero

        if (numpe .gt. 1) then
           do i = 1, i_num_bubbles
              do j = 1, 22
                 if (j.eq.9) then
!       Find out the minimum value for each bubble
                 phi_tmp = procs_dataset(i,j)
                 call MPI_ALLREDUCE (phi_tmp, phi_min, 1,
     &                MPI_DOUBLE_PRECISION, MPI_MIN, 
     &                MPI_COMM_WORLD, ierr)
                 unive_dataset(i,j) = phi_min
                 elseif(j.eq.10) then   !find max. 
                 phi_tmp = procs_dataset(i,j)
                 call MPI_ALLREDUCE (phi_tmp, phi_max, 1,
     &                MPI_DOUBLE_PRECISION, MPI_MAX,
     &                MPI_COMM_WORLD, ierr)
                 unive_dataset(i,j) = phi_max
                 else
                 one_procs = procs_dataset(i,j)
                 call MPI_ALLREDUCE (one_procs, all_procs, 1,
     &                MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
                 unive_dataset(i,j) = all_procs
                 endif
              enddo !index j
!------------------------------------------------------------------------
!       Here, we determine the x bounds, y bounds, z bounds for each
!       bubble
!------------------------------------------------------------------------
              do j = 1, 3
                 one_procs = procs_coordDn(i,j)
                 call MPI_ALLREDUCE (one_procs, all_procs, 1,
     &                MPI_DOUBLE_PRECISION, MPI_MIN,
     &                MPI_COMM_WORLD, ierr)
                 unive_coordDn(i,j) = all_procs

                 one_procs = procs_coordUp(i,j)
                 call MPI_ALLREDUCE (one_procs, all_procs, 1,
     &                MPI_DOUBLE_PRECISION, MPI_MAX,
     &                MPI_COMM_WORLD, ierr)
                 unive_coordUp(i,j) = all_procs
!       This range of elments coordinates inside a bubble will be used
!       is crossing the periodic planes
                 bubbl_coordDf(i,j) = unive_coordUp(i,j) -
     &                                unive_coordDn(i,j) 
              enddo
!------------------------------------------------------------------------
!       Find the two points to calculate local shear
!------------------------------------------------------------------------
!       The one with minimum d2wall
              Shear_dummy01(1) = Shear_NodeMin(i,1)
              Shear_dummy01(2) = myrank  ! myrank is to be coerced to a real*8

              call MPI_ALLREDUCE (Shear_dummy01, Shear_dummy02, 1,
     &                MPI_2DOUBLE_PRECISION, MPI_MINLOC,
     &                MPI_COMM_WORLD, ierr)

              Shear(i,1) = Shear_dummy02(1)
              Shear_dummy03 = Shear_NodeMin(i,2)
!       Broadcast the x vel (for local liq) associated with min d2wall
              call MPI_Bcast(Shear_dummy03,1,MPI_DOUBLE_PRECISION,
     &                       INT(Shear_dummy02(2)),MPI_COMM_WORLD,ierr)
              Shear(i,2) = Shear_dummy03
!       The one with maximum d2wall
              Shear_dummy01(1) = Shear_NodeMax(i,1)
              Shear_dummy01(2) = myrank  ! myrank is coerced to a real*8

              call MPI_ALLREDUCE (Shear_dummy01, Shear_dummy02, 1,
     &                MPI_2DOUBLE_PRECISION, MPI_MAXLOC,
     &                MPI_COMM_WORLD, ierr)

              Shear(i,3) = Shear_dummy02(1)
              Shear_dummy03 = Shear_NodeMax(i,2)
!       Broadcast the x vel (for local liq) associated with max d2wall
              call MPI_Bcast(Shear_dummy03,1,MPI_DOUBLE_PRECISION,
     &                       INT(Shear_dummy02(2)),MPI_COMM_WORLD,ierr)
              Shear(i,4) = Shear_dummy03

           enddo !i_num_bubbles
!------------------------------------------------------------------------
!       Calculate the average properties
!------------------------------------------------------------------------

        if(myrank.eq.master) then
        avg_info = zero

        do i = 1, i_num_bubbles
!       The gas entity smaller than 5x5x5 is not considered as
!       a valid bubble
           if(unive_dataset(i,16).lt.125) cycle
!       The average coordinates for each bubble
           if(unive_dataset(i,4).gt.0.0d0) then
           avg_info(i,1:3)  = unive_dataset(i,1:3)  /unive_dataset(i,4)
           avg_info(i,11:13)= unive_dataset(i,12:14)/unive_dataset(i,4)
           endif
           avg_info(i,14)   = unive_dataset(i,10)

           eq_rad = (3.0*unive_dataset(i,4)/(4.0*pi))**(1.0/3.0)
           avg_info(i,4)    = eq_rad
!       The bubble mass
           avg_info(i,5)    = unive_dataset(i,5)

!       The deformability factor is defined by the ratio b/w min level
!       set value with equivalent radius, which is 1.0 for spherical
!       bubbles ideally.
           if(eq_rad.ne.0.0d0)
     &     avg_info(i,6)    = abs(unive_dataset(i,9))/eq_rad

!       The local liquid velocity components near the bubble
           if(unive_dataset(i,11).gt.0.0d0)
     &     avg_info(i,7:9)  = unive_dataset(i,6:8)/unive_dataset(i,11)

!       The local liquid shear around the bubble
           diffinY = Shear(i,3) - Shear(i,1)
           if(diffinY.ne.0.0d0)
     &     avg_info(i,10)   = (Shear(i,4) - Shear(i,2))/diffinY

           if(iBK.eq.1) then
           breakupSeeder(i,4) = unive_dataset(i,4)
           breakupSeeder(i,5) = unive_dataset(i,15) 
           endif
!       Adjust the bubble centers coordinates for those crossing
!       periodic planes (1% tolerance is allowed)
           if(abs(bubbl_coordDf(i,1)-XLEN).lt.0.01d0*XLEN) then
              if(avg_info(i,1) .ge. Xmid ) then
                 avg_info(i,1) = DomainSize(2) - eq_rad * 
     &          (avg_info(i,1) - Xmid)/(0.5d0*XLEN-eq_rad)
              else
                 avg_info(i,1) = DomainSize(1) + eq_rad *
     &          (Xmid - avg_info(i,1))/(0.5d0*XLEN-eq_rad)
              endif
           endif

           if(abs(bubbl_coordDf(i,2)-YLEN).lt.0.01d0*YLEN) then
              if(avg_info(i,2) .ge. Ymid ) then
                 avg_info(i,2) = DomainSize(4) - eq_rad *
     &          (avg_info(i,2) - Ymid)/(0.5d0*YLEN-eq_rad)
              else
                 avg_info(i,2) = DomainSize(3) + eq_rad *
     &          (Ymid - avg_info(i,2))/(0.5d0*YLEN-eq_rad)
              endif
           endif

           if(abs(bubbl_coordDf(i,3)-ZLEN).lt.0.01d0*ZLEN) then
              if(avg_info(i,3) .ge. Zmid ) then
                 avg_info(i,3) = DomainSize(6) - eq_rad *
     &          (avg_info(i,3) - Zmid)/(0.5d0*ZLEN-eq_rad)
              else
                 avg_info(i,3) = DomainSize(5) + eq_rad *
     &          (Zmid - avg_info(i,3))/(0.5d0*ZLEN-eq_rad)
              endif
           endif
!-------------------------------------------------------------------------------
!       This part is for temperature gradient collection
!                                               -Mengnan Li
!-------------------------------------------------------------------------------
!... Store temperature and temperature gradient around the bubble
!           if(unive_dataset(i,18).ne.0.0)then
           avg_info(i,18) = unive_dataset(i,18)/unive_dataset(i,20)
           avg_info(i,19) = unive_dataset(i,19)/unive_dataset(i,20)
!... Store the number of cell outside the bubble
           avg_info(i,20) = unive_dataset(i,20)
           avg_info(i,21) = unive_dataset(i,21)
!... Store the total volume of the outside bubble collection shell
           avg_info(i,22) = unive_dataset(i,22)
!           endif
        if(bubboil.eq.1.0.or.bubgrow.eq.1.0) then
           write(*,*)'No.',i,'Bubble'
           write(*,*)'# of elem(out),# of elem(in)',
     &               avg_info(i,20),avg_info(i,21)
        endif
        if(bubboil.eq.1.0.or.bubgrow.eq.1.0) then
           write(*,*)'temp,tempG',
     &               avg_info(i,18),avg_info(i,19)
        endif
!-----------------------------------------------------------------------
!        write(*,*)'breakupSeeder in bubMPI =', breakupSeeder(i,:)
        enddo  !i_num_bubbles
!------------------------------------------------------------------------

        endif !myrank

        else
           write(*,*) 'Bubble tracking was not coded for 1-procs case!'
        endif
       if (numpe.gt.1) then
!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       do i=1,i_num_bubbles

        if(myrank.eq.master)then
        bubble_vol_temp   =  unive_dataset(i,4)
        bubble_tempG_temp =  avg_info(i,19)
        numshell_temp1    =  avg_info(i,20)
        numshell_temp2    =  avg_info(i,21)
        bubbleID_temp     =  unive_dataset(i,9)
        endif
!        write(*,*) 'I am here', bubble_vol(i)
          call MPI_Bcast(bubble_tempG_temp,1,MPI_DOUBLE_PRECISION,
     &                       master,MPI_COMM_WORLD,ierr)
!          for average temperature gradient
          call MPI_Bcast(numshell_temp1,1,MPI_DOUBLE_PRECISION,
     &                       master,MPI_COMM_WORLD,ierr)
!          for number of element outside the bubble
          call MPI_Bcast(numshell_temp2,1,MPI_DOUBLE_PRECISION,
     &                       master,MPI_COMM_WORLD,ierr)
!          for number of element inside the bubble
          call MPI_Bcast(bubble_vol_temp,1,MPI_DOUBLE_PRECISION,
     &                   master,MPI_COMM_WORLD,ierr)
!          for each bubble volume
          call MPI_Bcast(bubbleID_temp,1,MPI_DOUBLE_PRECISION,
     &                   master,MPI_COMM_WORLD,ierr)

          bubble_vol(i)  =  bubble_vol_temp
          bubble_tempG(i)=  bubble_tempG_temp
          numshell(i,1)  =  numshell_temp1
          numshell(i,2)  =  numshell_temp2
          bubbleID(i)    =  bubbleID_temp
!        write(*,*)'bubble_tempG',bubble_tempG(i)
       enddo
!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       endif

        deallocate ( unive_dataset )
        deallocate ( unive_coordDn ) 
        deallocate ( unive_coordUp ) 
        deallocate ( bubbl_coordDf ) 
        deallocate ( Shear         )  


        end     !BubMPIprocess ends
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
        subroutine banmaUpdate(xl, yl, banma, ien, bml)
!----------------------------------------------------------------------
!       Called in asigmr.f
!       This subroutine will update the marker field in each time
!       iteration.
!       Banma is the Chinese pinyin translation of Zebra, an African
!       wild horse with black-and-white stripes (markers) and an
!       erect mane.
!       The code will first find out the max of banma values in each
!       element and update the rest based on their levelset value on a
!       certain processor, the markers are updated from one element to
!       another, and this update is conducting samutaniously on
!       different processors
!
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"

        dimension banma(nshg,1),              bml(npro,nshl,1),
     &            yl(npro,nshl,ndofl),        xl(npro,nenl,nsd)
        real*8  bmtmp, xltmp(3)
        real*8  bmmax, Mrk_dist, Bub_radi
        integer:: imrk

!**********************************************************************
        IF(lstep.eq.ts_hold) THEN
!c... localization of marker field
        call local (banma,  bml,     ien,     1,      'gather  ')
        do i = 1,npro
           bmmax = 0.0d0
           do n = 1,nshl
              if(bml(i,n,1).gt.bmmax) bmmax = bml(i,n,1)
           enddo
           do n = 1,nshl
!c... if the nodal point is in the liquid
              if(yl(i,n,6).gt.0.0d0) then
                bml(i,n,1) = 0.0d0
!c... for the nodal points in the gas and near interface region...
              else
                bml(i,n,1) = bmmax
              endif
           enddo
        enddo
!c.... assemble the marker field
        call local (banma,    bml,   ien,    1,  'globaliz')
!**********************************************************************
        ELSEIF(lstep.ne.ts_hold.and.iClrLiq.eq.1) THEN
!       Do not use this part in the very first timestep because the
!       bub_cent is needed 
!       Localization of marker field
        call local (banma,  bml,     ien,     1,      'gather  ')
!       clean up ID field outside the bubble core region
        do i = 1,npro
           do n = 1,nshl
              if(yl(i,n,6).gt.0.0d0) bml(i,n,1) = 0.0d0
           enddo
        enddo

        do i = 1,npro
           bmmax = 0.0d0
           xltmp = zero
           do n = 1,nshl
              if(bml(i,n,1).gt.bmmax) bmmax = bml(i,n,1)
           enddo
!       Update the maker field differently in bubble elments and near
!       interface shell
           do n = 1,nshl
!       If the node is in the bubble region
              if(yl(i,n,6).lt.0.0d0) then
                bml(i,n,1) = bmmax
!       If the node is in the near interface liquid shell
              elseif(yl(i,n,6).lt.phi_outer.and.
     &               yl(i,n,6).gt.phi_inner) then
                xltmp(1:3) = xl(i,n,1:3)
                call ReColor(xltmp,bmtmp)
                bml(i,n,1) = bmtmp

              endif

           enddo        !nshl
        enddo           !npro

c.... assemble the marker field
        call local (banma,    bml,   ien,    1,  'globaliz')
!**********************************************************************
        ELSEIF(lstep.ne.ts_hold.and.iClrLiq.eq.0) THEN
        call local (banma,  bml,     ien,     1,      'gather  ')
        do i = 1,npro
           bmmax = 0.0d0
           xltmp = zero
           do n = 1,nshl
              if(bml(i,n,1).gt.bmmax) bmmax = bml(i,n,1)
           enddo
!       Update the maker field differently in bubble elments and near
!       interface shell
           do n = 1,nshl
!       If the node is in the liquid (eg. levelset > epsilon)
              if(yl(i,n,6).gt.0.0d0) then
                bml(i,n,1) = 0.0d0
!       If the node is in the bubble region
              else
                bml(i,n,1) = bmmax
              endif

           enddo        !nshl
        enddo           !npro

c.... assemble the marker field
        call local (banma,    bml,   ien,    1,  'globaliz')

        ENDIF   !ts_hold
        end     !banmaUpdate ends
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
        subroutine BubCollect(u1,       u2,     u3,     Sclr, dist2w,
     &                        xx,       yl,     bml,    elemvol_local,
     &                        rho,      Tempb,  gytemp,  CurvInfo)
!----------------------------------------------------------------------
!       Called in e3ivar.f
!       This subroutine is dealing with the bubble information
!       collection at the very bottom level.
!       In bub_info(i_num_bubbls, 17)
!       bubble-wise: x,y,z coord, vel, elem vol, mass, levelset;
!       local liq  : x vel, y vel, z vel, d2wall;
!       BT field   : marker(bubble ID)
!       Interface curvature is weighted by volume of the interface
!       element. 
!----------------------------------------------------------------------
        use  bub_track
        include "common.h"

        dimension u1(npro),     u2(npro),  u3(npro)
        dimension dist2w(npro), rho(npro), CurvInfo(npro) 
        dimension yl(npro,nshl,ndof),      xx(npro,nsd),
     &            Sclr(npro)
        dimension bml(npro,nshl,1)
        real*8    elemvol_local(ibksiz)
        real*8    bmmax 
        real*8    denswght
        dimension Tempb(npro), gytemp(npro,5)
        real*8    Tempb, gytemp, epsilon_lst_tmp
                  ! Temperature Collection Mengnan 9/6/15

        rholiq=datmat(1,1,1)
        rhogas=datmat(1,1,2) 

        bub_info = zero

        do i = 1, npro
!... collect LS value & markers for bubble and liquid shell
           if(Sclr(i) .le. 3.0d0*epsilonBT) then
              if(rholiq.eq.rhogas) then
                 denswght = 1.0d0
              else
                 denswght  =(rholiq-rho(i))/(rholiq-rhogas)
              endif
              bmmax     = 0.0
              do n = 1, nshl
                if(bml(i,n,1).gt.bmmax) bmmax = bml(i,n,1)
              enddo
              bub_info(i,11) = bmmax
!       Find out the minimum level set value
              do n = 1, nshl
                if(yl(i,n,6).lt.bub_info(i,10))
     &             bub_info(i,10) = yl(i,n,6)
              enddo
!       collect the local liquid velocity and the y coord of liquid
!       shell elments around the bubble
            if(Sclr(i).gt.2.0d0*epsilonBT)then
              elseif(Sclr(i).gt.epsilonBT) then
                bub_info(i,6)  = u1(i)
                bub_info(i,7)  = u2(i)
                bub_info(i,8)  = u3(i)
                bub_info(i,9)  = dist2w(i)
                bub_info(i,17) = elemvol_local(i) 
              elseif(Sclr(i).gt.0.0d0) then
                bub_info(i,5)  = elemvol_local(i)*denswght*rhogas
                bub_info(i,15) = CurvInfo(i)
!----------------------------------------------------------------------------
              if(bubboil.eq.1.or.bubgrow.eq.1)then
                bub_info(i,20)= 1.0d0
                bub_info(i,18)= Tempb(i)
                bub_info(i,19)= gytemp(i,5)
                bub_info(i,22)= elemvol_local(i)
!                if(isnan(gytemp(i,5)))then
!                     write(*,*)i
!                else
!                     write(*,*)'all good'
!                endif   
              endif
!----------------------------------------------------------------------------
              elseif(Sclr(i).gt.-epsilonBT) then
                bub_info(i,1)  = xx(i,1)
                bub_info(i,2)  = xx(i,2)
                bub_info(i,3)  = xx(i,3)
                bub_info(i,4)  = elemvol_local(i)
                bub_info(i,5)  = elemvol_local(i)*denswght*rhogas
                bub_info(i,12) = u1(i)*elemvol_local(i)
                bub_info(i,13) = u2(i)*elemvol_local(i)
                bub_info(i,14) = u3(i)*elemvol_local(i)
                bub_info(i,15) = CurvInfo(i)
              else
!... collect the bubble information in details
                bub_info(i,1)  = xx(i,1)
                bub_info(i,2)  = xx(i,2)
                bub_info(i,3)  = xx(i,3)
                bub_info(i,4)  = elemvol_local(i)
                bub_info(i,5)  = elemvol_local(i)*denswght*rhogas
                bub_info(i,12) = u1(i)*elemvol_local(i)
                bub_info(i,13) = u2(i)*elemvol_local(i)
                bub_info(i,14) = u3(i)*elemvol_local(i)
                bub_info(i,16) = real(myrank)

               if(Sclr(i).ge.-2.0d0*epsilonBT)then
                bub_info(i,21) = 1.0d0
!               count element in the shell inside the bubble for boiling
               endif
                


              endif ! Sclr(i)
!-----------------------------------------------------------------------
!               if(Sclr(i).le.-1.0d0*epsilonBT
!     &                .and.Sclr(i).ge.-2.0d0*epsilonBT)then
!                bub_info(i,21) = 1.0d0
!               count element in the shell inside the bubble for boiling
!               endif
!-----------------------------------------------------------------------

           endif !Sclr(i) for the bubble region and liquid shell
!              if(bub_info(i,21).ne.0.0)write(*,*)'bub21', bub_info(i,21)
!              if(bub_info(i,20).ne.0.0)write(*,*)'bub20', bub_info(i,20)
!              if(bub_info(i,11).ne.0.0)write(*,*)'bub11', bub_info(i,11)

        enddo !npro

        end     !BubCollect ends
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
        subroutine BubPrintOut(vf)
!----------------------------------------------------------------------
!       Called in itrdrv.f
!       This subroutine will print out the bubble information into
!       external files which can be further processed
!
!       bub_cent        :x coord, y coord, z coord, radius, bubble ID
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"

        integer:: iNghost, Nbubdef
        real*8  vf

!        write(*,'(1x,A,F14.7,A)') 'Current void: ',vf*100.0d0,'%'
        Nbubdef = 0             !number of deformed bubbles
        iNghost = 0
        
        do i = 1, i_num_bubbles
        IF(avg_info(i,4).gt.0.0d0) THEN
           bub_cent(i,1:4) = avg_info(i,1:4)
           bub_cent(i,5)   = REAL(i) 

!       Find out the ghost bubbles and save them into bub_cent
           if(avg_info(i,1)-DomainSize(1).lt.GhostRatio*XLEN) then
              iNghost = iNghost + 1
              bub_cent(i_num_bubbles+iNghost,1)   = bub_cent(i,1)+XLEN
              bub_cent(i_num_bubbles+iNghost,2:3) = bub_cent(i,2:3)
              bub_cent(i_num_bubbles+iNghost,4)   = bub_cent(i,4)
              bub_cent(i_num_bubbles+iNghost,5)   = bub_cent(i,5)
           endif

           if(avg_info(i,2)-DomainSize(3).lt.GhostRatio*YLEN) then
              iNghost = iNghost + 1
              bub_cent(i_num_bubbles+iNghost,1)   = bub_cent(i,1)
              bub_cent(i_num_bubbles+iNghost,2)   = bub_cent(i,2)+YLEN
              bub_cent(i_num_bubbles+iNghost,3)   = bub_cent(i,3)
              bub_cent(i_num_bubbles+iNghost,4)   = bub_cent(i,4)
              bub_cent(i_num_bubbles+iNghost,5)   = bub_cent(i,5)
           endif

           if(avg_info(i,3)-DomainSize(5).lt.GhostRatio*ZLEN) then
              iNghost = iNghost + 1
              bub_cent(i_num_bubbles+iNghost,1:2) = bub_cent(i,1:2)
              bub_cent(i_num_bubbles+iNghost,3)   = bub_cent(i,3)+ZLEN
              bub_cent(i_num_bubbles+iNghost,4)   = bub_cent(i,4)
              bub_cent(i_num_bubbles+iNghost,5)   = bub_cent(i,5)
           endif

           if(DomainSize(2)-avg_info(i,1).lt.GhostRatio*XLEN) then
              iNghost = iNghost + 1
              bub_cent(i_num_bubbles+iNghost,1)   = bub_cent(i,1)-XLEN
              bub_cent(i_num_bubbles+iNghost,2:3) = bub_cent(i,2:3)
              bub_cent(i_num_bubbles+iNghost,4)   = bub_cent(i,4)
              bub_cent(i_num_bubbles+iNghost,5)   = bub_cent(i,5)
           endif

           if(DomainSize(4)-avg_info(i,2).lt.GhostRatio*YLEN) then
              iNghost = iNghost + 1
              bub_cent(i_num_bubbles+iNghost,1)   = bub_cent(i,1)
              bub_cent(i_num_bubbles+iNghost,2)   = bub_cent(i,2)-YLEN
              bub_cent(i_num_bubbles+iNghost,3)   = bub_cent(i,3)
              bub_cent(i_num_bubbles+iNghost,4)   = bub_cent(i,4)
              bub_cent(i_num_bubbles+iNghost,5)   = bub_cent(i,5)
           endif

           if(DomainSize(6)-avg_info(i,3).lt.GhostRatio*ZLEN) then
              iNghost = iNghost + 1
              bub_cent(i_num_bubbles+iNghost,1:2) = bub_cent(i,1:2)
              bub_cent(i_num_bubbles+iNghost,3)   = bub_cent(i,3)-ZLEN
              bub_cent(i_num_bubbles+iNghost,4)   = bub_cent(i,4)
              bub_cent(i_num_bubbles+iNghost,5)   = bub_cent(i,5)
           endif
        ENDIF !The bubble exists physically

!       the bubble is recognized to be deformed when the factor is less
!       than certain limit
           if(avg_info(i,4).gt.0.0d0 .and. avg_info(i,6).lt.0.8d0) 
     &        Nbubdef = Nbubdef + 1
        end do

        open(unit=740,file='../bubStat/alpha.dat',status="old",
     &       action="write", position="append", iostat=ierror)
        if(ierror /= 0) STOP "Error creating file 740"
        write(740,833) lstep+1, Nbubtot, Nbubdef, C_int_adjust, vf,
     &  time
        close(740)

 833  format(I6, 1x, I4, 1x, I4, 3x, ES14.7, 1x, ES14.7, 1x, ES14.7)

        end     !BubPrintOut ends
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
        subroutine CountAlphaLines(LineCount, inistep)
!----------------------------------------------------------------------
!       This subroutine will return the intial timestep in alpha.dat and
!       total number of lines of which
!
!----------------------------------------------------------------------
        use, intrinsic :: iso_fortran_env

        integer, intent(out) :: LineCount
        integer :: Read_Code
        character (len=200) :: line
        character (len=200) :: filename
        integer :: inistep, na1, na2
        real*8  :: dummya(3)

        open(unit=740,file='../bubStat/alpha.dat',status="old",
     &       action="read",iostat=ierror)
        read(740,823) inistep, na1, na2, dummya(:)
        close(740)

        open(unit=51,file='../bubStat/alpha.dat',
     &       status="old",access='sequential',
     &       form='formatted',action='read')

        LineCount = 0

        ReadLoop: do

        read (51, '(A)', iostat=Read_Code)  line

        if(Read_Code /= 0) then
           if(Read_Code == iostat_end ) then
              exit ReadLoop    ! end of file --> line count found
           else
              write ( *, '( / "read error: ", I0 )' )  Read_Code
              stop
           endif
        endif

        LineCount = LineCount + 1

!        write (*, '( I0, ": ''", A, "''" )' )  LineCount, trim (line)
        if(len_trim(line)==0)then
           write(*,'("The above is an empty or all blank line.")')
        endif

        end do ReadLoop

!        write (*, *) "found", LineCount, " lines"
        close(51)

 823  format(I6, 1x, I4, 1x, I4, 3x, ES14.7, 1x, ES14.7, 1x, ES14.7)

        end
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
        subroutine reCountBub()
!----------------------------------------------------------------------
!       Called in itrdrv.f
!       This subroutine will count how many bubble exist in current time
!       iteration and determine the number of ghost bubbles
!
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"

!       the bubble is recognized when its volume is larger than zero
        do i = 1, i_num_bubbles
           if(avg_info(i,4).gt.0.0d0) then
              Nbubtot = Nbubtot + 1
              if(avg_info(i,1)-DomainSize(1).lt.GhostRatio*XLEN)
     &              Nghost = Nghost + 1
              if(avg_info(i,2)-DomainSize(3).lt.GhostRatio*YLEN)
     &              Nghost = Nghost + 1
              if(avg_info(i,3)-DomainSize(5).lt.GhostRatio*ZLEN)
     &              Nghost = Nghost + 1

              if(DomainSize(2)-avg_info(i,1).lt.GhostRatio*XLEN)
     &              Nghost = Nghost + 1
              if(DomainSize(4)-avg_info(i,2).lt.GhostRatio*YLEN)
     &              Nghost = Nghost + 1
              if(DomainSize(6)-avg_info(i,3).lt.GhostRatio*ZLEN)
     &              Nghost = Nghost + 1
           endif !The bubble exists physically
        enddo


        end     !reCountBub ends
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
        subroutine ReColor(xltmp,bmtmp)
!----------------------------------------------------------------------
!       This subroutine is used to recolor the node coordinates and
!       current bubble centers 
!
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"

        real*8  bmtmp, xltmp(3)
        real*8  Mdist, MdistMin
        integer ib

        bmtmp = 0.0d0
        Mdist = 0.0d0
        MdistMin = 1.0E6

        do ib = 1,size(bub_cent,1) 
!       The bubble must exist physically and in the neighborhood       
           if(bub_cent(ib,4).gt.0.0d0 .and.  
     &        abs(xltmp(1)-bub_cent(ib,1)).le.NbrhdRatio*XLEN) then     

              Mdist = sqrt((xltmp(1)-bub_cent(ib,1))**2 +
     &                  (xltmp(2)-bub_cent(ib,2))**2 +
     &                  (xltmp(3)-bub_cent(ib,3))**2)-bub_cent(ib,4)

              if(MdistMin.gt.Mdist) then
                 MdistMin = Mdist
                 bmtmp = bub_cent(ib,5)         !update the marker
              endif
           endif
        enddo

        end     !ReColor ends
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
        subroutine breakupDetector(ibreakupFlag)
!----------------------------------------------------------------------
!       This subroutine is used to detect the suspicious bubble breakup
!       events 
!       This subroutine is only called by master core
!       breakupSeeder(ib,1)     : current bubble ID
!       breakupSeeder(ib,2)     : action flag
!       breakupSeeder(ib,3)     : one parallel rank that where the
!                                 bubble has elements on
!       breakupSeeder(ib,4)     : total number of bubble point when
!                                 seeding
!       breakupSeeder(ib,5)     : volume of region occupied by auxiliary
!                                 ID
!       breakupSeeder(ib,6)     : auxiliary bubble ID if flag is not 0
!
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"

        integer ib, iPossibleBreakup
        integer iactionFlag, iauxBubID, MaxAuxID
        integer ibreakupFlag, iseeding


        ibreakupFlag    = 0
        iPossibleBreakup= 0
        MaxAuxID        = 0

        do ib = 1,i_num_bubbles
           if(MaxAuxID.lt.int(breakupSeeder(ib,6)))
     &        MaxAuxID =  int(breakupSeeder(ib,6))
        enddo 

!        write(*,*) 'i_num_bubbles in detector =', i_num_bubbles

        do ib = 1,i_num_bubbles
!           write(*,*) 'breakupSeeder in detector 1= ',
!     &                 breakupSeeder(ib,:)
           iseeding     = 0
           breakupSeeder(ib,1) = real(ib)
           iactionFlag  = int(breakupSeeder(ib,2))
           if(iactionFlag.eq.0.and.avg_info(ib,4).gt.0.0d0) then
           iauxBubID    = int(breakupSeeder(ib,1))
           if(mod(lstep,100).eq.0) iseeding     = 1
!       more sophisticated condition can be developed to determine when
!       setting iseeding to be 1. 
           if(iseeding.eq.1) then
!           if(lstep+1.eq.17403) then
!              write(*,*) 'ID seeding is switched on manually'
              iPossibleBreakup = iPossibleBreakup + 1
              iactionFlag      = 1
              iauxBubID        = 3*MaxAuxID + iPossibleBreakup
           endif

           breakupSeeder(ib,2) = real(iactionFlag)
           breakupSeeder(ib,3) = avg_info(ib,14)
           breakupSeeder(ib,6) = real(iauxBubID)

           endif !iactionFlag = 0
!           write(*,*) 'breakupSeeder in detector = ',
!     &                 breakupSeeder(ib,:)

        enddo  !ib


!       The ibreakupFlag is used to determine if we should call breakup
!       confirmer, basically if there is one non-zero breakup index, the
!       flag will be 1. 
        do ib = 1, i_num_bubbles
           if(breakupSeeder(ib,2).gt.0.0d0) then
              ibreakupFlag = 1
              exit
           endif
        enddo

        if(ibreakupFlag.eq.0) then
             write(*,*) 'There is no suspicious break-up event!'
        else
             write(*,*) 'Suspicious break-up event is detected!'
        endif

        end !breakupDetector ends
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
        subroutine breakupConfirmer(y, banma)
!----------------------------------------------------------------------
!       This subroutine is used to double check the suspicious bubble
!       breakup events
!       This subroutine is called by all the parallel cores
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        real*8    y(nshg,ndof)
        real*8    Ratio
        dimension banma(nshg,1)
        integer   i, ib, i_mrk
        integer   iactionFlag, iauxBubID, irank

        do ib = 1,i_num_bubbles
           iactionFlag  = int(breakupSeeder(ib,2))
           irank        = int(breakupSeeder(ib,3))
           iauxBubID    = int(breakupSeeder(ib,6))
           Ratio        = 0.0d0

           if(iactionFlag.eq.1) then
              iactionFlag = 2

              if(myrank.ne.master .and. myrank.eq.irank) then 
                 do i = 1, nshg
                 i_mrk = int(banma(i,1)) 
                 if(y(i,6).lt.0.0d0 .and. i_mrk.eq.ib) then
                 banma(i,1) = real(iauxBubID)
                 if(iactionFlag.eq.2) exit

                 endif
                 enddo
              endif

           elseif(iactionFlag.ge.2.and.iactionFlag.le.14) then
                  iactionFlag = iactionFlag + 1

           elseif(iactionFlag.eq.15) then
!       Confirm the exitance of new bubble and update the total number
!       of bubbles 
               iactionFlag = 0
               if(breakupSeeder(ib,4).gt.0.0d0) 
     &            Ratio = breakupSeeder(ib,5)/breakupSeeder(ib,4)
!               if(myrank.eq.master)write(*,*)'Ratio =',ib, Ratio
               if(Ratio.lt.0.1d0 .or. Ratio.gt.0.9d0 ) then  
!       case 1: fake trigger and no breakup is happening
                  do i = 1, nshg
                     if(int(banma(i,1)).eq.iauxBubID) 
     &               banma(i,1) = real(ib)
                  enddo 
               else
!       case 2: breakup event is captured and recorded
                  i_num_bubbles = i_num_bubbles + 1
                  do i = 1, nshg
                     if(int(banma(i,1)).eq.iauxBubID)
     &               banma(i,1) = real(i_num_bubbles)  
                  enddo
                  if(myrank.eq.master) 
     &            write(*,'(A,I5,A)')' New Bubble',i_num_bubbles,
     &                               ' is recognized!'
               endif 

               breakupSeeder(ib,6) = breakupSeeder(ib,1)
           endif

           breakupSeeder(ib,2)  = real(iactionFlag)
        enddo

!        if(myrank.eq.master) then
!        write(*,*) 'breakupSeeder in confirmor = ',
!     &                  breakupSeeder(:,:)
!        endif

        end !breakupConfirmer ends
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
        subroutine banmaCorrector(banma, ibreakupFlag, yold)
!----------------------------------------------------------------------
!       The subroutine will correct the marker field in restart files
!       when breakup algorithm is active
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        real*8    yold(nshg,ndof)
        dimension banma(nshg,1), bmALL(nshg,1)
        integer   i, i_mrk, ib
        integer   ibreakupFlag
        bmALL = banma

        if(ibreakupFlag.eq.0) then
           yold(:,7) = bmALL(:,1)
        else
!       Correct the banma field 
           do i=1,nshg
              i_mrk = int(bmALL(i,1))
              if(i_mrk.gt.i_num_bubbles) then
                 do ib=1,i_num_bubbles
                 if(i_mrk.eq.int(breakupSeeder(ib,6))) then
                    bmALL(i,1) = breakupSeeder(ib,1)
                    exit
                 endif
                 enddo
              endif
           enddo
           yold(:,7) = bmALL(:,1)
        endif

        end !banmaCorrector ends
c======================================================================
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c======================================================================
