!----------------------------------------------------------------------
!       This file contains several subroutines working for lift/drag
!       control mechanism developed by Jun Fang and Matt Thomas. And the
!       codes are further summarized and rewritten by Jun in a neater
!       manner. 
!       This file also contains the subroutines used in Jun Fang bubble
!       coalescence control algorithm. 
!       
!       Jun Fang, 2015
!----------------------------------------------------------------------

c======================================================================
c======================================================================
c======================================================================
        subroutine coalescenceDetection()
!----------------------------------------------------------------------
!       This subroutine is used to detect the number of coalescence
!       events in multi-bubble simulations and determine the centers of
!       control force application regions to prevent coalescence. 
!
!       ncoalEvent      : number of coalescence events
!       bubSurfDist     : the distance between two bubble interface
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        integer icoalEvent
        integer ib, jb, i
        real*8  bubSurfDist

        ncoalEvent      = 0
        do ib = 1, i_num_bubbles+Nghost-1
           if(bub_cent(ib,4).gt.0.0d0) then
              do jb = ib+1, i_num_bubbles+Nghost
                 if(bub_cent(jb,4).gt.0.0d0) then
                    bubSurfDist = sqrt(
     &                        (bub_cent(ib,1)-bub_cent(jb,1))**2.0d0 
     &                       +(bub_cent(ib,2)-bub_cent(jb,2))**2.0d0
     &                       +(bub_cent(ib,3)-bub_cent(jb,3))**2.0d0
     &                       ) - bub_cent(ib,4) - bub_cent(jb,4)
                        if(bubSurfDist.le.coalcon_dist*epsilonBT)
     &                 ncoalEvent = ncoalEvent + 1
                endif
              enddo
           endif
        enddo
        write(*,'(1x,A,I5)') 'ncoalEvent = ', ncoalEvent
        if(ncoalEvent.gt.coalest)
     &  write(*,*) 'Warning: estimated coalescn number is exceeded!'

        coalCenter(:,:) = 0.0d0

        if(ncoalEvent.ge.1) then !collecting coalescence centers' info
        icoalEvent      = 1
        do ib = 1, i_num_bubbles+Nghost-1
           if(bub_cent(ib,4).gt.0.0d0) then
              do jb = ib+1, i_num_bubbles+Nghost
                 if(bub_cent(jb,4).gt.0.0d0) then
                    bubSurfDist = sqrt(
     &                        (bub_cent(ib,1)-bub_cent(jb,1))**2.0d0
     &                       +(bub_cent(ib,2)-bub_cent(jb,2))**2.0d0
     &                       +(bub_cent(ib,3)-bub_cent(jb,3))**2.0d0
     &                       ) - bub_cent(ib,4) - bub_cent(jb,4)
                   if(bubSurfDist.le.coalcon_dist*epsilonBT.and.
     &                icoalEvent.le.coalest) then
                       coalCenter(icoalEvent,:) =
     &                 (bub_cent(ib,:)+bub_cent(jb,:))/2.0d0
                       icoalEvent = icoalEvent + 1
                    endif
                 endif
              enddo
           endif
        enddo
        endif !at least one potential coalescence event is detected. 
        if(icoalcon_verbo.eq.1) then
        do i = 1,ncoalEvent
           write(*,'(1x,A,I4,A,3ES14.4)')'Suspicious coalescence ',i,
     &             ' is found at ', coalCenter(i,:)

        enddo
        endif

        end     !coalescenceDetection ends
c======================================================================
c======================================================================
c======================================================================
        subroutine coalCtrlApp(xxtmp, Wetmp)
!----------------------------------------------------------------------
!       In this subroutine, the coalescence control force will be
!       applied in the form of surface tension force. 
!
!----------------------------------------------------------------------
        use bub_track   ! access to bubble information array
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        real*8   xxtmp(nsd), Wetmp
        real*8   ctrlRange
        integer  ib

        do ib = 1, ncoalEvent
           ctrlRange = sqrt(
     &                 (coalCenter(ib,1)-xxtmp(1))**2.0d0
     &                +(coalCenter(ib,2)-xxtmp(2))**2.0d0
     &                +(coalCenter(ib,3)-xxtmp(3))**2.0d0
     &                 )

           if(ctrlRange.le.5.0d0*epsilonBT) then
              Wetmp = CoalInvSigma
              exit
           endif
        enddo

        end     !coalCtrlApp ends
c======================================================================
c======================================================================
c======================================================================

!----------------------------------------------------------------------
! MAGNUS
! From this point onwards I will implement the subroutines necessary for
! the bubble controller to work

      subroutine CFrestar(ycf_old_log, numoldyhistind,
     &                      ixcf,        iycf,       izcf,
     &                      sumxcf_o,    sumycf_o,   sumzcf_o,
     &                      oldcf_dt,    oldcf_dtlset)
!
!-----------------------------------------------------------------------
!
! This subroutine is used to reload the information of control forces 
! and open necessary files when starting cases. The prototype is 
! Matt's stuff in major scripts, such as itrdrv.f, elmgmr.f and so on. 
!
! Jun Fang,  Matt. A. Thomas                             Fall, 2014 
!-----------------------------------------------------------------------
      use cf_arrays  ! m
      USE, INTRINSIC :: ISO_C_BINDING !for calling C++ routines
      use bub_track ! magnus, mpid
        
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer ierror, nxyzcflines
      integer i, j, k
      integer ixcf,   iycf,  izcf
      integer ioldyhistst,  ioldyhisten
      integer numoldyhistind
      integer indyhist
      integer bub_id

      real*8  sumxcf_o(i_num_bubbles), sumycf_o(i_num_bubbles), sumzcf_o(i_num_bubbles)
      real*8  xcf_old_dummy
      real*8, allocatable :: ycf_old_dummy(:,:)     !mpid, extend dimensions
      real*8  zcf_old_dummy
      real*8  cf_restart_array(32,i_num_bubbles) !mpid, extend dimensions
      real*8  oldcf_dt, oldcf_dtlset 

      logical ycf_old_log
      logical ResCFexist      !to continue matts cf after a break in phasta
      character*50:: lstepc
      character*2 :: id_char ! mpid, character for the bubble id
      logical dir_lf
      logical dir_con             ! mpid



!----------------------------------------------------------------------
! print the number of files required for mpid - magnus
      if (myrank.eq.master) then
        write(*,*) 'number of bubbles to control', i_num_bubbles
      end if
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!c.... Open files for matt thomas control force
      if(myrank.eq.master) write(*,*)
      if(myrank.eq.master) write(*,*)'iCForz ', iCForz
      if(myrank.eq.master) write(*,*)'iCForz_where ',iCForz_where
      if(myrank.eq.master) write(*,*)

!----------------------------------------------------------------------
!       Check the existence of folder lift-results, create a new one if
!       not
!----------------------------------------------------------------------

! ---------------------------------------------------------------------
! check controller directory exists, for mpid directory structure
      if(myrank.eq.master) then
              inquire(directory ='../controller', exist=dir_con)
              if (dir_con.eqv..False.) then
                     write(*,*) 'controller directory does not exist.', 
     &                             ' Create a new one!'
                     call system('mkdir ../controller')
              else   !controller dir alredy exists
                     write(*,*) 'controller directory already exists'
              end if !controller dir exists
       end if ! rank .eq. master
! ----------------------------------------------------------------------

      do bub_id=1, i_num_bubbles !for mpid a create directory for each bubble
!       write(*,*) 'bub_id is', bub_id
        write(id_char, '(I0)') bub_id   !convert bub_id to character
        if(myrank.eq.master)then
            inquire(directory='../controller/'//trim(id_char)//'_lift-results',
     &          exist=dir_lf)
            if(dir_lf.eqv..False.) then
            write(*,*) 'lift-results not found for bubble ', id_char,
     &                 '... Create a new one!'
            write(*,*)
            call system('mkdir ../controller/'//trim(id_char)//'_lift-results')
            else
            write(*,*) 'lift-results already exists for bubble', id_char
            write(*,*)
            end if

            inquire(directory='../controller/'//trim(id_char)//'_lift-results/restarts',exist=dir_lf)
            if(dir_lf.eqv..False.) then
            call system('mkdir ../controller/'//trim(id_char)//'_lift-results/restarts')
            endif

        end if
      end do    ! bubble count

      do bub_id=1, i_num_bubbles     !mpid bubble loop
        write(id_char, '(I0)') bub_id   !convert id int to char

      if(myrank.eq.master)then
         if(     lstep.ge.0 .and. lstep.le.9)then
           write(lstepc,'(i1.1)')lstep
         else if(lstep.ge.10 .and. lstep.le.99)then
           write(lstepc,'(i2.2)')lstep
         else if(lstep.ge.100 .and. lstep.le.999)then
           write(lstepc,'(i3.3)')lstep
         else if(lstep.ge.1000 .and. lstep.le.9999)then
           write(lstepc,'(i4.4)')lstep
         else if(lstep.ge.10000 .and. lstep.le.99999)then
           write(lstepc,'(i5.5)')lstep
         else if(lstep.ge.100000 .and. lstep.le.999999)then
           write(lstepc,'(i6.6)')lstep
         else if(lstep.ge.1000000 .and. lstep.le.9999999)then
           write(lstepc,'(i7.7)')lstep
         else if(lstep.ge.10000000 .and. lstep.le.99999999)then
           write(lstepc,'(i8.8)')lstep
         else if(lstep.ge.100000000 .and. lstep.le.999999999)then
           write(lstepc,'(i9.9)')lstep
         else if(lstep.ge.1000000000 .and. lstep.le.9999999999)then
           write(lstepc,'(i10.10)')lstep
         end if

         write(*,*)'Looking for restart_control_force'//trim(lstepc)//
     &            '.dat'

         inquire(file='../controller/'//trim(id_char)//'_lift-results/restarts/'//
     &             'restart_control_force'//
     &             trim(lstepc)//'.dat',
     &             exist=ResCFexist)

         if(ResCFexist .eqv. .true.) then       !ResCFexist is true

            i_res_cf = 1
            write(*,*)'for bubble ', id_char            
            write(*,*)'restart_control_force'//trim(lstepc)//
     &                  '.datexists'
            write(*,*)'Continuing control force based on ',
     &                  'restart_control_force'//trim(lstepc)//'.dat'
            write(*,*)'Pre-processing lift-results files to ensure ',
     &                  'lift result file continuity'


            call preprocess_lift_files(lstep)
       if(myrank.eq.master) write(*,*) 'preprocessed_lift_files!' !mpid debug


      OPEN(unit=7001+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/average_distance.dat',
     &  status="old",action="write",position="append",iostat=ierror)
      IF(ierror /= 0) STOP "Error opening file 741"

      OPEN(unit=7002+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/average_velocity.dat',
     &  status="old",action="write",position="append",iostat=ierror)
      IF(ierror /= 0) STOP "Error opening file 742"
     
      OPEN(unit=7003+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/total_control_force.dat',
     &  status="old",action="write",position="append",iostat=ierror)
      IF(ierror /= 0) STOP "Error opening file 743"
     
      OPEN(unit=7004+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/historical_average.dat',
     &  status="old",action="write",position="append",iostat=ierror)
      IF(ierror /= 0) STOP "Error opening file 744"
     
      OPEN(unit=7005+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/pos_vel_difference.dat',
     &  status="old",action="write",position="append",iostat=ierror)
      IF(ierror /= 0) STOP "Error opening file 745"
     
      OPEN(unit=7006+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/average_control_force.dat',
     &  status="old",action="write",position="append",iostat=ierror)
      IF(ierror /= 0) STOP "Error opening file 746"
     
      OPEN(unit=7007+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/newtons_control_force.dat',
     &  status="old",action="write",position="append",iostat=ierror)
      IF(ierror /= 0) STOP "Error opening file 747"

      !mpid, below change 845 --> 8001+bub_id*10
      OPEN(unit=8001+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/restarts/'//
     &                    'restart_control_force'//
     &                    trim(lstepc)//'.dat',
     &  status="old",action="read",iostat=ierror)
      IF(ierror /= 0) STOP "Error opening file 845"

         else if(ResCFexist .eqv. .false.) then !ResCFexist is false

            i_res_cf = 0

            if(     lstep.ge.0 .and. lstep.le.9)then
              write(lstepc,'(i1.1)')lstep
            else if(lstep.ge.10 .and. lstep.le.99)then
              write(lstepc,'(i2.2)')lstep
            else if(lstep.ge.100 .and. lstep.le.999)then
              write(lstepc,'(i3.3)')lstep
            else if(lstep.ge.1000 .and. lstep.le.9999)then
              write(lstepc,'(i4.4)')lstep
            else if(lstep.ge.10000 .and. lstep.le.99999)then
              write(lstepc,'(i5.5)')lstep
            else if(lstep.ge.100000 .and. lstep.le.999999)then
              write(lstepc,'(i6.6)')lstep
            else if(lstep.ge.1000000 .and. lstep.le.9999999)then
              write(lstepc,'(i7.7)')lstep
            else if(lstep.ge.10000000 .and. lstep.le.99999999)then
              write(lstepc,'(i8.8)')lstep
            else if(lstep.ge.100000000 .and. lstep.le.999999999)then
              write(lstepc,'(i9.9)')lstep
            else if(lstep.ge.1000000000 .and. lstep.le.9999999999)then
              write(lstepc,'(i10.10)')lstep
            end if

      OPEN(unit=8002+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/restarts/'//
     &             'restart_control_force'//trim(lstepc)//'.dat',
     &             status="new",action="write",iostat=ierror)
      IF(ierror /= 0) STOP "file 846 creation error"

      if(myrank.eq.master)
     &  write(*,*)'restart_control_force.dat does not exist'
      
      if(myrank.eq.master)
     &  write(*,*)'Creating restart_control_force'//trim(lstepc)//'.dat'
    
      if(myrank.eq.master)
     &  write(*,*)'Before stopping phasta, ',
     &  'wait until a restart cf file is populated ',
     &  '(if you desire to restart ',
     &  'the simulation from that timestep.) '
      
      if(myrank.eq.master)
     &  write(*,*)'If this cf restart file is not ',
     &  'populated before exiting phasta'
      if(myrank.eq.master)write(*,*)'delete it and re-run'
      CLOSE(8002+bub_id*10)

!      if(myrank.eq.master) write(*,*) 'first open statement' !debug
      OPEN(unit=7001+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/average_distance.dat',
     &   status="new",action="write",iostat=ierror)
      IF(ierror /= 0) STOP "file 741 creation error"
!      if(myrank.eq.master) write(*,*) 'first open statement' !debug
     
      OPEN(unit=7002+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/average_velocity.dat',
     &   status="new",action="write",iostat=ierror)
      IF(ierror /= 0) STOP "file 742 creation error"
     
      OPEN(unit=7003+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/total_control_force.dat',
     &   status="new",action="write",iostat=ierror)
      IF(ierror /= 0) STOP "file 743 creation error"
      
      OPEN(unit=7004+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/historical_average.dat',
     &   status="new",action="write",iostat=ierror)
      IF(ierror /= 0) STOP "file 744 creation error"
      
      OPEN(unit=7005+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/pos_vel_difference.dat',
     &   status="new",action="write",iostat=ierror)
      IF(ierror /= 0) STOP "file 745 creation error"
      
      OPEN(unit=7006+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/average_control_force.dat',
     &   status="new",action="write",iostat=ierror)
      IF(ierror /= 0) STOP "file 746 creation error"
     
      OPEN(unit=7007+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/newtons_control_force.dat',
     &   status="new",action="write",iostat=ierror)
      IF(ierror /= 0) STOP "file 747 creation error"
     
      OPEN(unit=7008+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/xyzcf.dat',
     &   status="new",action="write",iostat=ierror)
      IF(ierror /= 0) STOP "file 748 creation error"

        end if                                          !resCFexist
        end if                                          !myrank
      !if(myrank.eq.master) write(*,*) 'about to call MPI_BARRIER' !debug
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !if(myrank.eq.master) write(*,*) 'called MPI_BARRIER' !debug
       !if(myrank.eq.master) write(*,*) 'about to call MPI_Bcast' !debug
      if (numpe > 1) call MPI_Bcast(i_res_cf,1,MPI_INTEGER,master,
     &                                MPI_COMM_WORLD,ierr)
      !if(myrank.eq.master) write(*,*) 'called MPI_Bcast' !debug
!--------------------------------------------------------------------
!...initialize stuff for matts control algorithm, and read the history
!for control force capability
      
      if(i_res_cf.eq.1)then
         if(myrank.eq.master)then
!              if(myrank.eq.master) write(*,*) 'about to read 845'
            read(8001+bub_id*10,733)avgycf(bub_id),
     &      avgycforceold(bub_id),dy_new(bub_id),avgyvelold(bub_id),
     &      ddyvel(bub_id),avgydistold(bub_id),ydistideal(bub_id),sumycf_o(id),
     &      avgxcf(bub_id),avgxcforceold(bub_id),dx_new(bub_id),avgxvelold(bub_id),
     &      ddxvel(bub_id),avgxdistold(bub_id),xdistideal(bub_id),sumxcf_o(bub_id),
     &      avgzcf(bub_id),avgzcforceold(bub_id),dz_new(bub_id),avgzvelold(bub_id),
     &      ddzvel(bub_id),avgzdistold(bub_id),zdistideal(bub_id),sumzcf_o(bub_id),
     &      epsilon_lsd,oldcf_dt,oldcf_dtlset,CFLfl_max,
     &      CFLbuint_max,CFLls_max,vf_now,C_int_adjust
            close(8001+bub_id*10)
!            if(myrank.eq.master) write(*,*) 'closed 845'

            cf_restart_array(1,bub_id)  = avgycf(bub_id)
            cf_restart_array(2,bub_id)  = avgycforceold(bub_id)
            cf_restart_array(3,bub_id)  = dy_new(bub_id)
            cf_restart_array(4,bub_id)  = avgyvelold(bub_id)
            cf_restart_array(5,bub_id)  = ddyvel(bub_id)
            cf_restart_array(6,bub_id)  = avgydistold(bub_id)
            cf_restart_array(7,bub_id)  = ydistideal(bub_id)
            cf_restart_array(8,bub_id)  = sumycf_o(bub_id)
            cf_restart_array(9,bub_id)  = avgxcf(bub_id)
            cf_restart_array(10,bub_id) = avgxcforceold(bub_id)
            cf_restart_array(11,bub_id) = dx_new(bub_id)
            cf_restart_array(12,bub_id) = avgxvelold(bub_id)
            cf_restart_array(13,bub_id) = ddxvel(bub_id)
            cf_restart_array(14,bub_id) = avgxdistold(bub_id)
            cf_restart_array(15,bub_id) = xdistideal(bub_id)
            cf_restart_array(16,bub_id) = sumxcf_o(bub_id)
            cf_restart_array(17,bub_id) = avgzcf(bub_id)
            cf_restart_array(18,bub_id) = avgzcforceold(bub_id)
            cf_restart_array(19,bub_id) = dz_new(bub_id)
            cf_restart_array(20,bub_id) = avgzvelold(bub_id)
            cf_restart_array(21,bub_id) = ddzvel(bub_id)
            cf_restart_array(22,bub_id) = avgzdistold(bub_id)
            cf_restart_array(23,bub_id) = zdistideal(bub_id)
            cf_restart_array(24,bub_id) = sumzcf_o(bub_id)
            cf_restart_array(25,bub_id) = epsilon_lsd
            cf_restart_array(26,bub_id) = oldcf_dt
            cf_restart_array(27,bub_id) = oldcf_dtlset
            cf_restart_array(28,bub_id) = CFLfl_max
            cf_restart_array(29,bub_id) = CFLbuint_max
            cf_restart_array(30,bub_id) = CFLls_max
            cf_restart_array(31,bub_id) = vf_now
            cf_restart_array(32,bub_id) = C_int_adjust
         end if !master

         if (numpe > 1) then
            call MPI_Bcast(cf_restart_array,32,MPI_DOUBLE_PRECISION,
     &                       master,MPI_COMM_WORLD,ierr)
         end if
!...the following commands make sure that all the processors get the
!same initial values of control force variables
         
         if(myrank.ne.master) then
!           write(*,*) 'about to allocate pid variables...'
            avgycf(bub_id)            = cf_restart_array(1,bub_id)
            avgycforceold(bub_id)     = cf_restart_array(2,bub_id)
            dy_new(bub_id)            = cf_restart_array(3,bub_id)
            avgyvelold(bub_id)        = cf_restart_array(4,bub_id)
            ddyvel(bub_id)            = cf_restart_array(5,bub_id)
            avgydistold(bub_id)       = cf_restart_array(6,bub_id)
            ydistideal(bub_id)        = cf_restart_array(7,bub_id)
            sumycf_o(bub_id)          = cf_restart_array(8,bub_id)
            avgxcf(bub_id)            = cf_restart_array(9,bub_id)
            avgxcforceold(bub_id)     = cf_restart_array(10,bub_id)
            dx_new(bub_id)            = cf_restart_array(11,bub_id)
            avgxvelold(bub_id)        = cf_restart_array(12,bub_id)
            ddxvel(bub_id)            = cf_restart_array(13,bub_id)
            avgxdistold(bub_id)       = cf_restart_array(14,bub_id)
            xdistideal(bub_id)        = cf_restart_array(15,bub_id)
            sumxcf_o(bub_id)          = cf_restart_array(16,bub_id)
            avgzcf(bub_id)            = cf_restart_array(17,bub_id)
            avgzcforceold(bub_id)     = cf_restart_array(18,bub_id)
            dz_new(bub_id)            = cf_restart_array(19,bub_id)
            avgzvelold(bub_id)        = cf_restart_array(20,bub_id)
            ddzvel(bub_id)            = cf_restart_array(21,bub_id)
            avgzdistold(bub_id)       = cf_restart_array(22,bub_id)
            zdistideal(bub_id)        = cf_restart_array(23,bub_id)
            sumzcf_o(bub_id)          = cf_restart_array(24,bub_id)
            epsilon_lsd       = cf_restart_array(25,bub_id)
            oldcf_dt          = cf_restart_array(26,bub_id)
            oldcf_dtlset      = cf_restart_array(27,bub_id)
            CFLfl_max         = cf_restart_array(28,bub_id)
            CFLbuint_max      = cf_restart_array(29,bub_id)
            CFLls_max         = cf_restart_array(30,bub_id)
            vf_now           = cf_restart_array(31,bub_id)
            C_int_adjust      = cf_restart_array(32,bub_id)
!            write(*,*) 'allocated mpid variables'
         end if

         if(myrank.eq.master)then
          write(*,*)
          write(*,*) 'For bubble ', id_char
          write(*,*)'Using eps_lsd from cf restart file, ',
     &       'epsilon_lsd = ',  epsilon_lsd
          write(*,*)'Using dt from last cf restart file, ',
     &       'oldcf_dt = ',     oldcf_dt
          write(*,*)'Using dtlset from cf restart file, ',
     &       'oldcf_dtlset = ', oldcf_dtlset
          write(*,*)'Using CFLfl_max from cf restart file, ',
     &       'CFLfl_max = ',    CFLfl_max
          write(*,*)'Using CFLbuint_max from cf restart file, ',
     &       'CFLbuint_max = ', CFLbuint_max
          write(*,*)'Using CFLls_max from cf restart file, ',
     &       'CFLls_max = ',    CFLls_max
          write(*,*)'Using vf_now from cf restart file, ',
     &       'vf_now = ',   vf_now
          write(*,*)'Using C_int_adjust from cf restart file, ',
     &       'C_int_adjust = ', C_int_adjust
          write(*,*)
         end if

      else !i_res_cf is 0 for the following algorithm
         avgycf(bub_id)               = 0.0d0
         avgycforceold(bub_id)        = 0.0d0
         dy_new(bub_id)               = 0.0d0
         avgyvelold(bub_id)           = 0.0d0
         ddyvel(bub_id)               = 0.0d0
         avgydistold(bub_id)          = 0.0d0
         ydistideal(bub_id)           = 0.0d0
         sumycf_o(bub_id)             = 0.0d0
         avgxcf(bub_id)               = 0.0d0
         avgxcforceold(bub_id)        = 0.0d0
         dx_new(bub_id)               = 0.0d0
         avgxvelold(bub_id)           = 0.0d0
         ddxvel(bub_id)               = 0.0d0
         avgxdistold(bub_id)          = 0.0d0
         xdistideal(bub_id)           = 0.0d0
         sumxcf_o(bub_id)             = 0.0d0
         avgzcf(bub_id)               = 0.0d0
         avgzcforceold(bub_id)        = 0.0d0
         dz_new(bub_id)               = 0.0d0
         avgzvelold(bub_id)           = 0.0d0
         ddzvel(bub_id)               = 0.0d0
         avgzdistold(bub_id)          = 0.0d0
         zdistideal(bub_id)           = 0.0d0
         sumzcf_o(bub_id)             = 0.0d0
      end if !i_res_cf ends here

!...allocate and initialize for matt's control force

      if(i_res_cf.eq.1)then

         ycf_old_log = .true.

!...need to find out how many lines in xyzcf.dat file
         if(myrank.eq.master)then
          call lcounter(nxyzcflines,'../controller/'//trim(id_char)//'_lift-results/xyzcf.dat')
          write(*,*)'The number of lines in xyzcf.dat is ',nxyzcflines

          OPEN(unit=7008+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/xyzcf.dat',
     &      status="old",action="read",iostat=ierror)
          IF(ierror /= 0) STOP "Error opening file 748, earlier"

         end if

      if (numpe > 1) then
         call MPI_Bcast(nxyzcflines,1,MPI_INTEGER,
     &                    master,MPI_COMM_WORLD,ierr)
      endif
!The variable nxyzcflines has been successfully broadcast to different
!processors, which was tested at 9/8/2014 by Jun. 
!        if(myrank.eq.master) write(*,*) 'nxyzcflines is broadcast!'
!        if(myrank.ne.master) write(*,*)'rank',myrank,
!     &                       'value',nxyzcflines

!...get x,y,zcf from last simulation before the restart right now since we 
!   only worry about yhist term averaging over certain number of timesteps 
!   (i.e. x and z term are always averaged over total simulation timestep, 
!   therefore x,zcf_old aren't even used) then we don't allocate x,zcf_old 
!   arrays and we just read into dummy variables for them.

!...The array ycf_old has to be declared in itrdrv.f
!           allocate (ycf_old(numoldyhistind))
!        if(myrank.ne.master) write(*,*) 'rank: ',myrank,
!     &                                  'Before reading the cf records'
         if(myrank.eq.master)then
!...only store values of xyzcf from last sim that will be used
!   in this simulation
          ioldyhistst = (lstep+1) - numts_histyavg + 1
          ioldyhisten = (lstep+1) - 1
          write(*,*) 'ioldyhistst =', ioldyhistst
          write(*,*) 'ioldyhisten =', ioldyhisten

          indyhist=0
          do i=1,nxyzcflines
           if((i.ge.ioldyhistst).and.(i.le.ioldyhisten))then
            indyhist = indyhist + 1
            read(7008+bub_id*10,738)xcf_old_dummy,ycf_old(indyhist,bub_id),zcf_old_dummy
           else
            read(7008+bub_id*10,738)xcf_old_dummy,
     &        ycf_old_dummy(indyhist,bub_id) ,zcf_old_dummy
           endif

          enddo !i=1,nxyzcflines
          close(7008+bub_id*10)
         endif !(myrank.eq.master)

!        if(myrank.eq.master) write(*,*) 'ycf_old is to be broadcast!'
         if (numpe > 1) then
          call MPI_BCAST(ycf_old,numoldyhistind,MPI_DOUBLE_PRECISION,
     &                       master,MPI_COMM_WORLD,ierr)
         end if
!        if(myrank.eq.master) write(*,*) 'ycf_old is broadcast!'

         if(myrank.eq.master)then

          OPEN(unit=7009+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/xyzcf_old.dat',
     &      status="replace",action="write",iostat=ierror)
          IF(ierror /= 0) STOP "Error opening file 749"

          do i=1,numoldyhistind
           write(7009+bub_id*10,738)xcf_old_dummy,ycf_old(i,bub_id),zcf_old_dummy
          end do
          close(7009+bub_id*10)

          OPEN(unit=7008+bub_id*10,file='../controller/'//trim(id_char)//'_lift-results/xyzcf.dat',
     &      status="old",action="write",position="append",iostat=ierror)
          IF(ierror /= 0) STOP "Error opening file 748, later"

         end if !myrank

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       else !i_res_cf

          avgxcf(bub_id) = 0.0d0
          avgycf(bub_id) = 0.0d0
          avgzcf(bub_id) = 0.0d0

          dx_new(bub_id) = 0.0d0
          dy_new(bub_id) = 0.0d0
          dz_new(bub_id) = 0.0d0

      end if !i_res_cf

         ixcf=1
         iycf=1
         izcf=1
!...got to calculate where drag sign flips for matt's control force
       y_drag_flip = 0.0d0  - vel_centre / shear_rate
       if(myrank.eq.master)then
        write(*,*)'Shear rate is ',shear_rate
        write(*,*)'Centerline velocity is ', vel_centre
        write(*,*)'drag force flips signs at ',y_drag_flip
        write(*,*)
      end if

 733  format(32(E22.15, 3x))
 736  format(E22.15, 2x, E22.15, 2x, E22.15, 2x, E22.15)
 737  format(I6, 2x, E22.15, 2x, E22.15, 2x, E22.15)
 738  format(E22.15, 2x, E22.15, 2x, E22.15)
 739  format(E22.15, 2x, E22.15, 2x, E22.15, 2x,
     1           E22.15, 2x, E22.15, 2x, E22.15)

        end do                                          !mpid, bubble loop

        end !CFrestar initialization
c=====================================================================
c=====================================================================
c=====================================================================
      subroutine CFcalculator()
c!
c!----------------------------------------------------------------------
c!Called in elmgmr.f
c!This subroutien is used to calculate the control forces based on 
c!necessary variables, the corresponding constants can be modified in
c!solver.inp.  
c!
c! Jun, Fall 2014
c!----------------------------------------------------------------------
      use bub_track
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

!       real*8, allocatable :: xdistancesum(:)

!       allocate(xdistancesum(i_num_bubbles))

!       initialize values for control force before loop over
!       element-blocks
      xdistancesum    = 0.0d0
      ydistancesum    = 0.0d0
      zdistancesum    = 0.0d0
      xvelsum         = 0.0d0
      yvelsum         = 0.0d0
      zvelsum         = 0.0d0
      xcforcesum      = 0.0d0
      ycforcesum      = 0.0d0
      zcforcesum      = 0.0d0
      xforcenewtsum   = 0.0d0
      yforcenewtsum   = 0.0d0
      zforcenewtsum   = 0.0d0
      velwghtsum      = 0.0d0
      bubvolsum       = 0.0d0
      denssum         = 0.0d0
      nzinBsum        = 0

       do j = 1, i_num_bubbles
        xdistancesum(j)    = 0.0d0
        ydistancesum(j)    = 0.0d0
        zdistancesum(j)    = 0.0d0
        xvelsum(j)         = 0.0d0
        yvelsum(j)         = 0.0d0
        zvelsum(j)         = 0.0d0
        xcforcesum(j)      = 0.0d0
        ycforcesum(j)      = 0.0d0
        zcforcesum(j)      = 0.0d0
        xforcenewtsum(j)   = 0.0d0
        yforcenewtsum(j)    = 0.0d0
        zforcenewtsum(j)    = 0.0d0
        velwghtsum(j)       = 0.0d0
        bubvolsum(j)        = 0.0d0
        denssum(j)          = 0.0d0
        nzinBsum(j)         = 0

        dx_new(j) = avgxdistold(j) - xdistideal(j)
        dy_new(j) = avgydistold(j) - ydistideal(j)
        dz_new(j) = avgzdistold(j) - zdistideal(j)
       end do
! magnus - old code:
      !dx_new = avgxdistold - xdistideal
      !dy_new = avgydistold - ydistideal
      !dz_new = avgzdistold - zdistideal

! magnus - ignore this implementation to the whole domain  
      if(iCForz_where .eq. 1) then !apply cf to whole domain
        if (myrank.eq.master) write(*,*) 'The mpid controller does not work when applied to ',
     &        'the whole domain'   
!         y_c_f = ycfcoeff(1) * avgycf +
!     &           ycfcoeff(2) * ( avgycforceold +
!     &           ycfcoeff(3) * dy_new +
!     &           ycfcoeff(4) * avgyvelold +
!     &           ycfcoeff(5) * ddyvel +
!     &           ycfcoeff(6) * dy_new*abs(dy_new) +
!     &           ycfcoeff(7) * dy_new*dy_new*dy_new +
!     &           ycfcoeff(8) * avgyvelold*abs(avgyvelold) +
!     &           ycfcoeff(9) * avgyvelold*avgyvelold*avgyvelold )
     
!         x_c_f = xcfcoeff(1) * avgxcf +
!     &           xcfcoeff(2) * ( avgxcforceold +
!     &           xcfcoeff(3) * dx_new +
!     &           xcfcoeff(4) * avgxvelold +
!     &           xcfcoeff(5) * ddxvel +
!     &           xcfcoeff(6) * dx_new*abs(dx_new) +
!     &           xcfcoeff(7) * dx_new*dx_new*dx_new +
!     &           xcfcoeff(8) * avgxvelold*abs(avgxvelold) +
!     &           xcfcoeff(9) * avgxvelold*avgxvelold*avgxvelold ) +
!     &           xcfcoeff(10)* (dy_new-y_drag_flip) * 
!     &                          abs(dy_new-y_drag_flip)

!        z_c_f = zcfcoeff(1) * avgzcf +
!    &           zcfcoeff(2) * (avgzcforceold +
!     &           zcfcoeff(3) * dz_new +
!     &           zcfcoeff(4) * avgzvelold +
!     &           zcfcoeff(5) * ddzvel +
!     &           zcfcoeff(6) * dz_new*abs(dz_new) +
!     &           zcfcoeff(7) * dz_new*dz_new*dz_new +
!     &           zcfcoeff(8) * avgzvelold*abs(avgzvelold) +
!     &           zcfcoeff(9) * avgzvelold*avgzvelold*avgzvelold )
   
      else !apply cf to only bubble
       do bub_id=1, i_num_bubbles
         y_c_f(bub_id) = ycfcoeff(1) * avgycf(bub_id) +
     &                ycfcoeff(2) * ( avgycforceold(bub_id) -
     &                ycfcoeff(3) * dy_new(bub_id) -
     &                ycfcoeff(4) * avgyvelold(bub_id) -
     &                ycfcoeff(5) * ddyvel(bub_id) -
     &                ycfcoeff(6) * dy_new(bub_id)*abs(dy_new(bub_id)) -
     &                ycfcoeff(7) * dy_new(bub_id)*dy_new(bub_id)*dy_new(bub_id) -
     &                ycfcoeff(8) * avgyvelold(bub_id)*abs(avgyvelold(bub_id)) -
     &                ycfcoeff(9) * avgyvelold(bub_id)*avgyvelold(bub_id)*avgyvelold(bub_id))
   
            x_c_f(bub_id) = xcfcoeff(1) * avgxcf(bub_id) +
     &                xcfcoeff(2) * ( avgxcforceold(bub_id) -
     &                xcfcoeff(3) * dx_new(bub_id) -
     &                xcfcoeff(4) * avgxvelold(bub_id) -
     &                xcfcoeff(5) * ddxvel(bub_id) -
     &                xcfcoeff(6) * dx_new(bub_id)*abs(dx_new(bub_id)) -
     &                xcfcoeff(7) * dx_new(bub_id)*dx_new(bub_id)*dx_new(bub_id) -
     &                xcfcoeff(8) * avgxvelold(bub_id)*abs(avgxvelold(bub_id)) -
     &                xcfcoeff(9) * avgxvelold(bub_id)*avgxvelold(bub_id)*avgxvelold(bub_id)) -
     &                xcfcoeff(10)* (dy_new(bub_id)-y_drag_flip) * 
     &                              abs(dy_new(bub_id)-y_drag_flip)

            z_c_f(bub_id) = zcfcoeff(1) * avgzcf(bub_id) +
     &                zcfcoeff(2) * (avgzcforceold(bub_id) -
     &                zcfcoeff(3) * dz_new(bub_id) -
     &                zcfcoeff(4) * avgzvelold(bub_id) -
     &                zcfcoeff(5) * ddzvel(bub_id) -
     &                zcfcoeff(6) * dz_new(bub_id)*abs(dz_new(bub_id)) -
     &                zcfcoeff(7) * dz_new(bub_id)*dz_new(bub_id)*dz_new(bub_id) -
     &                zcfcoeff(8) * avgzvelold(bub_id)*abs(avgzvelold(bub_id)) -
     &                zcfcoeff(9) * avgzvelold(bub_id)*avgzvelold(bub_id)*avgzvelold(bub_id))
       end do        ! cycle through bubble ids
      end if !iCForz_where

      end !end control force calculating
c======================================================================
c======================================================================
c======================================================================
      subroutine CFASSY()
!----------------------------------------------------------------------
!       Called in elmgmr.f
!       the following loop is the one for all the blocks distributed
!       onto a single processor, inter-processor communication will come
!       after this.
!
!----------------------------------------------------------------------
      use cf_arrays   ! to access lift control arrays
      include "common.h"
   
   
      do i = 1, npro
       do j=1, i_num_bubbles       ! mpid
         xdistancesum(j) = xdistancesum(j) + cf_var(i,1,j)
         ydistancesum(j) = ydistancesum(j) + cf_var(i,2,j)
         zdistancesum(j) = zdistancesum(j) + cf_var(i,3,j)
         xvelsum(j)      = xvelsum(j) + cf_var(i,4,j)
         yvelsum(j)      = yvelsum(j) + cf_var(i,5,j)
         zvelsum(j)      = zvelsum(j) + cf_var(i,6,j)
         xcforcesum(j)   = xcforcesum(j) + cf_var(i,7,j)
         ycforcesum(j)   = ycforcesum(j) + cf_var(i,8,j)
         zcforcesum(j)   = zcforcesum(j) + cf_var(i,9,j)
         xforcenewtsum(j)= xforcenewtsum(j) + cf_var(i,10,j)
         yforcenewtsum(j)= yforcenewtsum(j) + cf_var(i,11,j)
         zforcenewtsum(j)= zforcenewtsum(j) + cf_var(i,12,j)
         velwghtsum(j)   = velwghtsum(j) + cf_var(i,13,j)
         bubvolsum(j)    = bubvolsum(j) + cf_var(i,14,j)
         denssum(j)      = denssum(j) + cf_var(i,15,j)
         if(cf_var(i,2,j).ne.0.0d0)then
            nzinBsum(j)  = nzinBsum(j) + 1
         end if
       end do ! mpid controller
      end do  ! elements
   
      end     !control force assembly
c======================================================================
c======================================================================
c======================================================================
      subroutine CFMPIprocess()
!
!----------------------------------------------------------------------
!       Called in elmgmr.f
!       This subroutine is used to process the control force information
!       from each mesh blocks and processors. There is a lot of MPI
!       communications.
!
!----------------------------------------------------------------------
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
   
      do j=1, i_num_bubbles
         if (numpe > 1) then
            call MPI_ALLREDUCE (xdistancesum(j), totalxdist(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (ydistancesum(j), totalydist(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (zdistancesum(j), totalzdist(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (nzinBsum(j), ntotnzinB(j), 1,
     &           MPI_INTEGER,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (velwghtsum(j), totalvelwght(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (xvelsum(j), totalxvel(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (yvelsum(j), totalyvel(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (zvelsum(j), totalzvel(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (xcforcesum(j), totalxcforce(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (ycforcesum(j), totalycforce(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (zcforcesum(j), totalzcforce(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (bubvolsum(j), totbubvol(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (denssum(j), totbubdens(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (yforcenewtsum(j), totyfnewtsum(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (xforcenewtsum(j), totxfnewtsum(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE (zforcenewtsum(j), totzfnewtsum(j), 1,
     &           MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)

            ycfnewtons(j) = totyfnewtsum(j)
            xcfnewtons(j) = totxfnewtsum(j)
            zcfnewtons(j) = totzfnewtsum(j)

            totalycforce(j) = totalycforce(j) / totbubvol(j)
            totalxcforce(j) = totalxcforce(j) / totbubvol(j)
            totalzcforce(j) = totalzcforce(j) / totbubvol(j)

            avgydistance(j) = totalydist(j) / real(ntotnzinB(j))
            avgxdistance(j) = totalxdist(j) / real(ntotnzinB(j))
            avgzdistance(j) = totalzdist(j) / real(ntotnzinB(j))

            avgyvel(j) = totalyvel(j) / totalvelwght(j)
            avgxvel(j) = totalxvel(j) / totalvelwght(j)
            avgzvel(j) = totalzvel(j) / totalvelwght(j)

            avgycforce(j) = totalycforce(j)
            avgxcforce(j) = totalxcforce(j)
            avgzcforce(j) = totalzcforce(j)

         else !for single processor case

            avgydistance(j) = ydistancesum(j) / real(nzinBsum(j))
            avgxdistance(j) = xdistancesum(j) / real(nzinBsum(j))
            avgzdistance(j) = zdistancesum(j) / real(nzinBsum(j))
            avgyvel(j) = yvelsum(j) / velwghtsum(j)
            avgxvel(j) = xvelsum(j) / velwghtsum(j)
            avgzvel(j) = zvelsum(j) / velwghtsum(j)
            avgycforce(j) = ycforcesum(j) / bubvolsum(j)
            avgxcforce(j) = xcforcesum(j) / bubvolsum(j)
            avgzcforce(j) = zcforcesum(j) / bubvolsum(j)

            write(*,*)'Control force implemented!'
            write(*,*)'mpi is not used'
            
            write(*,*)'For bubble ', j
            write(*,*)'average y distance: ', avgydistance(j)
            write(*,*)'average x distance: ', avgxdistance(j)
            write(*,*)'average z distance: ', avgzdistance(j)

            write(*,*)'average y velocity: ', avgyvel(j)
            write(*,*)'average x velocity: ', avgxvel(j)
            write(*,*)'average z velocity: ', avgzvel(j)

            write(*,*)'average y control force: ', avgycforce(j)
            write(*,*)'average x control force: ', avgxcforce(j)
            write(*,*)'average z control force: ', avgzcforce(j)
         endif
       end do ! mpid, bubble cycle

      end !end MPI processing
c======================================================================
c======================================================================
c======================================================================

      subroutine CFCollect(xx,        u1,     u2,     u3,
     &                       Sclr,      epsilon_ls_tmp, elemvol_local,
     &                       rho,       sforce)
!----------------------------------------------------------------------
!       Called in e3ivar.f
!       Control Force bubble information collection algorithm: control
!       force doesn't start updating until the 3rd timestep. at the
!       first timestep, there is no avgdistold to compare distideal too;
!       at the second timestep, dy_new will just equal zero because the
!       avgydistold of timestep 2 is the avgdist of timestep 1 which is
!       also the distideal
!
!----------------------------------------------------------------------
      use cf_arrays   ! to access lift control arrays
      use bub_track    !magnus, mpid
      include "common.h"

      dimension u1(npro),  u2(npro),  u3(npro)
      dimension xx(npro,nsd),         Sclr(npro),
     &            sforce(npro,3),       rho(npro)

      real*8    elemvol_local(ibksiz)
      real*8    epsilon_ls_tmp,       denswght

      integer bub_id


      cf_var = zero
      rholiq=datmat(1,1,1)
      rhogas=datmat(1,1,2)        

! below placed for mpid method testing
      !if (myrank .eq. master) write(*,*) 'block to check bubble element belongs to'
!      do j=1, npro
!        id= int(bub_info(j,11))
!        if (id .ne. 0 ) write(*,*) 'nonzero id found: ', id
        !do bub=1, i_num_bubbles
            !write(*,*) 'bub id is ', id
            !if (bub .eq. id) then
                !write(*,*)'element belongs to bubble on rank', myrank
            !end if
        !enddo
!      enddo

      do i=1,npro
         if(Sclr(i).le.epsilon_ls_tmp) then
            id = int(bub_info(i,11))    !mpid
            do bub_id=1, i_num_bubbles

              if (id .eq. bub_id) then !mpid
                     cf_var(i,1,id) = xx(i,1)
                     cf_var(i,2,id) = xx(i,2)
                     cf_var(i,3,id) = xx(i,3)
                     cf_var(i,4,id) = u1(i)*(rholiq-rho(i))*elemvol_local(i)
                     cf_var(i,5,id) = u2(i)*(rholiq-rho(i))*elemvol_local(i)
                     cf_var(i,6,id) = u3(i)*(rholiq-rho(i))*elemvol_local(i)
                     cf_var(i,7,id) = x_c_f(id)*elemvol_local(i)
                     cf_var(i,8,id) = y_c_f(id)*elemvol_local(i)
                     cf_var(i,9,id) = z_c_f(id)*elemvol_local(i)
                     cf_var(i,13,id)= (rholiq-rho(i))*elemvol_local(i)
                     cf_var(i,14,id)= elemvol_local(i)
                     cf_var(i,15,id)= rho(i)*elemvol_local(i)*
     &                     ((rholiq-rho(i))/(rholiq-rhogas))
              !velocity in bubble need be averaged by
              !(rho_l-rho(i))*vol(i)
              !densty in bub wghtd by elment vol and denswght
              end if !mpid, check element id
              end do
         end if
      end do

      !mpid note, ignore application to the whole domain for now
      if(iCForz_where .eq. 1) then !apply cf to whole domain
       write(*,*) 'mpid controller does not work with the whole domain option'
!----------------------------------------------------------------------
!"This block is ignored for the mpid controller, as the controller
!must be applied to each individual bubble" - magnus, Oct 2024
!         do i=1,npro
!            sforce(i,1) = sforce(i,1) + x_c_f
!            sforce(i,2) = sforce(i,2) + y_c_f
!            sforce(i,3) = sforce(i,3) + z_c_f
!            if(Sclr(i).le.epsilon_ls_tmp)then
            !Extract the Drag Force, Lift Force, and Z Force in [N]
!              cf_var(i,10)=(rholiq-rho(i))*
!     &                        elemvol_local(i)*(datmat(1,5,1)+x_c_f)
!              cf_var(i,11)=(rholiq-rho(i))*
!     &                        elemvol_local(i)*(datmat(2,5,1)+y_c_f)
!              cf_var(i,12)=(rholiq-rho(i))*
!     &                        elemvol_local(i)*(datmat(3,5,1)+z_c_f)
!            end if
!         end do
!----------------------------------------------------------------------
      else !apply cf inside bubble
         do i = 1, npro
            if(Sclr(i).le.epsilon_ls_tmp)then
                id = int(bub_info(i,11))    !mpid
                do bub_id=1, i_num_bubbles
                     if (id .eq. bub_id) then !mpid
                      denswght=(rholiq-rho(i))/(rholiq-rhogas)
                      sforce(i,1)= sforce(i,1) + x_c_f(id)*denswght/rho(i)
                      sforce(i,2)= sforce(i,2) + y_c_f(id)*denswght/rho(i)
                      sforce(i,3)= sforce(i,3) + z_c_f(id)*denswght/rho(i)
                      !Extract the Force in [N]
                      cf_var(i,10,id)=x_c_f(id)*denswght*elemvol_local(i)
                      cf_var(i,11,id)=y_c_f(id)*denswght*elemvol_local(i)
                      cf_var(i,12,id)=z_c_f(id)*denswght*elemvol_local(i)
                     end if
                end do
            end if
         end do
      end if !iCForz_where

      end     !CFCollect ends
c======================================================================
c======================================================================
c======================================================================

      subroutine CFpostprocs(istp,       numoldyhistind,
     &                         ixcf,       iycf,         izcf,
     &                         sumxcf_o,   sumycf_o,     sumzcf_o,     
     &                         sumycf_old, ycf_old_log)
c!
c!----------------------------------------------------------------------
c!Called in itrdrv.f
c!This subroutine is dealing with the post-processing for control forces
c!
c! Jun, Fall 2014
c!----------------------------------------------------------------------
      use cf_arrays !access to control force arrays
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      real*8 sumxcf, sumycf
      real*8 sumzcf
      real*8 sumxcf_o(i_num_bubbles), sumycf_o(i_num_bubbles)
      real*8 sumzcf_o(i_num_bubbles)
      real*8 sumycf_old

      integer numoldtsyhist
      integer numnewtsyhist
      integer ixcf, iycf, izcf
      integer bub_id

      logical ycf_old_log
      character*50:: lstepc
      character*2 :: id_char

!        if(myrank.eq.master) write(*,*) 'numoldyhistind = ',
!     &                                   numoldyhistind
      do j=1, i_num_bubbles
      if((i_res_cf .eq. 0).and.(istp .eq. 1))then

          ydistideal(j) = avgydistance(j)
          xdistideal(j) = avgxdistance(j)
          zdistideal(j) = avgzdistance(j)

      endif

      avgydistold(j) = avgydistance(j)
      avgxdistold(j) = avgxdistance(j)
      avgzdistold(j) = avgzdistance(j)

      if((i_res_cf .eq. 0).and.(istp .eq. 1))then
          ddyvel(j) = 0.0d0
          ddxvel(j) = 0.0d0
          ddzvel(j) = 0.0d0
      else
          ddyvel(j) = avgyvel(j) - avgyvelold(j)
          ddxvel(j) = avgxvel(j) - avgxvelold(j)
          ddzvel(j) = avgzvel(j) - avgzvelold(j)
      end if

      avgyvelold(j) = avgyvel(j)
      avgxvelold(j) = avgxvel(j)
      avgzvelold(j) = avgzvel(j)
      avgycforceold(j) = avgycforce(j)
      avgxcforceold(j) = avgxcforce(j)
      avgzcforceold(j) = avgzcforce(j)

      xcf(istp) = avgxcforce(j) -
     &      (xcfcoeff(10)*(dy_new(j)-y_drag_flip)*abs(dy_new(j)-y_drag_flip)) !remove coupled term from hist

      if((mod(istp-1,ntout).eq.0).and.(istp.ne.1))then
          ixcf=ixcf+ntout !this is done so that sumxcf_o does not double count values. pattern looks like
      endif!1:1,1:2,..,1:1*ntout,
           !1*ntout+1:1*ntout+1,1*ntout+1:1*ntout+2,..,1*ntout+1:2*ntout
           !2*ntout+1:2*ntout+1,2*ntout+1:2*ntout+2,..,2*ntout+1:3*ntout
      sumxcf = sum(xcf(ixcf:istp)) + sumxcf_o(j)
      avgxcf(j) = sumxcf / real(lstep,8)

!...best if numts_histyavg is < numxyzcflines. that way when doing
!restart using cf algorithm, you don't need to worry about retaining 
!data from the past several simulations. so if numts_histyavg = 500, 
!and you break simulation every 100 ts then you'll have to worry about
! retaining ycf data that spans over several nxyzcf.dat files, which 
!is currently not an available option for running restart control force
      ycf(istp,j) = avgycforce(j)
!        if(myrank.eq.master)write(*,*)'i_res_cf=',i_res_cf
!...Get start and end index of ycf from previous simulation
      if((i_res_cf.eq.1).and.(ycf_old_log.eqv..true.))then
          iyhistst = istp !lstep - numts_histyavg + 1
          if(istp.eq.1)then
             iyhisten = numoldyhistind  !lstep - 1
          end if
!        if(myrank.eq.master)write(*,*)'iyhistst, iyhisten',
!     &                                 iyhistst,iyhisten
!...Get the sum of the ycf array from previous simulation
          if(iyhistst.le.iyhisten)then
             sumycf_old=sum(ycf_old(iyhistst:iyhisten,j))
          end if
!        if(myrank.eq.master)write(*,*)'sumycf_old ',sumycf_old
          if(iyhistst.gt.iyhisten)then
             ycf_old_log = .false.
          end if
      end if !i_res_cf & ycf_old_log

!...compute current historical term based on current timestep and 
!overlap of previous simulation
      if(i_res_cf.eq.1)then
         if(ycf_old_log.eqv..true.)then
            numoldtsyhist = iyhisten - iyhistst + 1
            numnewtsyhist = numts_histyavg - numoldtsyhist
            sumycf = sumycf_old + sum(ycf(1:numnewtsyhist,j))
            avgycf(j) = sumycf / real(numts_histyavg,8)
!        if(myrank.eq.master)write(*,*)'numoldtsyhist, numnewtsyhist ',
!     &                                 numoldtsyhist,numnewtsyhist
!        if(myrank.eq.master)write(*,*)'sumycf ',sumycf
!        if(myrank.eq.master)write(*,*)'avgycf ',avgycf
!        if(myrank.eq.master)write(*,*)
         elseif(ycf_old_log.eqv..false.)then
                 avgycf(j) = sum(ycf(istp-numts_histyavg+1:istp,j)) /
     &                      real(numts_histyavg,8)
!        if(myrank.eq.master)write(*,*)'avgycf ',avgycf
!        if(myrank.eq.master)write(*,*)
         endif !ycf_old_log.eqv..true.

      elseif(i_res_cf.eq.0)then

             if(istp .ge. numts_histyavg)then
                avgycf(j) = sum(ycf(istp-numts_histyavg+1:istp,j)) /
     &                     real(numts_histyavg,8)

             else
                avgycf(j) = sum(ycf) / real(lstep,8)
             endif

      endif !i_res_cf

      zcf(istp) = avgzcforce(j)
      if((mod(istp-1,ntout).eq.0).and.(istp.ne.1))then
          izcf=izcf+ntout
      end if
      sumzcf = sum(zcf(izcf:istp)) + sumzcf_o(j)
      avgzcf(j) = sumzcf / real(lstep,8)

      totalycforceold(j) = totalycforce(j)
      totalxcforceold(j) = totalxcforce(j)
      totalzcforceold(j) = totalzcforce(j)
      
      if(myrank.eq.master)then
         write(7001+int(j)*10,736) time-delt(itseq), avgxdistance(j), 
     &                    avgydistance(j), avgzdistance(j)
         write(7002+int(j)*10,737) lstep, avgxvel(j), avgyvel(j), avgzvel(j)
         write(7003+int(j)*10,738) totalxcforce(j), totalycforce(j), totalzcforce(j)
         write(7004+int(j)*10,738) avgxcf(j), avgycf(j), avgzcf(j)
         write(7005+int(j)*10,739) dx_new(j),dy_new(j),dz_new(j),ddxvel(j),ddyvel(j),ddzvel(j)
         write(7006+int(j)*10,738) avgxcforce(j), avgycforce(j), avgzcforce(j)
         write(7007+int(j)*10,740) lstep, time-delt(itseq), xcfnewtons(j), 
     &                    ycfnewtons(j), zcfnewtons(j)
         write(7008+int(j)*10,738) xcf(istp), ycf(istp,j), zcf(istp)
      end if
!write the restart files for control forces
      if (mod(lstep, ntout) .eq. 0) then

          sumxcf_o(j) = sum(xcf(istp-ntout+1:istp)) + sumxcf_o(j)
          sumzcf_o(j) = sum(zcf(istp-ntout+1:istp)) + sumzcf_o(j)


          if(     lstep.ge.0 .and. lstep.le.9)then
              write(lstepc,'(i1.1)')lstep
          else if(lstep.ge.10 .and. lstep.le.99)then
              write(lstepc,'(i2.2)')lstep
          else if(lstep.ge.100 .and. lstep.le.999)then
              write(lstepc,'(i3.3)')lstep
          else if(lstep.ge.1000 .and. lstep.le.9999)then
              write(lstepc,'(i4.4)')lstep
          else if(lstep.ge.10000 .and. lstep.le.99999)then
              write(lstepc,'(i5.5)')lstep
          else if(lstep.ge.100000 .and. lstep.le.999999)then
              write(lstepc,'(i6.6)')lstep
          else if(lstep.ge.1000000 .and. lstep.le.9999999)then
              write(lstepc,'(i7.7)')lstep
          else if(lstep.ge.10000000 .and. lstep.le.99999999)then
              write(lstepc,'(i8.8)')lstep
          else if(lstep.ge.100000000 .and. lstep.le.999999999)then
              write(lstepc,'(i9.9)')lstep
          else if(lstep.ge.1000000000 .and. lstep.le.9999999999)then
              write(lstepc,'(i10.10)')
          end if

          if (myrank.eq.master) then
              bub_id = int(j)
             write(id_char, '(I0)') bub_id 
             OPEN(unit=8001+j*10,file='../controller/'//trim(id_char)//'_lift-results/restarts/'//
     &                    'restart_control_force'//
     &                    trim(lstepc)//'.dat',
     &         status="replace",action="write",iostat=ierror)

             IF(ierror /= 0) STOP "Error opening file 845 -- postproc"

        write(8001+j*10,733)avgycf(j),avgycforceold(j),dy_new(j),avgyvelold(j),
     &                ddyvel(j),avgydistold(j),ydistideal(j),sumycf_o(j),
     &                avgxcf(j),avgxcforceold(j),dx_new(j),avgxvelold(j),
     &                ddxvel(j),avgxdistold(j),xdistideal(j),sumxcf_o(j),
     &                avgzcf(j),avgzcforceold(j),dz_new(j),avgzvelold(j),
     &                ddzvel(j),avgzdistold(j),zdistideal(j),sumzcf_o(j),
     &                epsilon_lsd,delt(itseq),dtlset,CFLfl_max,
     &                CFLbuint_max,CFLls_max,vf_now,C_int_adjust
        write(*,*)' writing restart_control_force'//
     &                         trim(lstepc)//'.dat'
             CLOSE(8001+j*10)
          endif !myrank
      endif !mod(lstep, ntout)
       end do ! temporary end of mpid

 733  format(32(ES22.15, 3x))
 736  format(ES22.15, 2x, ES22.15, 2x, ES22.15, 2x, ES22.15)
 737  format(I6, 2x, ES22.15, 2x, ES22.15, 2x, ES22.15)
 738  format(ES22.15, 2x, ES22.15, 2x, ES22.15)
 739  format(ES22.15, 2x, ES22.15, 2x, ES22.15, 2x,
     1           ES22.15, 2x, ES22.15, 2x, ES22.15)
 740  format(I6, 2x, ES22.15, 2x, ES22.15, 2x, ES22.15, 2x, ES22.15)

      end !end post-processing
c======================================================================
c======================================================================
c======================================================================
        subroutine preprocess_lift_files(i_ts)
!----------------------------------------------------------------------
!       Called in subroutine CFrestar
!       Check to see which files need to be shortened based on file
!       length and initial time step of the restarted simulation. Then
!       replace the exisiting lift-files with the shortened data sets so
!       that the last line number in the replaced lift-file is equal to
!       the initial time step of the restarted simulation. 
! 
!----------------------------------------------------------------------
        real*8,allocatable :: avg_dis(:,:)
        real*8,allocatable :: avg_vel(:,:)
        real*8,allocatable :: tot_cfr(:,:)
        real*8,allocatable :: his_avg(:,:)
        real*8,allocatable :: pvd_iff(:,:)
        real*8,allocatable :: avg_cfr(:,:)
        real*8,allocatable :: nwt_cfr(:,:)
        real*8,allocatable :: xyz_cfr(:,:)

        integer, allocatable :: avg_vts(:)

        integer :: nl_ad
        integer :: nl_av
        integer :: nl_tc
        integer :: nl_ha
        integer :: nl_pv
        integer :: nl_ac
        integer :: nl_nc
        integer :: nl_cf

        integer :: i
        integer :: i_ts
        integer :: ierror, status

        character*2 :: id_char

        !There are 8 files to pre-process. Getting lengths of files
        do bub_id=1, i_num_bubbles
            
            write(id_char, '(I0)') bub_id
            call lcounter(nl_ad,'../controller/'//trim(id_char)//'_lift-results/average_distance.dat')
            call lcounter(nl_av,'../controller/'//trim(id_char)//'_lift-results/average_velocity.dat')
            call lcounter(nl_tc,'../controller/'//trim(id_char)//'_lift-results/total_control_force.dat')
            call lcounter(nl_ha,'../controller/'//trim(id_char)//'_lift-results/historical_average.dat')
            call lcounter(nl_pv,'../controller/'//trim(id_char)//'_lift-results/pos_vel_difference.dat')
            call lcounter(nl_ac,'../controller/'//trim(id_char)//'_lift-results/average_control_force.dat')
            call lcounter(nl_nc,'../controller/'//trim(id_char)//'_lift-results/newtons_control_force.dat')
            call lcounter(nl_cf,'../controller/'//trim(id_char)//'_lift-results/xyzcf.dat')

            if(i_ts.lt.nl_ad)then
            write(*,715)'Condensing average_distance.dat to ',i_ts,
     &               ' number of lines from ',nl_ad
            allocate(avg_dis(i_ts,4),stat=status)
            if(status /= 0) stop 'Error allocating avg_dis'

            open(unit=701,file='../controller/'//trim(id_char)//'_lift-results/average_distance.dat',
     &        status="old",action="read",iostat=ierror)
            if(ierror /= 0) stop "file 701 creation error"
            do i=1,i_ts
            read(701,711)avg_dis(i,1:4)
            end do
            close(701)

            open(unit=701,file='../controller/'//trim(id_char)//'_lift-results/average_distance.dat',
     &        status="replace",action="write",iostat=ierror)
            if(ierror /= 0) stop "file 701 creation error"
            do i=1,i_ts
            write(701,711)avg_dis(i,1:4)
            end do
            close(701)
            deallocate(avg_dis)
            end if

            if(i_ts.lt.nl_av)then
            write(*,715)'Condensing average_velocity.dat to ',i_ts,
     &               ' number of lines from ',nl_av
            allocate(avg_vel(i_ts,3),stat=status)
            if(status /= 0) stop 'Error allocating avg_vel'
            allocate(avg_vts(i_ts),stat=status)
            if(status /= 0) stop 'Error allocating avg_vts'

            open(unit=702,file='../controller/'//trim(id_char)//'_lift-results/average_velocity.dat',
     &        status="old",action="read",iostat=ierror)
            if(ierror /= 0) stop "file 702 creation error"
            do i=1,i_ts
            read(702)avg_vts(i),avg_vel(i,1:3)
            end do
            close(702)

            open(unit=702,file='../controller/'//trim(id_char)//'_lift-results/average_velocity.dat',
     &        status="replace",action="write",iostat=ierror)
            if(ierror /= 0) stop "file 702 creation error"
            do i=1,i_ts
            write(702,712)avg_vts(i),avg_vel(i,1:3)
            end do
            close(702)
            deallocate(avg_vts)
            deallocate(avg_vel)
            end if

            if(i_ts.lt.nl_tc)then
            write(*,715)'Condensing total_control_force.dat to ',i_ts,
     &               ' number of lines from ',nl_tc
            allocate(tot_cfr(i_ts,3),stat=status)
            if(status /= 0) stop 'Error allocating tot_cfr'

            open(unit=703,file='../controller/'//trim(id_char)//'_lift-results/total_control_force.dat',
     &        status="old",action="read",iostat=ierror)
            if(ierror /= 0) stop "file 703 creation error"
            do i=1,i_ts
            read(703,713)tot_cfr(i,1:3)
            end do
            close(703)

            open(unit=703,file='../controller/'//trim(id_char)//'_lift-results/total_control_force.dat',
     &        status="replace",action="write",iostat=ierror)
            if(ierror /= 0) stop "file 703 creation error"
            do i=1,i_ts
            write(703,713)tot_cfr(i,1:3)
            end do
            close(703)
            deallocate(tot_cfr)
            end if

            if(i_ts.lt.nl_ha)then
            write(*,715)'Condensing historical_average.dat to ',i_ts,
     &               ' number of lines from ',nl_ha
            allocate(his_avg(i_ts,3),stat=status)
            if(status /= 0) stop 'Error allocating his_avg'

            open(unit=704,file='../controller/'//trim(id_char)//'_lift-results/historical_average.dat',
     &        status="old",action="read",iostat=ierror)
            if(ierror /= 0) stop "file 704 creation error"
            do i=1,i_ts
            read(704,713)his_avg(i,1:3)
            end do
            close(704)

            open(unit=704,file='../controller/'//trim(id_char)//'_lift-results/historical_average.dat',
     &        status="replace",action="write",iostat=ierror)
            if(ierror /= 0) stop "file 704 creation error"
            do i=1,i_ts
            write(704,713)his_avg(i,1:3)
            end do
            close(704)
            deallocate(his_avg)
            end if

            if(i_ts.lt.nl_pv)then
            write(*,715)'Condensing pos_vel_difference.dat to ',i_ts,
     &               ' number of lines from ',nl_pv
            allocate(pvd_iff(i_ts,6),stat=status)
            if(status /= 0) stop 'Error allocating pvd_iff'

            open(unit=705,file='../controller/'//trim(id_char)//'_lift-results/pos_vel_difference.dat',
     &       status="old",action="read",iostat=ierror)
            if(ierror /= 0) stop "file 705 creation error"
            do i=1,i_ts
            read(705,714)pvd_iff(i,1:6)
            end do
            close(705)

            open(unit=705,file='../controller/'//trim(id_char)//'_lift-results/pos_vel_difference.dat',
     &       status="replace",action="write",iostat=ierror)
            if(ierror /= 0) stop "file 705 creation error"
            do i=1,i_ts
            write(705,714)pvd_iff(i,1:6)
            end do
            close(705)
            deallocate(pvd_iff)
            end if

            if(i_ts.lt.nl_ac)then
            write(*,715)'Condensing average_control_force.dat to ',i_ts,
     &               ' number of lines from ',nl_ac
            allocate(avg_cfr(i_ts,3),stat=status)
            if(status /= 0) stop 'Error allocating avg_cfr'

            open(unit=706,file='../controller/'//trim(id_char)//'_lift-results/average_control_force.dat',
     &       status="old",action="read",iostat=ierror)
            if(ierror /= 0) stop "file 706 creation error"
            do i=1,i_ts
            read(706,713)avg_cfr(i,1:3)
            end do
            close(706)

            open(unit=706,file='../controller/'//trim(id_char)//'_lift-results/average_control_force.dat',
     &       status="replace",action="write",iostat=ierror)
            if(ierror /= 0) stop "file 706 creation error"
            do i=1,i_ts
            write(706,713)avg_cfr(i,1:3)
            end do
            close(706)
            deallocate(avg_cfr)
            end if

            if(i_ts.lt.nl_nc)then
            write(*,715)'Condensing newtons_control_force.dat to ',i_ts,
     &               ' number of lines from ',nl_nc
            allocate(nwt_cfr(i_ts,3),stat=status)
            if(status /= 0) stop 'Error allocating nwt_cfr'

            open(unit=707,file='../controller/'//trim(id_char)//'_lift-results/newtons_control_force.dat',
     &        status="old",action="read",iostat=ierror)
            if(ierror /= 0) stop "file 707 creation error"
            do i=1,i_ts
            read(707,716)nwt_cfr(i,1:3)
            end do
            close(707)

            open(unit=707,file='../controller/'//trim(id_char)//'_lift-results/newtons_control_force.dat',
     &        status="replace",action="write",iostat=ierror)
            if(ierror /= 0) stop "file 707 creation error"
            do i=1,i_ts
            write(707,716)nwt_cfr(i,1:3)
            end do
            close(707)
            deallocate(nwt_cfr)
            end if

            if(i_ts.lt.nl_cf)then
            write(*,715)'Condensing xyzcf.dat to ',i_ts,
     &               ' number of lines from ',nl_cf
            allocate(xyz_cfr(i_ts,3),stat=status)
            if(status /= 0) stop 'Error allocating xyz_cfr'

            open(unit=708,file='../controller/'//trim(id_char)//'_lift-results/xyzcf.dat',
     &        status="old",action="read",iostat=ierror)
            if(ierror /= 0) stop "file 708 creation error"
            do i=1,i_ts
            read(708,713)xyz_cfr(i,1:3)
            end do
            close(708)

            open(unit=708,file='../controller/'//trim(id_char)//'_lift-results/xyzcf.dat',
     &        status="replace",action="write",iostat=ierror)
            if(ierror /= 0) stop "file 708 creation error"
            do i=1,i_ts
            write(708,713)xyz_cfr(i,1:3)
            end do
            close(708)
            deallocate(xyz_cfr)
            end if

711   format(E22.15, 2x, E22.15, 2x, E22.15, 2x, E22.15)
712   format(I6, 2x, E22.15, 2x, E22.15, 2x, E22.15)
713   format(E22.15, 2x, E22.15, 2x, E22.15)
714   format(E22.15, 2x, E22.15, 2x, E22.15, 2x,
     &       E22.15, 2x, E22.15, 2x, E22.15)
715   format(A, I6, A, I6)
716   format(I6, 2x, E22.15, 2x, E22.15, 2x, E22.15, 2x, E22.15)

      end do    !mpid, bubble loop
        return
        end
c======================================================================
c======================================================================
c======================================================================

      subroutine lcounter(LineCount,filename)
!----------------------------------------------------------------------
!       Called in CFrestar and preprocess_lift_files
!       This subroutine will count the number of lines in a certain data
!       file
! 
!----------------------------------------------------------------------
         use, intrinsic :: iso_fortran_env

         integer, intent(out) :: LineCount
         integer :: Read_Code
         character (len=200) :: line
         character (len=200) :: filename

         open(unit=51,file=trim(filename),
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

      end
c======================================================================
c======================================================================
c======================================================================