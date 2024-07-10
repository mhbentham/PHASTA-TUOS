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

