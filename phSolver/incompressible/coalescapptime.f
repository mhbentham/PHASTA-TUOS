         subroutine CoalescAppTime(avgxcoordf, avgycoordf, avgzcoordf,
     &                             avgxcoordold2, avgycoordold2,
     &                             avgzcoordold2, app_time, itrtimestp)
c
c!----------------------------------------------------------------------
c
c! This routine assigns the new center coordinates for the coalescence
c! control and tracks the amount of time the coalescence force has
c! been active.
c
c! Matt Talley, Winter 2014.
c!----------------------------------------------------------------------
c
        use pvsQbi  ! brings in NABI
        use stats   !
        use pointer_data  ! brings in the pointers for the blocked arrays
        use local_mass
        use spat_var_eps
        use timedata  ! for iblkts usage

        include "common.h"
        include "mpif.h"

        real*8 avgxcoordf(coalest), avgycoordf(coalest),
     &         avgzcoordf(coalest), avgxcoordold2(coalest),
     &         avgycoordold2(coalest), avgzcoordold2(coalest),
     &         avgcoorddist(coalest,coalest),
     &         app_time(coalest,2)

        real*8 itrtimestp, coalesc_time

        integer event_tag(coalest,coalest)

c!.... Initialize variables
        event_tag(:,:) = 0
        app_time(:,1) = 0.0d0
        avgcoorddist(:,:) = 1.0d4
        coalesc_time = 0.0d0

c!.... Calculate the maximum coalescence time
        coalesc_time = sqrt((((coalbubrad)**3) * datmat(1,1,1))
     &                 / (16.0d0 *(1/Bo))) * log(1.0d-4/1.0d-8)

        do k1 = 1, coalest
           avgxcoordold(k1) = avgxcoordf(k1)
           avgycoordold(k1) = avgycoordf(k1)
           avgzcoordold(k1) = avgzcoordf(k1)

c!....Track the amount of time the Coalescence Control Algorithm has
c!....active for each different event

           if (app_time(k1,2).le.coalesc_time) then
              if (avgxcoordold2(k1).gt.-1.0d3) then
                 do k2 = 1, coalest
!                    if (avgxcoordold(k2).gt.-1.0d3) then
                       avgcoorddist(k1,k2) = sqrt((avgxcoordold(k2) -
     &                 avgxcoordold2(k1))**2 + (avgycoordold(k2) -
     &                 avgycoordold2(k1))**2 + (avgzcoordold(k2) -
     &                 avgzcoordold2(k1))**2)
!                    endif
                 enddo ! k2

                 do k2 = 1, coalest
                    if ((avgcoorddist(k1,k2).lt.coalbubrad).and.
     &                 (event_tag(k1,k2).eq.0)) then
                       app_time(k2,1) = app_time(k1,2) + itrtimestp/2.0d0
                       event_tag(k1,:) = 1
                       event_tag(:,k2) = 1
                       coalcon_rem(k1) = 0

                       if (myrank.eq.master) write(*,*) 'Coalescence',
     &                 ' Event #: ',k1,' to ',k2
                       if (myrank.eq.master) write(*,*) 'x average',
     &                 ' position:', avgxcoordold(k2)
                       if (myrank.eq.master) write(*,*) 'y average',
     &                 ' position:', avgycoordold(k2)
                       if (myrank.eq.master) write(*,*) 'z average',
     &                 ' position:', avgzcoordold(k2)
                    endif
                 enddo ! k2

                 do k2 = 1, coalest
                    if (event_tag(k1,k2).eq.0) then
                       if ((avgcoorddist(k1,k2).gt.coalbubrad).and.
     &                    (avgcoorddist(k1,k2).lt.1.0d4)) then
                          app_time(k1,2) = 0.0d0
                          coalcon_rem(k1) = 1
                          event_tag(k1,:) = 1

                          if (myrank.eq.master) write(*,*) 'Old',
     &                    ' Coalescence Event #: ',k1,' has ended',
     &                    ' because they bounced off one another'
                       endif
                    endif
                 enddo ! k2
              endif

           else

              if (avgxcoordold2(k1).gt.-1.0d3) then
                 do k2 = 1, coalest
!                    if (avgxcoordold(k2).gt.-1.0d3) then
                       avgcoorddist(k1,k2) = sqrt((avgxcoordold(k2) -
     &                 avgxcoordold2(k1))**2 + (avgycoordold(k2) -
     &                 avgycoordold2(k1))**2 + (avgzcoordold(k2) -
     &                 avgzcoordold2(k1))**2)
!                    endif
                 enddo ! k2

                 do k2 = 1, coalest
                    if ((avgcoorddist(k1,k2).lt.coalbubrad).and.
     &                 (event_tag(k1,k2).eq.0)) then
                       app_time(k2,1) = app_time(k1,2) + itrtimestp/2.0d0
                       event_tag(k1,:) = 1
                       event_tag(:,k2) = 1
                       coalcon_rem(k1) = 1

                       if (myrank.eq.master) write(*,*) 'Coalescence',
     &                 ' Event #: ',k1,' to ',k2,' has exceeded the',
     &                 ' drainage time and the force is being removed.'
                    endif
                 enddo ! k2

                 do k2 = 1, coalest
                    if (event_tag(k1,k2).eq.0) then
                       if ((avgcoorddist(k1,k2).gt.coalbubrad).and.
     &                    (avgcoorddist(k1,k2).lt.1.0d4)) then
                          app_time(k1,2) = 0.0d0
                          coalcon_rem(k1) = 1
                          event_tag(k1,:) = 1

                          if (myrank.eq.master) write(*,*) 'Old',
     &                    ' Coalescence Event #: ',k1,' has ended'
                       endif
                    endif
                 enddo ! k2
              endif
           endif !(app_time)

           if (app_time(k1,2).le.coalesc_time) then
              if (avgxcoordold2(k1).le.-1.0d3) then
                 do k2 = 1, coalest
                    if ((event_tag(k1,k2).eq.0).and.
     &              (avgxcoordold(k2).gt.-1.0d3)) then

                       app_time(k2,1) = itrtimestp/2.0d0
                       event_tag(:,k2) = 1
                       coalcon_rem(k1) = 0

                       if (myrank.eq.master) write(*,*) 'New',
     &                 ' Coalescence Event #: ',k2
                       if (myrank.eq.master) write(*,*) 'x average',
     &                 ' position:', avgxcoordold(k2)
                       if (myrank.eq.master) write(*,*) 'y average',
     &                 ' position:', avgycoordold(k2)
                       if (myrank.eq.master) write(*,*) 'z average',
     &                 ' position:', avgzcoordold(k2)

                    endif
                 enddo !k2
              endif
           endif
        enddo ! k1

        app_time(:,2) = app_time(:,1)
        avgxcoordold2(:) = avgxcoordold(:)
        avgycoordold2(:) = avgycoordold(:)
        avgzcoordold2(:) = avgzcoordold(:)

!        if (myrank.eq.master) then
!           do k = 1, coalest
!              write(7906,207) k, app_time(k,1), app_time(k,2),
!     &        coalcon_rem(k), avgxcoordold(k), avgycoordold(k),
!     &        avgzcoordold(k), itrtimestp, coalesc_time
!           enddo
!        endif

!207   format(I4,1x,ES24.16,1x,ES24.16,1x,I4,1x,ES24.16,1x,ES24.16,1x,
!     &       ES24.16,1x,ES24.16,1x,ES24.16)

        end
