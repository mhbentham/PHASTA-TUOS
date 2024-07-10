         subroutine CoalescCon (xarray, yarray, zarray, coordtag,
     &                          bubradius, bubradius2, avgxcoordf,
     &                          avgycoordf, avgzcoordf)
c
c!----------------------------------------------------------------------
c
c! This routine computes the center coordinates for coalescence events
c! during the simulation.
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
c
        include "common.h"
        include "mpif.h"
c
        real*8 xarray(ibksiz), yarray(ibksiz), zarray(ibksiz)
        integer coordtag(ibksiz) !Passed arrays from e3ivar

        real*8 avgxcoord, avgycoord, avgzcoord, avgvectdist,
!Coalescence control center pt
     &         avgxcoordf(coalest), avgycoordf(coalest),
     &         avgzcoordf(coalest)

        real*8 totxcoordsum, totycoordsum, totzcoordsum, totvectdistsum,
     &         totxcoordsum_mult(coalest), totycoordsum_mult(coalest),
     &         totzcoordsum_mult(coalest) !Total sum of coordinates

        integer totcoordcount, totvectnumbsum,
     &          totcoordcount_mult(coalest) !Total number or coordinates

        real*8 xcoordsum, ycoordsum, zcoordsum, vectdistsum,
     &         xcoordsum_mult(coalest), ycoordsum_mult(coalest),
     &         zcoordsum_mult(coalest) !Sum from each processor

        integer totcoordnumb, vectnumbsum, diffnumbsum,
     &          totcoordnumb_mult(coalest) !Total number from each processor

        real*8 globalxcoord(ibksiz,nelblk), globalycoord(ibksiz,nelblk),
     &         globalzcoord(ibksiz,nelblk), vectdist(ibksiz,nelblk),
     &         xvectcoord(ibksiz,nelblk), yvectcoord(ibksiz,nelblk),
     &         zvectcoord(ibksiz,nelblk), vectangle(ibksiz,nelblk),
     &         vectdist2(ibksiz,nelblk) !Arrays from each processor

        real*8 dotmax(ibksiz,nelblk), vectdist_max_tmp, vect_max_xcoord,
     &         vect_max_ycoord, vect_max_zcoord, vect_max_xcoord_tmp,
     &         vect_max_ycoord_tmp, vect_max_zcoord_tmp

        real*8 phi_max, bubradius, bubradius2, length_bside_tri,
     &         angle1, angle2, hypot_len_1, hypot_len_2,
     &         avgcoordfdist(coalest,coalest)

        integer totcoordnumbarray(ibksiz,nelblk),
     &          vectnumb(ibksiz,nelblk)

        integer intone, vect_max_i, vect_max_iblk,
     &          sign_of_vect_max_xcoord, sign_of_vect_max_ycoord,
     &          sign_of_vect_max_zcoord, sign_of_vect_max_xcoord_tmp,
     &          sign_of_vect_max_ycoord_tmp, sign_of_vect_max_zcoord_tmp

        integer coalesc_tag(ibksiz,nelblk), consol_tag(coalest),
     &          avgcoordf_erase_tag(coalest)

c! Initialize bubble coalecence control variables

        xcoordsum = 0.0d0 !Matt T.
        ycoordsum = 0.0d0
        zcoordsum = 0.0d0
        totcoordnumb = 0

        globalxcoord(:,:) = zero
        globalycoord(:,:) = zero
        globalzcoord(:,:) = zero
        totcoordnumbarray(:,:) = 0

        avgxcoord = -1.0d3
        avgycoord = -1.0d3
        avgzcoord = -1.0d3

        totxcoordsum = 0.0d0
        totycoordsum = 0.0d0
        totzcoordsum = 0.0d0
        totcoordcount = 0

        do iblk = 1, nelblk
           iel    = lcblk(1,iblk)
           npro   = lcblk(1,iblk+1) - iel

c!....Storing the coordinates into each globalcoord array and counting how many
c!....points are in each array. Then summing up each x, y, and z coordinate
c!....array to get a sum of each x, y, and z coordinate

           do i = 1, npro !Matt T.
              globalxcoord(i,iblk) = xarray(i)
              globalycoord(i,iblk) = yarray(i)
              globalzcoord(i,iblk) = zarray(i)
              totcoordnumbarray(i,iblk) = coordtag(i)
           enddo

           do m = 1, npro
              xcoordsum = xcoordsum + globalxcoord(m,iblk)
              ycoordsum = ycoordsum + globalycoord(m,iblk)
              zcoordsum = zcoordsum + globalzcoord(m,iblk)
              totcoordnumb = totcoordnumb + totcoordnumbarray(m,iblk)
           enddo
        enddo

         if (numpe.gt.1) then   !Matt T.
            call MPI_ALLREDUCE(xcoordsum, totxcoordsum, 1,
     &           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(ycoordsum, totycoordsum, 1,
     &           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(zcoordsum, totzcoordsum, 1,
     &           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(totcoordnumb, totcoordcount, 1,
     &           MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

            if (totcoordcount.gt.0) then
               avgxcoord = totxcoordsum / DBLE(totcoordcount) !Matt T.
               avgycoord = totycoordsum / DBLE(totcoordcount)
               avgzcoord = totzcoordsum / DBLE(totcoordcount)
            endif
         else
            if (totcoordnumb.gt.0) then
               avgxcoord = xcoordsum / DBLE(totcoordnumb) !Matt T.
               avgycoord = ycoordsum / DBLE(totcoordnumb)
               avgzcoord = zcoordsum / DBLE(totcoordnumb)
            endif
         endif

c!....Begin the check for multiple coalescence events

         vectdistsum = 0.0d0
         vectnumbsum = 0
         vectdist(:,:) = 0.0d0
         vectnumb(:,:) = 0

         totvectdistsum = 0
         avgvectdist = 0.0d0
         totvectnumbsum = 0
         intone = 1

        do iblk = 1, nelblk
           iel    = lcblk(1,iblk)
           npro   = lcblk(1,iblk+1) - iel

           do i = 1, npro
              if (totcoordnumbarray(i,iblk).eq.intone) then

                 vectdist(i,iblk) = sqrt((globalxcoord(i,iblk) -
     &           avgxcoord)**2 + (globalycoord(i,iblk) - avgycoord)**2
     &           + (globalzcoord(i,iblk) - avgzcoord)**2)

                 vectnumb(i,iblk) = 1

              endif
           enddo ! i

           do m = 1, npro
              vectdistsum = vectdistsum + vectdist(m,iblk)
              vectnumbsum = vectnumbsum + vectnumb(m,iblk)
           enddo ! m
        enddo

         if (numpe.gt.1) then
            call MPI_ALLREDUCE(vectdistsum, totvectdistsum, 1,
     &           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE(vectnumbsum, totvectnumbsum, 1,
     &           MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

            if (totvectnumbsum.gt.0) then
               avgvectdist = totvectdistsum / DBLE(totvectnumbsum)
            endif
         else
            if (vectnumbsum.gt.0) then
               avgvectdist = vectdistsum / DBLE(vectnumbsum)
            endif
         endif

c!....Initialize new variables for the angle between vectors

         coalesc_tag(:,:) = 0

         xvectcoord(:,:) = 0.0d0
         yvectcoord(:,:) = 0.0d0
         zvectcoord(:,:) = 0.0d0

         xcoordsum_mult(:) = 0.0d0
         ycoordsum_mult(:) = 0.0d0
         zcoordsum_mult(:) = 0.0d0
         totcoordnumb_mult(:) = 0

         totxcoordsum_mult(:) = 0.0d0
         totycoordsum_mult(:) = 0.0d0
         totzcoordsum_mult(:) = 0.0d0
         totcoordcount_mult(:) = 0

         if ((avgxcoord.gt.-1.0d3).and.(avgvectdist.ge.bubradius)) then

            do iblk = 1, nelblk
               iel    = lcblk(1,iblk)
               npro   = lcblk(1,iblk+1) - iel

               do i = 1, npro
                  if (totcoordnumbarray(i,iblk).eq.intone) then
                     xvectcoord(i,iblk) = globalxcoord(i,iblk) -
     &               avgxcoord
                     yvectcoord(i,iblk) = globalycoord(i,iblk) -
     &               avgycoord
                     zvectcoord(i,iblk) = globalzcoord(i,iblk) -
     &               avgzcoord
                  endif
               enddo
            enddo

            do k = 1, coalest

c!....Values to be re-initialized for each coalescence
               sign_of_vect_max_xcoord = 0
               sign_of_vect_max_ycoord = 0
               sign_of_vect_max_zcoord = 0
               sign_of_vect_max_xcoord_tmp = 0
               sign_of_vect_max_ycoord_tmp = 0
               sign_of_vect_max_zcoord_tmp = 0

               vect_max_xcoord = 0.0d0
               vect_max_ycoord = 0.0d0
               vect_max_zcoord = 0.0d0
               vect_max_xcoord_tmp = 0.0d0
               vect_max_ycoord_tmp = 0.0d0
               vect_max_zcoord_tmp = 0.0d0

               vectdist_max = 0.0d0
               vectdist_max_tmp = 0.0d0
               vect_max_i = 0
               vect_max_iblk = 0

               vectangle(:,:) = 0.0d0
               dotmax(:,:) = 0.0d0
               phi_max = 0.0d0
               length_bside_tri = 0.0d0
               vectdist2(:,:) = 0.0d0
               angle1 = 0.0d0
               angle2 = 0.0d0
               hypot_len_1 = 0.0d0
               hypot_len_2 = 0.0d0

               do iblk = 1, nelblk
                  iel    = lcblk(1,iblk)
                  npro   = lcblk(1,iblk+1) - iel

                  do i = 1, npro
                     if ((coalesc_tag(i,iblk).eq.0).and.
     &               (vectdist(i,iblk).gt.vectdist_max)) then
                        vectdist_max = vectdist(i,iblk)
                     endif
                  enddo
               enddo

               if (numpe.gt.1) then

                  call MPI_ALLREDUCE (vectdist_max,vectdist_max_tmp,1,
     &            MPI_DOUBLE_PRECISION,MPI_MAX, MPI_COMM_WORLD,ierr)
               else
                  vectdist_max_tmp = vectdist_max
               endif

               vectdist_max = vectdist_max_tmp

               do iblk = 1, nelblk
                  iel    = lcblk(1,iblk)
                  npro   = lcblk(1,iblk+1) - iel

                  do i = 1, npro
                     if ((coalesc_tag(i,iblk).eq.0).and.
     &         (abs(vectdist(i,iblk) - vectdist_max).lt.1.0d-24)) then
                        vect_max_i = i
                        vect_max_iblk = iblk !Values only known to one processor

                        exit

                     endif
                  enddo
               enddo

               if ((vect_max_i.gt.0).and.(vect_max_iblk.gt.0)) then

                  vect_max_xcoord = xvectcoord(vect_max_i,vect_max_iblk)
                  vect_max_ycoord = yvectcoord(vect_max_i,vect_max_iblk)
                  vect_max_zcoord = zvectcoord(vect_max_i,vect_max_iblk)

               endif

               if (vect_max_xcoord.lt.0.0d0) then
                  sign_of_vect_max_xcoord = 1
               endif

               if (vect_max_ycoord.lt.0.0d0) then
                  sign_of_vect_max_ycoord = 1
               endif

               if (vect_max_zcoord.lt.0.0d0) then
                  sign_of_vect_max_zcoord = 1
               endif

               call MPI_ALLREDUCE (abs(vect_max_xcoord),
     &              vect_max_xcoord_tmp,1,MPI_DOUBLE_PRECISION,MPI_MAX,
     &              MPI_COMM_WORLD,ierr)
               call MPI_ALLREDUCE (abs(vect_max_ycoord),
     &              vect_max_ycoord_tmp,1,MPI_DOUBLE_PRECISION,MPI_MAX,
     &              MPI_COMM_WORLD,ierr)
               call MPI_ALLREDUCE (abs(vect_max_zcoord),
     &              vect_max_zcoord_tmp,1,MPI_DOUBLE_PRECISION,MPI_MAX,
     &              MPI_COMM_WORLD,ierr)
               call MPI_ALLREDUCE (sign_of_vect_max_xcoord,
     &              sign_of_vect_max_xcoord_tmp,1,MPI_INTEGER,MPI_MAX,
     &              MPI_COMM_WORLD,ierr)
               call MPI_ALLREDUCE (sign_of_vect_max_ycoord,
     &              sign_of_vect_max_ycoord_tmp,1,MPI_INTEGER,MPI_MAX,
     &              MPI_COMM_WORLD,ierr)
               call MPI_ALLREDUCE (sign_of_vect_max_zcoord,
     &              sign_of_vect_max_zcoord_tmp,1,MPI_INTEGER,MPI_MAX,
     &              MPI_COMM_WORLD,ierr)

               if (sign_of_vect_max_xcoord_tmp.eq.intone) then
                  vect_max_xcoord = -vect_max_xcoord_tmp
               else
                  vect_max_xcoord = vect_max_xcoord_tmp
               endif

               if (sign_of_vect_max_ycoord_tmp.eq.intone) then
                  vect_max_ycoord = -vect_max_ycoord_tmp
               else
                  vect_max_ycoord = vect_max_ycoord_tmp
               endif

               if (sign_of_vect_max_zcoord_tmp.eq.intone) then
                  vect_max_zcoord_tmp = -vect_max_zcoord_tmp
               else
                  vect_max_zcoord = vect_max_zcoord_tmp
               endif

               do iblk = 1, nelblk
                  iel    = lcblk(1,iblk)
                  npro   = lcblk(1,iblk+1) - iel

                  do i = 1, npro
                     if ((coalesc_tag(i,iblk).eq.0).and.
     &                  (totcoordnumbarray(i,iblk).eq.intone)) then
                        dotmax(i,iblk) = xvectcoord(i,iblk)
     &                  * vect_max_xcoord
     &                  + yvectcoord(i,iblk) * vect_max_ycoord
     &                  + zvectcoord(i,iblk) * vect_max_zcoord

                        if ((dotmax(i,iblk) / (vectdist(i,iblk)
     &                  * vectdist_max)).gt.1.0d0) then
                           vectangle(i,iblk) = acos(1.0d0)
                        else if ((dotmax(i,iblk) / (vectdist(i,iblk)
     &                  * vectdist_max)).lt.-1.0d0) then
                           vectangle(i,iblk) = acos(-1.0d0)
                        else
                           vectangle(i,iblk) = acos(dotmax(i,iblk)
     &                     / (vectdist(i,iblk) * vectdist_max))
                        endif
                     endif
                  enddo
               enddo

c!....Determine the distance between each vector and vector_max
               do iblk = 1, nelblk
                  iel    = lcblk(1,iblk)
                  npro   = lcblk(1,iblk+1) - iel

                  do i = 1, npro
                     if ((coalesc_tag(i,iblk).eq.0).and.
     &                  (totcoordnumbarray(i,iblk).eq.intone)) then
                        vectdist2(i,iblk) = sqrt((xvectcoord(i,iblk) -
     &                     vect_max_xcoord)**2 + (yvectcoord(i,iblk) -
     &                     vect_max_ycoord)**2 + (zvectcoord(i,iblk) -
     &                     vect_max_zcoord)**2)
                     endif
                  enddo ! i
               enddo ! iblk

c!....Initialze the first tag based on the proper processor and tag the rest
               if ((vect_max_i.gt.0).and.(vect_max_iblk.gt.0)) then
                  coalesc_tag(vect_max_i,vect_max_iblk) = k
               endif

c!....Determine the maximum angle and assign event numbers
               length_bside_tri = sqrt(vectdist_max**2 -
     &         (2.0d0*bubradius)**2)

               if (vectdist_max**2.ge.(2.0d0*bubradius)**2) then

                  do iblk = 1, nelblk
                     iel    = lcblk(1,iblk)
                     npro   = lcblk(1,iblk+1) - iel

                     do i = 1, npro
                        if ((totcoordnumbarray(i,iblk).eq.intone).and.
     &                     (vectdist(i,iblk).ge.2.0d0*bubradius)) then
                           phi_max = acos(length_bside_tri /
     &                        vectdist_max)
                        endif

                        if ((totcoordnumbarray(i,iblk).eq.intone).and.
     &                     (vectdist(i,iblk).lt.2.0d0*bubradius)) then
                           phi_max = asin((2.0d0*bubradius) /
     &                        vectdist_max)
                        endif

                        if ((coalesc_tag(i,iblk).eq.0).and.
     &                     (totcoordnumbarray(i,iblk).eq.intone)) then
                           if ((abs(vectdist(i,iblk)-vectdist_max).le.
     &                        (2.0d0*bubradius)).and.
     &                        (vectangle(i,iblk).le.phi_max)) then

                              coalesc_tag(i,iblk) = k
                           endif
                        endif
                     enddo ! i
                  enddo ! iblk

               else

c!....Calculate the maximum angle at the 1st epsilon contour
                  if (vectdist_max.gt.0.0d0) then
                     angle1 = acos(bubradius2 / vectdist_max)
                     hypot_len_1 = vectdist_max*sin(angle1)
                     hypot_len_2 = 2.0d0*bubradius - hypot_len_1
                     angle2 = atan(hypot_len_2 / bubradius2)

                     phi_max = angle1 + angle2

                     do iblk = 1, nelblk
                        iel    = lcblk(1,iblk)
                        npro   = lcblk(1,iblk+1) - iel

                        do i = 1, npro
                           if ((coalesc_tag(i,iblk).eq.0).and.
     &                        (totcoordnumbarray(i,iblk).eq.intone)) then
                              if ((vectdist2(i,iblk).le.
     &                           (2.0d0*bubradius)).and.
     &                           (vectangle(i,iblk).le.phi_max)) then

                                 coalesc_tag(i,iblk) = k
                              endif
                           endif
                        enddo ! i
                     enddo ! iblk
                  endif
               endif

               do iblk = 1, nelblk
                  iel    = lcblk(1,iblk)
                  npro   = lcblk(1,iblk+1) - iel

                  do m = 1, npro
                     if (coalesc_tag(m,iblk).eq.k) then
                        xcoordsum_mult(k) = xcoordsum_mult(k) +
     &                  globalxcoord(m,iblk)
                        ycoordsum_mult(k) = ycoordsum_mult(k) +
     &                  globalycoord(m,iblk)
                        zcoordsum_mult(k) = zcoordsum_mult(k) +
     &                  globalzcoord(m,iblk)
                        totcoordnumb_mult(k) = totcoordnumb_mult(k) +
     &                  totcoordnumbarray(m,iblk)
                     endif
                  enddo ! m
               enddo ! iblk

               if (numpe > 1) then   !Matt T.
                  call MPI_ALLREDUCE(xcoordsum_mult(k),
     &                 totxcoordsum_mult(k),1,MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, MPI_COMM_WORLD,ierr)
                  call MPI_ALLREDUCE(ycoordsum_mult(k),
     &                 totycoordsum_mult(k),1,MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, MPI_COMM_WORLD,ierr)
                  call MPI_ALLREDUCE(zcoordsum_mult(k),
     &                 totzcoordsum_mult(k),1,MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, MPI_COMM_WORLD,ierr)
                  call MPI_ALLREDUCE(totcoordnumb_mult(k),
     &                 totcoordcount_mult(k),1,MPI_INTEGER,MPI_SUM,
     &                 MPI_COMM_WORLD, ierr)

                  if (totcoordcount_mult(k).gt.0) then
                     avgxcoordf(k) = totxcoordsum_mult(k) /
     &               DBLE(totcoordcount_mult(k)) !Matt T.
                     avgycoordf(k) = totycoordsum_mult(k) /
     &               DBLE(totcoordcount_mult(k))
                     avgzcoordf(k) = totzcoordsum_mult(k) /
     &               DBLE(totcoordcount_mult(k))
                  endif

               else

                  if (totcoordnumb_mult(k).gt.0) then
                     avgxcoordf(k) = xcoordsum_mult(k) /
     &               DBLE(totcoordnumb_mult(k)) !Matt T.
                     avgycoordf(k) = ycoordsum_mult(k) /
     &               DBLE(totcoordnumb_mult(k))
                     avgzcoordf(k) = zcoordsum_mult(k) /
     &               DBLE(totcoordnumb_mult(k))
                  endif

               endif ! numpe

            enddo ! coalest

c!....Consolidate any average points that are too close to one another
            avgcoordfdist(:,:) = 2.0d0*bubradius
            consol_tag(:) = 0
            avgcoordf_erase_tag(:) = 0

            do k1 = 1, (coalest-1)
               if (avgxcoordf(k1).gt.-1.0d3) then
                  do k2 = (k1+1), coalest
                     avgcoordfdist(k1,k2) = sqrt((avgxcoordf(k1) -
     &               avgxcoordf(k2))**2 + (avgycoordf(k1) -
     &               avgycoordf(k2))**2 + (avgzcoordf(k1) -
     &               avgzcoordf(k2))**2)
                  enddo
               endif
            enddo

            do k1 = 1, (coalest-1)
               do k2 = (k1+1), coalest
                  if (avgcoordfdist(k1,k2).lt.(2.0d0*bubradius)) then
                     do iblk = 1, nelblk
                        iel    = lcblk(1,iblk)
                        npro   = lcblk(1,iblk+1) - iel

                        do i = 1, npro
                           if (coalesc_tag(i,iblk).eq.k2) then
                              coalesc_tag(i,iblk) = k1
                           endif
                        enddo ! i
                     enddo ! iblk

                     consol_tag(k1) = 1
                     avgcoordf_erase_tag(k2) = 1

                  endif ! avgcoordfdist
               enddo ! k2
            enddo ! k1

            do k1 = 1, (coalest-1)
               do k2 = (k1+1), coalest
                  if (avgcoordf_erase_tag(k2).eq.intone) then
                     avgxcoordf(k2) = -1.0d3
                     avgycoordf(k2) = -1.0d3
                     avgzcoordf(k2) = -1.0d3
                  endif
               enddo ! k2

               if (consol_tag(k1).eq.intone) then

                  xcoordsum_mult(k1) = 0.0d0
                  ycoordsum_mult(k1) = 0.0d0
                  zcoordsum_mult(k1) = 0.0d0
                  totcoordnumb_mult(k1) = 0.0d0

                  totxcoordsum_mult(k1) = 0.0d0
                  totycoordsum_mult(k1) = 0.0d0
                  totzcoordsum_mult(k1) = 0.0d0
                  totcoordcount_mult = 0.0d0

                  do iblk = 1, nelblk
                     iel    = lcblk(1,iblk)
                     npro   = lcblk(1,iblk+1) - iel

                     do m = 1, npro
                        if (coalesc_tag(m,iblk).eq.k1) then
                           xcoordsum_mult(k1) = xcoordsum_mult(k1) +
     &                     globalxcoord(m,iblk)
                           ycoordsum_mult(k1) = ycoordsum_mult(k1) +
     &                     globalycoord(m,iblk)
                           zcoordsum_mult(k1) = zcoordsum_mult(k1) +
     &                     globalzcoord(m,iblk)
                           totcoordnumb_mult(k1) = totcoordnumb_mult(k1)
     &                     + totcoordnumbarray(m,iblk)
                        endif
                     enddo !m
                  enddo !iblk

                  if (numpe > 1) then   !Matt T.
                     call MPI_ALLREDUCE(xcoordsum_mult(k1),
     &                    totxcoordsum_mult(k1),1,MPI_DOUBLE_PRECISION,
     &                    MPI_SUM, MPI_COMM_WORLD,ierr)
                     call MPI_ALLREDUCE(ycoordsum_mult(k1),
     &                    totycoordsum_mult(k1),1,MPI_DOUBLE_PRECISION,
     &                    MPI_SUM, MPI_COMM_WORLD,ierr)
                     call MPI_ALLREDUCE(zcoordsum_mult(k1),
     &                    totzcoordsum_mult(k1),1,MPI_DOUBLE_PRECISION,
     &                    MPI_SUM, MPI_COMM_WORLD,ierr)
                     call MPI_ALLREDUCE(totcoordnumb_mult(k1),
     &                    totcoordcount_mult(k1),1,MPI_INTEGER,MPI_SUM,
     &                    MPI_COMM_WORLD, ierr)

                     if (totcoordcount_mult(k1).gt.0) then
                        avgxcoordf(k1) = totxcoordsum_mult(k1) /
     &                  DBLE(totcoordcount_mult(k1))
                        avgycoordf(k1) = totycoordsum_mult(k1) /
     &                  DBLE(totcoordcount_mult(k1))
                        avgzcoordf(k1) = totzcoordsum_mult(k1) /
     &                  DBLE(totcoordcount_mult(k1))
                     endif

                  else

                     if (totcoordnumb_mult(k1).gt.0) then
                        avgxcoordf(k1) = xcoordsum_mult(k1) /
     &                  DBLE(totcoordnumb_mult(k1))
                        avgycoordf(k1) = ycoordsum_mult(k1) /
     &                  DBLE(totcoordnumb_mult(k1))
                        avgzcoordf(k1) = zcoordsum_mult(k1) /
     &                  DBLE(totcoordnumb_mult(k1))
                     endif

                  endif ! numpe
               endif ! consoltag
            enddo ! k1

         else
            avgxcoordf(1) = avgxcoord
            avgycoordf(1) = avgycoord
            avgzcoordf(1) = avgzcoord
         endif !avgxcoord and avgvectdist

         end
