!  PROGRAM: Average
!
!  PURPOSE:   Window-Averaging DNS data provided by PHASTA in for use in NPHASE
! See also Version_info.txt
! Version 12.0: Full (derivatives in all directions) TKE Source terms are now computed (valid not only in BL approximation). 
!		Is necessary for rough wall flows in particular. Igor - 09/26/2011.
!		a) This is done using the "plane" feature, and coordinating the "plane" locations for the proper derivatives can be taken numerically Igor: 11/10/2011.
!		b) Module structure is adopted
! Version 12.1: numvar = 15 is being phased out (19 is fully supported)
! Version 12.2: Added the output for: Mean velocity gradients; turbulent viscosity etc. IAB: 05/11/2012
! Version 12.3: added the non-uniform homogenenuity treatment (e.g. for Pipe), IAB: 04/21/2015
!**********************************************************************************************************

	program Average

	implicit none

        include "vars.f"
	include "common.f"

! read the current case path:
	call get_path(ipath)

	open(102, file = trim(ipath)//'/info.dat')

! xyzts.dat header must be read here !:
        open(2, file = trim(ipath)//'/xyzts.dat')
        read(2, *) np, nskip0, tol, nd1, nd2, nd3
	close(2)

	call get_input

! Read the point data:
	if (nArea.eq.1) then
        open(2, file = trim(ipath)//'/xyz_area.dat')
	else
        open(2, file = trim(ipath)//'/xyzts.dat')
	end if
        read(2, *) np, nskip0, tol, nd1, nd2, nd3

	if (nArea.eq.1) write(*,*) ' Cell area is available'

	if (nFixChan.eq.1) 
     1   write(*,*) ' Data order will be fixed for some channel flows '
! Adjust the Ntime using the Nskip:
	write(*,*) 'Orig Ntime, Nskip, New Ntime: ', Ntime, Nskip, Ntime/Nskip
	Ntime = Ntime/nskip
c	pause
! Allocate the coordinates arrays:
	 allocate(x(np))
         allocate(y(np))
         allocate(z(np))
	if (nArea.eq.1) then
	 allocate(A(np))
         do i = 1, np
           read(2,*) x(i), y(i), z(i), A(i)
         end do
	else
	 do i = 1, np
	   read(2,*) x(i), y(i), z(i)
	 end do
	end if
	close(2)

        call date_and_time(values=time_array_1)
        Global1 = time_array_1 (5) * 3600 + time_array_1 (6) * 60
     1       + time_array_1 (7) + 0.001 * time_array_1 (8)

! Curl should be allocated using homogeneous directions as well:
         if (nCurl.gt.0) then 
	   npnh = np   ! non-hom np
	 end if

        if (Nregions.eq.0) then
           Nhom = Nhom0
        end if

        if (Nregions.eq.0) then
	 if (Nhom.gt.1) then
	   write(*,*) 'Processing ', np, ' points (X',Nhom,' homogeneous) = ', np*Nhom
	 else
	   write(*,*) 'Processing ', np, ' points'
	 end if
        else
! This is the Pipe-type case with variable number of homogeneous directions:
         do iregion = 1, Nregions          
          Nrhom(Nregions - iregion + 1) = int(real(Nhom0)*2.0**(iregion-1))
          Nlhom(iregion) = int(real(NrLayers)*2.0**(iregion-1))
         end do 
         Nhom = Nrhom(1)   ! Definition for the initial time estimate
         MSNlhom = 0    ! Total number of layers
         MSNrhom = 0    ! Total number of hom points
           write(*,*) ' Multi-region homogeneous directions: '
         NLcount(0) = 0
         do iregion = 1, Nregions
           write(*,234) iregion, Nrhom(iregion), Nlhom(iregion)
  234    format(' Region ', I3, '; Number of hom. points: ',I3,
     1    ';  Number of layers: ', I3)
          MSNlhom = MSNlhom + Nlhom(iregion)
          NLcount(iregion) = MSNlhom
          MSNrhom = MSNrhom + Nrhom(iregion)*Nlhom(iregion)
!          write(*,*) ' Layer count = ', NLcount(iregion)
         end do
        write(*,*) 'Total layers: ', MSNlhom, '  total hom. points: ', MSNrhom
        end if

! This is np normalization here: Must be done differently for special cases(Pipe)
        if (Nregions.eq.0) then
           np = int(real(np)/(real(Nhom))/real(Nzp))
        else
          np = int(real(np)/(real(MSNrhom))/real(Nzp)*real(MSNlhom)) 
        end if
        
        write(*,*) 'Normalized np = ', np


	if (Nzp.gt.1) then
	  write(*,*) 'Processing ', Nzp, ' planes of data.'
	end if
	write(*,*) 'nd2 = ', nd2
	if (nd2.eq.14) write(*,*) 'Processing single-phase case '
        if (nd2.eq.15) write(*,*) 'Processing two-phase case '
        if (nd2.eq.19) write(*,*) 'Processing two-phase '
     1       //'case with distance field gradient'

	write(*,*) 'Number of generalized homogeneous directions: ', Nhom

	write(*,*) 'Number of runs to process (ensemble averaging): ', Nrun
	write(*,*) 'It is assumed that all the runs cover the same time range '

! Open the data file:
        reclength = 2*8+3+15*(max(nd2,15))         ! record length
	do ir = 1, Nrun
	  fname = trim(ipath)//'/varts_run'//MyChar2(ir)//'.dat'
          open(19+ir, file = trim(fname)
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
	  write(*,*) 'File #', 19+ir, ' is opened, file name ', trim(fname)
	end do
! We will go through all the timesteps and record the time evolution:
	 Ttime = 0.0E0		! total time of the data
	 allocate(Ctime(0:Ntime))
         allocate(iphase(Ntime,Nhom,Nrun))
	 Ctime(:) = 0.0E0	! current time in the data
	 do j = 1, Ntime
	 if (mod(j,1000).eq.0) write(*,*) 'j, Ntime = ', j, Ntime, np, Nhom, (j-1)*np*Nhom+1
! Recnum should be estimated here:
        if (Nregions.eq.0) then
           Recnum = (j-1)*np*Nhom+1
        else
           Recnum = (j-1)*int(real(np*MSNrhom)/real(MSNlhom))+1
        end if
	 if (nd2.eq.19) then
!            Recnum = Nzp*((j-1)*np*Nhom+(i-1)*Nhom+ih - 1)+mzp
	  read(20, '(2I8, I3, 19E15.7)', REC = Recnum) 
     1     itime(j), jj, iphase(1,1,1), raw1(1:19)
	 else
          read(20, '(2I8, I3, 15E15.7)', REC = Recnum)
     1     itime(j), jj, iphase(1,1,1), raw1(1:15)
	 end if
!         write(*,*) 'u = ', j, raw1(1:3)
	  if (j.gt.0) then 
		Ctime(j) = Ctime(j-1) + raw1(5)
	  else
		Ctime(j) = raw1(5)
	  end if
   	 end do
	 Ttime=Ctime(Ntime)
	 write(*,*) 'Total time of the simulation is : ', Ttime, ' s'
! We have to compute the iphtime and phtime arrays (those are the averaging window locations)
! Estimate the number of averaging windows:
         Nave = nint(Ttime/Nphasedt) + 1
	 allocate(phtime(3,Nave))
         allocate(iphtime(3,Nave))
	 if (nArea.eq.1) then
	   allocate(frate(5,Nave))
	 end if
	 phtime(1, 1) = Ctime(1)
	 phtime(2, 1) = Ctime(1) + Averdt
	 phtime(3, 1) = 0.5*(phtime(1,1)+phtime(2,1))
	 do j = 2, Nave
	   phtime(3,j) = phtime(3,j-1) + Nphasedt
	   phtime(1,j) = phtime(3,j) - 0.5*Averdt
           phtime(2,j) = phtime(3,j) + 0.5*Averdt
	 end do 
   55    continue
! Trim the extra windows:
	 if (phtime(2,Nave).gt.Ttime.and.Nave.gt.1) then
	  write(*,*) 'Removing the last window since ', phtime(2,Nave), ' > ', Ttime
	  Nave = Nave - 1
	  write(*,*) 'Nave = ', Nave
	  goto 55
	 end if
         if (phtime(2,Nave).gt.Ttime.and.Nave.eq.1) then
	  write(*,*) 'Requested window size (', Averdt,
     1    's) exceeds the total simulation time, ', Ttime, 's'
          write(*,*) 
     1    'The only averaing window will cover the whole simulation.'
          phtime(1, 1) = Ctime(1)
          phtime(2, 1) = Ttime
          phtime(3, 1) = 0.5*(phtime(1,1)+phtime(2,1))
	 end if

! Expand the last window to use up the data (optional):
	 phtime(2,Nave) = Ttime
         phtime(3,Nave) = 0.5*(phtime(1,Nave)+phtime(2,Nave))
	 
! Generate the iphtime array:
!$omp   parallel do
	do k = 1, Nave
	do jj = 1, 3
! Go through the timesteps:
	 do j = 2, Ntime
	  if (phtime(jj,k).ge.Ctime(j-1).and.phtime(jj,k).le.Ctime(j)) then
		iphtime(jj,k) = j
		goto 56
	  end if
	 end do ! j
	  write(*,*) ' ******************** WARNING ! ************************'
	  write(*,*) 'There is a problem with ', jj, k, ' phtime location. '
	  write(*,*) 'phtime = ', phtime(jj,k)
	  write(*,*) ' *******************************************************' 
   56    continue
	end do ! jj
	end do ! k

	 write(*,*) 'The following time windows will be processed:'
	 write(102,*) 'The following time windows will be processed:'
	do k = 1, Nave
	 write(*,57) phtime(1,k), phtime(3,k), phtime(2,k), phtime(2,k)-phtime(1,k), 
     1       Nstart+Nskip*iphtime(1,k), 
     1       NStart+Nskip*iphtime(3,k), NStart+Nskip*iphtime(2,k)
     1     , Nskip*(iphtime(2,k) - iphtime(1,k))
         write(102,57) phtime(1,k), phtime(3,k), phtime(2,k), phtime(2,k)-phtime(1,k),
     1       Nstart+Nskip*iphtime(1,k),
     1       NStart+Nskip*iphtime(3,k), NStart+Nskip*iphtime(2,k)
     1     , Nskip*(iphtime(2,k) - iphtime(1,k))
   57    format('Time: ', 4F12.6, '  Steps: ', 4I12)
	end do
        close(102)

! Allocate the main arrays:
         if (Nregions.eq.0) then
            allocate(raw(max(nd2,19),Ntime,2*Nhom,Nrun))
         else
            allocate(raw(max(nd2,19),Ntime,2*Nrhom(1),Nrun))
         end if
         allocate(Mean(19, Nave))   ! The last value will be pressure (which is first in the data)
         allocate(Mean2(19, Nave))
         allocate(Dsym(20, np))
         allocate(Prod(Nave,np))
         allocate(DiffT(Nave,np))
         allocate(DiffP(Nave,np))
         allocate(DiffPL(Nave,np))
         allocate(TKE(2, Nave, np))
         allocate(Tdiff(2, Nave))
         allocate(EPS(2, Nave))
         allocate(Alpha(Nave))
         allocate(Adiff(2, Nave))
         allocate(Stress(6, Nave))
         allocate(TStress(6, Nave))
         if (ncurl.gt.0) allocate(Curl(Nave, 3, npnh))
         if (nd2.eq.19) allocate(WorkI(3, Nave, np))

        do mzp = 1, Nzp   ! Planes loop

	 write(*,*) 'Processing data plane #', mzp
! Open the flow rate data file for the current plane:
	 if (nArea.eq.1) then
	  open(8, file=trim(ipath)//'/flow_rates_'//
     1     '_plane'//MyChar2(mzp)//'.dat')
	 write(8,'(A115)') '   tsn:    time: ''Liquid volume flow rate'': '
     1     //'''Gas volume flow rate'': Total: Alpha:        Ul:'
     2     //'            Ug:      '
	end if
! Loop over averaging windows:
	do k = Nave, 1, -1
!        do k = 1, Nave
! Record the current time:
	  if (nArea.eq.1) then 
              frate(1, k) = phtime(3,k)
              frate(2:5,k) = 0.0E0
          end if
	 write(*,*) 'Window #', k, '   isym = ', isym
	 open(5, file=trim(ipath)//'/inflow_'//MyChar3(k)//
     1     '_plane'//MyChar2(mzp)//'.dat')
         if (Nrun.gt.1) open(6, file=trim(ipath)//'/ensdiff_'//MyChar3(k)//
     1     '_plane'//MyChar2(mzp)//'.dat')
         open(7, file=trim(ipath)//'/stress_'//MyChar3(k)//
     1     '_plane'//MyChar2(mzp)//'.dat')
         if (ncurl.gt.0) open(11, file=trim(ipath)//'/Curl_'//MyChar3(k)//
     1     '_plane'//MyChar2(mzp)//'.dat')
	 write(5,*) ' Number of points, time range, ',
     1     'time instant, number of averaged samples: '
         write(5,58) np, phtime(1:3,k), 
     1     (iphtime(2, k)-iphtime(1, k))*Nhom*Nrun
         if (Nrun.gt.1) write(6,*) ' Number of points, time range, '
     1     ,'time instant, number of averaged samples: '
         if (Nrun.gt.1) write(6,58) np, phtime(1:3,k), 
     1     (iphtime(2, k)-iphtime(1, k))*Nhom*Nrun
         write(7,*) ' Number of points, time range, ',
     1     'time instant, number of averaged samples: '
         write(7,58) np, phtime(1:3,k), 
     1     (iphtime(2, k)-iphtime(1, k))*Nhom*Nrun
         if (ncurl.gt.0) write(11,*) ' Number of points, time range, ',
     1     'time instant, number of averaged samples: '
         if (ncurl.gt.0) write(11,58) npnh/nzp, phtime(1:3,k),   ! Check this!
     1     (iphtime(2, k)-iphtime(1, k))*Nrun
   58    format(I6, F12.7,' : ', F12.7, ';  ', F12.7, I10)
	 if (nd2.ge.15) then
	 write(5,'(A275)') 'x:            y:           z:    '
     1     //'            U:            V:      '
     2     //'        W:           TKE:          EPS:       '
     3     //'   Alpha:          U2:     '
     4     //'      V2:          W2:       TKE2:      EPS2: '
     5     //'   Pressure:     Pressure2:   Tau_xy:   '
     6     //'      dU/dy:      nu_T/nu_cl:      Cmu:  ' 
	 else
         write(5,'(A220)') '          x:               y:              z:      '
     1     //'           U:              V:        '
     2     //'        W:             TKE:            EPS:        '
     3     //'   Pressure:         Tau_xy:     '
     4     //'      dU/dy:         nu_T/nu_cl:         Cmu:  '
	 end if
       if (Nrun.gt.1) 
     1    write(6,'(A96)') 'x:            y:           z:        '
     1       //'    Tdiff:         Adiff:   Liquid FT:    Gas FT: '
         write(7,'(A240)') 'x:            y:           z:          '//
     1     '  u1u1:         u2u2:        u3u3:         u1u2:     '
     2     //'    u1u3:         u2u3:  '//
     2   '      u1''u1'':       u2''u2'':       u3''u3'':     '
     3     //'  u1''u2'':       u1''u3'':      u2''u3'': '
         if (ncurl.gt.0) write(11,'(A101)') 'x:            y:        '//
     1         '         z:   '
     1       //'         Curl_X:          Curl_Y:          Curl_Z:'

! ********* BIG LOOP OVER DIFFERENT POINTS *********************
	i = 0
  901   continue
        if (nFixChan.eq.1.and.i.eq.2) goto 902
  904   continue
	i = i + 1
!	do  i = 1, np
        if (nFixChan.eq.1.and.i.eq.2) goto 904    ! Skipped the second point
  902   continue
! Read the data for the current point in the current window:
	call date_and_time(values=time_array_0)
         if (i.eq.1) write(*,*) 'Expected number of points ',  Nzp*np*Nhom
      	start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 
     1         + time_array_0 (7) + 0.001 * time_array_0 (8)
! Move advanced Nhom here:
! Next loop should have variable size depending on the region. Region is defined
! by "i" value. Can we dynamically change what Nhom is based on "i"? We can!
            if (Nregions.gt.0) then
             do iregion = 1, Nregions
               if (i.gt.NLcount(iregion-1).and.i.le.NLcount(iregion)) then
                        Nhom = Nrhom(iregion)
                  NSpecial = Nhom*(i-1-NLcount(iregion-1))
                  if (iregion.gt.1) then
                   do k1 = 1, iregion-1
                        NSpecial = NSpecial + Nlhom(k1)*Nrhom(k1)
                   end do
                  end if
               end if
             end do
            end if
          do ir = 1, Nrun
            do j = iphtime(1, k), iphtime(2, k)
              do ih = 1, Nhom    !  We need to learn how to vary this for the Pipe cases
! This expression should depend on which region are we in
        if (Nregions.eq.0) then
            Recnum = Nzp*((j-1)*np*Nhom+(i-1)*Nhom+ih - 1)+mzp
        else
           Recnum = Nzp*((j-1)*int(real(np*MSNrhom)/real(MSNlhom))+NSpecial+ih - 1)+mzp
        end if
!         write(*,*) 'recnum, i, j, ih, NS ', Recnum, i, j, ih, NSpecial
         if (nd2.eq.19) then
            read(19+ir, '(2I8, I3, 19E15.7)', REC = Recnum)
     1      itime(j), jj, iphase(j, ih, ir), raw(1:19, j, ih, ir)
         else
            read(19+ir, '(2I8, I3, 15E15.7)', REC = Recnum)
     1      itime(j), jj, iphase(j, ih, ir), raw(1:15, j, ih, ir)
          end if
               end do ! ih
	    end do ! j
         end do ! ir
	call date_and_time(values=time_array_1)
	finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60 
     1       + time_array_1 (7) + 0.001 * time_array_1 (8)
	 write(*,*) ' Read ', i, ' in ', finish_time - start_time, ' out of ', np
! The data is read. Now we need to properly average it:
! First we take care of the mean velocity part.
	Mean(:,k) = 0.0E0
	Mean2(:,k) = 0.0E0
        timeVF1 = 0.0E0  ! Separate time for Volume Fraction computation
        timeVF2 = 0.0E0
	time1 = 0.0E0
	time2 = 0.0E0
        time1p = 0.0E0
        time2p = 0.0E0
	Time_Int1 = 0.0E0
	Time_Int2 = 0.0E0
	Alpha(k) = 0.0E0
!	write(*,*) 1
! single phase distance field correction:
	if (nd2.eq.14) then 
           raw(15,:,:,:) = 1.0   
           raw(19,:,:,:) = 1.0
        end if
        do ih2 = 0, isym - ifhalf  ! These are two components of symmetry in channel flow cases
	 if (mod(Nhom,2).eq.1) then
		Nhom2 = Nhom - 1
	 else
		Nhom2 = Nhom
	 end if 
	do ih = 1, Nhom2, isym+1
	do ir = 1, Nrun
	 do j = iphtime(1, k), iphtime(2, k)
!        write(*,*) ih2, ih, ir, j
! Compute the VF time first (interface separated by the zero level set):
         if (nd2.eq.14) then
                timeVF1 = timeVF1 + raw(5,j,ih+ih2,ir)
               else
          if (raw(19,j,ih+ih2,ir).ge.0.5E0) then
		timeVF1 = timeVF1 + raw(5,j,ih+ih2,ir)
	  else
                timeVF2 = timeVF2 + raw(5,j,ih+ih2,ir)
!                write(*,*) ' Gas detected ', timeVF1,  raw(5,j,ih+ih2,ir)
	  end if
        end if ! nd2
!        write(*,*) 'time is computed'
! This is a more involved approach where we collect the velocities only inside the interface to compute the second phase velocity
!            write(*,*) 'raw19, 5 = ', raw(19,j,ih+ih2,ir), raw(5,j,ih+ih2,ir)
           if (nd2.eq.14) raw(19,j,ih+ih2,ir) = 1.0E0
           if (raw(19,j,ih+ih2,ir).ge.0.99E0) then		! Checking the distance field; This is liquid phase
		time1 = time1 + raw(5,j,ih+ih2,ir)  ! accumulate time
                time1p = time1p + raw(5,iphtime(3,k),ih+ih2,ir)  ! accumulate time for pressure term
!            write(*,*) 'raw19, 5 = ', raw(19,j,ih+ih2,ir), raw(5,j,ih+ih2,ir)
! Stream wise U - velocity is not symmetric: 
!        write(*,*) 'Step 2'
	        Mean(1, k) = Mean(1, k) 
     2          + raw(2, j, ih+ih2,ir)*raw(5, j, ih+ih2,ir)  ! multiplied by dt
!         if (i.eq.4.and.j.eq.iphtime(1, k)+1) write(*,*) ih+ih2, j, raw(2, j,ih+ih2,ir), Mean(1, k)
!         if (i.eq.4) write(*,*) ih+ih2, j, raw(2,j,ih+ih2,ir), Mean(1, k)
! Normal to the wall velocity (V) is NOT symmetric:
!        write(*,*) 'Step 3'
                Mean(2, k) = Mean(2, k) 
     1          + (1)**ih2 * raw(3, j, ih+ih2,ir)*raw(5, j, ih+ih2,ir)  ! multiplied by dt
! Span wise (W) velocity, time, dU/dx, dV/dx and dW/dx is not:
!        write(*,*) 'Step 4'
                Mean(3:7, k) = Mean(3:7, k) 
     1          + raw(4:8, j, ih+ih2,ir)*raw(5, j, ih+ih2,ir)  ! multiplied by dt
! dU/dy, dV/dy and dW/dy is:
!        write(*,*) 'Step 5'
                Mean(8:10, k) = Mean(8:10, k) 
     1          + (-1)**isym * (-1)**ih2 * raw(9:11, j, ih+ih2,ir)
     2            *raw(5, j, ih+ih2,ir)  ! multiplied by dt
! dU/dz, dV/dz and dW/dz is not: 
!        write(*,*) 'Step 6'
                Mean(11:13, k) = Mean(11:13, k) 
     1          + raw(12:14, j, ih+ih2,ir)*raw(5, j, ih+ih2,ir)  ! multiplied by dt
               Mean(14, k) = Mean(14, k) + raw(1, j, ih+ih2,ir)*raw(5, j, ih+ih2,ir)
	   end if
!	   if (raw(15,j,ih+ih2,ir).lt.eps_ls.and.raw(15,j,ih+ih2,ir).gt.-eps_ls) then   ! This is the interface only (use for the dispersed phase velocity computation)
!           if (raw(19,j,ih+ih2,ir).le.1.0E-4) then   ! This is inside the interface only, pure, non mixed gas (use for the dispersed phase velocity computation)
           if (raw(15,j,ih+ih2,ir).le.-eps_ls) then   ! This is inside the interface only, pure, non mixed gas (use for the dispersed phase velocity computation)
                time2 = time2 + raw(5,j,ih+ih2,ir)  ! accumulate time
!            write(*,*) 'raw19, 5 = ', raw(19,j,ih+ih2,ir), raw(5,j,ih+ih2,ir), ih, ih2, j
                time2p = time2p + raw(5,iphtime(3,k),ih+ih2,ir)  ! accumulate time
                Mean2(1:4, k) = Mean2(1:4, k) + raw(2:5, j, ih+ih2,ir)*raw(5, j, ih+ih2,ir)  
                Mean2(5:13, k) = Mean2(5:13, k) + (-1)**ih2 * raw(6:14, j, ih+ih2,ir)*raw(5, j, ih+ih2,ir)  
                Mean2(14, k) = Mean2(14, k) + raw(1, j, ih+ih2,ir)*raw(5, j, ih+ih2,ir)
	   end if
	 end do   ! j
	end do ! ir
	end do ! ih
	end do ! ih2
! Average over time:
	if (time1.gt.0.0E0) Mean(1:14, k) = Mean(1:14, k)/time1
!          if (i.eq.4) write(*,*) 'U(4) = ',Mean(1, k), time1
        if (time2.gt.0.0E0) Mean2(1:14, k) = Mean2(1:14, k)/time2
! Test condition: if U_gas < 5% of U_l, then U_g = U_l:   Igor, August 2012
        if (Mean2(1,k).le.0.05E0*Mean(1,k)) Mean2(1:14,k)= Mean(1:14,k)     
! Compute second phase (negative distance field) volume fraction (Alpha):
	if ((timeVF1+timeVF2).gt.0.0E0) then
		Alpha(k) = timeVF2/(timeVF1+timeVF2)
	else
	  Alpha(k) = -1.0E0 
	end if
!         write(*,*) 'k, time2, alpha: ', k, time2, Alpha(k)
!        write(*,*) 'Step 9', nArea, ncurl
! Here based on the mean velocity and volume fraction we should be able to compute the volumetric
! flow rate (using the cell area information)
	if (nArea.eq.1) then
          ishift = (i-1)*Nhom*Nzp+1
	  frate(2,k) = frate(2,k) + Mean(1,k)*(1.0E0 - Alpha(k))*A(ishift) 
          frate(3,k) = frate(3,k) + Mean2(1,k)*Alpha(k)*A(ishift)
 	  frate(4,k) = frate(4,k) + Alpha(k)*A(ishift)
          frate(5,k) = frate(5,k) + A(ishift)
!	  write(*,*) 'i, A(ishift), frate = ', ishift, A(ishift), frate(5,k)
	end if
! Now it is time to compute the turbulent kinetic energy:
	if (ncurl.gt.0) Curl(k, :, :) = 0.0E0
        TKE(1:2,k,i) = 0.0E0
        Tdiff(1:2,k) = 0.0E0
        Adiff(1:2,k) = 0.0E0
        Stress(1:6,k) = 0.0E0
        TStress(1:6,k) = 0.0E0
        Prod(k,i) = 0.0E0
        DiffT(k,i) = 0.0E0
        DiffP(k,i) = 0.0E0
        DiffPL(k,i) = 0.0E0
        if (nd2.eq.19) WorkI(:,k,i) = 0.0E0
        flprev = 0.0E0
        do ih2 = 0, isym - ifhalf  ! These are two components of symmetry in channel flow cases
        do ih = 1, Nhom2, isym+1
	 do j = iphtime(1, k), iphtime(2, k)
          do ir = 1, Nrun
! Separate the fields at this point:
!        write(*,*) 'Step 11'
	    if (raw(15,j,ih+ih2,ir).gt.5.0*eps_ls.and.time1.gt.0) then	! Liquid phase
		flprev = fluct
		fluct = 0.0E0
		ttfl = 0.E0  ! Turb. Tran. Fluctiation 
!$omp   parallel do
		do k1 = 1, 3
!		    fluct = fluct + 0.5E0*((raw(k1+1, j, ih+ih2, ir)-Mean(k1,k))**2.0E0)*raw(5, j, ih+ih2, ir)
                    fluct = fluct + 0.5E0*(raw(k1+1, j, ih+ih2, ir)-Mean(k1,k))
     1           *(raw(k1+1, j, ih+ih2, ir)-Mean(k1,k))*raw(5, j, ih+ih2, ir)

! Turbulent transport contribution:
		    ttfl = ttfl + 0.5E0*((raw(k1+1, j, ih+ih2, ir)-Mean(k1,k))**2.0E0)
     1                 *(raw(3, j, ih+ih2, ir)-Mean(2,k))*raw(5, j, ih+ih2, ir)
		end do  ! k1, 3	  
!        if (i.lt.2) write(*,*) 'Step 12'
	        TKE(1,k,i) = TKE(1,k,i) + fluct/time1   ! divided by this phase time
	        DiffT(k,i) = DiffT(k,i) - (-1)**ih2 * ttfl/time1	! Only the turbulent transport is added so far.
!        if (i.lt.2) then 
!        write(*,*) k,i,DiffT(k,i)
!        end if
! Add the pressure diffusion component:
             DiffP(k,i) = DiffP(k,i) - 1.0/rho*((-1)**ih2 * raw(3, j, ih+ih2, ir) - Mean(2,k))   ! See the rho definition from above (depends on a case)
     1           *(raw(1, j, ih+ih2, ir) - Mean(14,k))*raw(5, j, ih+ih2, ir)/time1
! This term is equivalent to the pressure diffusion term, computed using "Lahey's" formula: d/dy(u*u'*v') = -1/rho*d/dy(p'v')
             DiffPL(k,i) = DiffPL(k,i) + ((-1)**isym * (-1)**ih2 * raw(3, j, ih+ih2, ir) - Mean(2,k))   ! See the rho definition from above (depends on a case)
     1           *(raw(2, j, ih+ih2, ir) - Mean(1,k))*(raw(2, j, ih+ih2, ir) )*raw(5, j, ih+ih2, ir)/time1
! Ver. 11.8:  Double loop for all components (zero in F.D. smooth wall channel):
	   do i3 = 1, 3  ! Mean velocity component in the mean derivative
	     do j3 = 1, 3  ! Derivative direction component in the mean derivative
	     Prod(k,i) = Prod(k,i) - (raw(i3+1, j, ih+ih2, ir) - Mean(i3,k))
     1           * (raw(j3+1, j, ih+ih2, ir) - Mean(j3,k))*raw(5, j, ih+ih2, ir)/time1    ! Fix the symmetry issue here !!!
     2           * Mean(1+3*j3+i3,k)   ! Should be correct ....         
	     end do
	   end do


! Only liquid phase is used to compute the stress:
	     Stress(1,k) = Stress(1,k) + raw(2, j, ih+ih2, ir)**2.0E0*raw(5, j, ih+ih2, ir)/time1
             Stress(2,k) = Stress(2,k) + raw(3, j, ih+ih2, ir)**2.0E0*raw(5, j, ih+ih2, ir)/time1
             Stress(3,k) = Stress(3,k) + raw(4, j, ih+ih2, ir)**2.0E0*raw(5, j, ih+ih2, ir)/time1
             Stress(4,k) = Stress(4,k) + (-1)**isym * (-1)**ih2 *
     1  raw(2, j, ih+ih2, ir)*raw(3, j, ih+ih2, ir)*raw(5, j, ih+ih2, ir)/time1
             Stress(5,k) = Stress(5,k) + (-1)**isym * (-1)**ih2 *
     1  raw(2, j, ih+ih2, ir)*raw(4, j, ih+ih2, ir)*raw(5, j, ih+ih2, ir)/time1
             Stress(6,k) = Stress(6,k) + (-1)**isym * (-1)**ih2 *
     1  raw(3, j, ih+ih2, ir)*raw(4, j, ih+ih2, ir)*raw(5, j, ih+ih2, ir)/time1
!        write(*,*) 'Step 13'
! Turbulent stress computation:
             TStress(1,k) = TStress(1,k) + (raw(2, j, ih+ih2, ir)-Mean(1,k))
     1           *(raw(2, j, ih+ih2, ir)-Mean(1,k))*raw(5, j, ih+ih2, ir)/time1
             TStress(2,k) = TStress(2,k) + (raw(3, j, ih+ih2, ir)-Mean(2,k))
     1           *(raw(3, j, ih+ih2, ir)-Mean(2,k))*raw(5, j, ih+ih2, ir)/time1
             TStress(3,k) = TStress(3,k) + (raw(4, j, ih+ih2, ir)-Mean(3,k))
     1           *(raw(4, j, ih+ih2, ir)-Mean(3,k))*raw(5, j, ih+ih2, ir)/time1
             TStress(4,k) = TStress(4,k) + (-1)**isym * (-1)**ih2 *(raw(2, j, ih+ih2, ir)-Mean(1,k))
     1           *(raw(3, j, ih+ih2, ir)-Mean(2,k))*raw(5, j, ih+ih2, ir)/time1
             TStress(5,k) = TStress(5,k) + (-1)**isym * (-1)**ih2 *(raw(2, j, ih+ih2, ir)-Mean(1,k))
     1           *(raw(4, j, ih+ih2, ir)-Mean(3,k))*raw(5, j, ih+ih2, ir)/time1
             TStress(6,k) = TStress(6,k) + (-1)**isym * (-1)**ih2 *(raw(3, j, ih+ih2, ir)-Mean(2,k))
     1           *(raw(4, j, ih+ih2, ir)-Mean(3,k))*raw(5, j, ih+ih2, ir)/time1
	     end if
! compute Curl here: ****************************************************
           if (nCurl.gt.0) then
!	write(*,*) 'k,i,ih,ih2=', k, i, ih, ih2 
             Curl(k, 1, (i-1)*Nhom+ih+ih2) = Curl(k, 1, (i-1)*Nhom+ih+ih2) +
     1     (raw(11, j, ih+ih2, ir) - raw(13, j, ih+ih2, ir))*raw(5, j, ih+ih2, ir)/(time1/real(Nhom))
             Curl(k, 2, (i-1)*Nhom+ih+ih2) = Curl(k, 2, (i-1)*Nhom+ih+ih2) +
     1     (raw(12, j, ih+ih2, ir)- raw(8, j, ih+ih2, ir))*raw(5, j, ih+ih2, ir)/(time1/real(Nhom))
             Curl(k, 3, (i-1)*Nhom+ih+ih2) = Curl(k, 3, (i-1)*Nhom+ih+ih2) +
     1     (raw(7, j, ih+ih2, ir) - raw(9, j, ih+ih2, ir))*raw(5, j, ih+ih2, ir)/(time1/real(Nhom))
           end if
!        write(*,*) 'Step 14'
! End of curl computation **********************************************
           if (raw(15,j,ih+ih2,ir).le.-eps_ls.and.time2.gt.0) then ! Gas phase:
                flprev = fluct
                fluct = 0.0E0
                do k1 = 1, 3
                    fluct = fluct + 0.5E0*((raw(k1+1, j, ih+ih2, ir)-Mean2(k1,k))**2.0E0)*raw(5, j, ih+ih2, ir)		! Second field Mean2 is used here!
                end do  ! k1, 3   
		TKE(2,k,i) = TKE(2,k,i) + fluct/time2 
	   end if   ! Gas phase TKE condition
           if (nd2.eq.19) then   ! Interfacial work term computation :
! Part one:  triple correlation term on the liquid side of the interface
            if (raw(19,j,ih+ih2,ir).lt.1.0E0
     1     .and.raw(19,j,ih+ih2,ir).gt.0.5E0) then  
!$omp   parallel do
            do k1 = 1, 3
             do k2 = 1, 3
               WorkI(1, k, i) = WorkI(1, k, i) 
     1       + (raw(1+k1, j, ih+ih2, ir)-Mean(k1,k))
     2        *(raw(1+k2, j, ih+ih2, ir)-Mean(k2,k))
     3        *(raw(1+k1, j, ih+ih2, ir)-Mean(k1,k))
     4        *raw(15+k2, j, ih+ih2, ir)*raw(5, j, ih+ih2, ir)
             end do ! k2
            end do ! k1
              Time_Int1 = Time_Int1 + raw(5, j, ih+ih2, ir)   ! Accumulate the time spent in the liquid part of the interface
            end if   ! part one
! Part two: Gas side of the interface:
            if (raw(19,j,ih+ih2,ir).lt.0.5E0
     1     .and.raw(19,j,ih+ih2,ir).gt.0.0E0) then   
!$omp   parallel do
            do k1 = 1, 3
             do k2 = 1, 3
               WorkI(2, k, i) = WorkI(2, k, i)			! Density is NOT present in those terms
     1       + (raw(1+k1, j, ih+ih2, ir)-Mean2(k1,k))
     2        *(raw(1+k2, j, ih+ih2, ir)-Mean2(k2,k))
     3        *(raw(1+k1, j, ih+ih2, ir)-Mean2(k1,k))
     4        *raw(15+k2, j, ih+ih2, ir)*raw(5, j, ih+ih2, ir)
             end do ! k2
            end do ! k1
              Time_Int2 = Time_Int2 + raw(5, j, ih+ih2, ir)   ! Accumulate the time spent in the liquid part of the interface
            end if   ! part two

	   end if ! nd2.eq.19   Interfacial work term
	  end do   ! ir
! At this point we have recorded the instantenous energies of fluctuations of the last 2 runs.
! Let's compute the non-dimensional difference and record it
	   if (Nrun.gt.1) then
	     CTdiff = abs(fluct - flprev)/(0.5*(fluct + flprev))*raw(5, j, ih+ih2, Nrun)
	     CAdiff = abs(raw(15,j,ih+ih2,Nrun) - raw(15,j,ih+ih2,Nrun-1))*raw(5, j, ih+ih2, Nrun)
!	     if (k.eq.1) write(*,10) CTdiff, fluct, flprev, raw(1+1, j, ih+ih2, 1), raw(1+1, j, ih+ih2, 2)
! Let's accumulate the values:
             Tdiff(1,k) = Tdiff(1,k) + CTdiff/(time1+time2)*2.0
 	     Adiff(1,k) = Adiff(1,k) + CAdiff/(time1+time2)*2.0 
	   end if
	 end do  ! j
	end do  ! ih
	end do  ! ih2
!        write(*,*) 'Step 15'
! Normalize with time:
	if (nd2.eq.19) then
         if (Time_Int1.gt.0.0E0) WorkI(1,k,i) = WorkI(1,k,i)/Time_Int1
         if (Time_Int2.gt.0.0E0) WorkI(2,k,i) = WorkI(2,k,i)/Time_Int2
	end if
! Here we have to average the dissipation rate and substract the mean dissipation from it:
        Eps(1:2,k) = 0.0E0
        do ih2 = 0, isym - ifhalf  ! These are two components of symmetry in channel flow cases
	do ih = 1, Nhom2, isym + 1 
	do ir = 1, Nrun
!$omp   parallel do
         do j = iphtime(1, k), iphtime(2, k)
           if (raw(15,j,ih,ir).gt.eps_ls.and.time1.gt.0) then   ! Liquid phase
                fluct = 0.0E0
                do k1 = 5, 7
                  fluct = fluct + ((raw(k1+1, j, ih+ih2, ir)-Mean(k1,k))**2.0E0)*raw(5, j, ih+ih2, ir)
                end do  ! k1
                do k1 = 8, 10
                  fluct = fluct + (((-1)**isym * (-1)**ih2 * raw(k1+1, j, ih+ih2, ir)-Mean(k1,k))**2.0E0)*raw(5, j, ih+ih2, ir)
                end do  ! k1
                do k1 = 11, 13
                  fluct = fluct + ((raw(k1+1, j, ih+ih2, ir)-Mean(k1,k))**2.0E0)*raw(5, j, ih+ih2, ir)
                end do  ! k1
           Eps(1,k) = Eps(1,k) + nu_cl*fluct/time1
	   end if
           if (raw(15,j,ih,ir).le.-eps_ls.and.time2.gt.0) then   ! Gas phase 
                fluct = 0.0E0
                do k1 = 5, 7
                  fluct = fluct + ((raw(k1+1, j, ih+ih2, ir)-Mean2(k1,k))**2.0E0)*raw(5, j, ih+ih2, ir)
                end do  ! k1
                do k1 = 8, 10
                  fluct = fluct + (((-1)**isym * (-1)**ih2 * raw(k1+1, j, ih+ih2, ir)-Mean2(k1,k))**2.0E0)*raw(5, j, ih+ih2, ir)
                end do  ! k1
                do k1 = 11, 13
                  fluct = fluct + ((raw(k1+1, j, ih+ih2, ir)-Mean2(k1,k))**2.0E0)*raw(5, j, ih+ih2, ir)
                end do  ! k1
	     Eps(2,k) = Eps(2,k) + nu_cl2*fluct/time2
	   end if
         end do   ! j
	end do   ! ir
	end do   ! ih
	end do   ! ih2
!        write(*,*) 'Step 16'
! In case of one phase present we should "equalize" the values:
	 if (Alpha(k).lt.1E-06) then
	  Mean2(1:3,k) = Mean(1:3,k)
	  TKE(2,k,i) = TKE(1,k,i)
	  Eps(2,k) = Eps(1,k)
	 end if

         if (1.0E0-Alpha(k).lt.1E-06) then
          Mean(1:3,k) = Mean2(1:3,k)
          TKE(1,k,i) = TKE(2,k,i)
          Eps(1,k) = Eps(2,k)
         end if
	 if (Eps(2,k).lt.1E-06) Eps(2,k) = 0.0E0

! Compute the number of flow-through for each phase:
	 FT(1, k) = (phtime(2,k)-phtime(1,k))*max(Mean(1,k),Mean(2,k),Mean(3,k))/Lx 
	 FT(2, k) = (phtime(2,k)-phtime(1,k))*max(Mean2(1,k),Mean2(2,k),Mean2(3,k))/Lx

! Print the result in the single file:
! Another correction for the proper coordinates:
          if (nRegions.eq.0) then
	    Mcoord = (i-1)*Nhom*Nzp + mzp
           else
             do iregion = 1, Nregions
               if (i.gt.NLcount(iregion-1).and.i.le.NLcount(iregion)) then
                        Nhom = Nrhom(iregion)
                  NSpecial = Nhom*(i-1-NLcount(iregion-1))
                  if (iregion.gt.1) then
                   do k1 = 1, iregion-1
                        NSpecial = NSpecial + Nlhom(k1)*Nrhom(k1)
                   end do
                  end if
               end if
             end do
            Mcoord = Nzp*NSpecial+mzp
          end if

! Version 12.2 (IAB, 05/2012): Print out the following additional information in
! the "inflow" data file:
!       tau_xy  
!       dU/dy
!       Estimated nu_T
!       C_mu based on k-eps assumption

         tau_xy = dabs(TStress(4,k))   ! Might be the source of the problem: check symmetry in TStress above !
         dUdy = dabs(Mean(8,k))
         if (dabs(dUdy).gt.1.0D-07) then 
                nuT = tau_xy/dUdy 
                Cmu = nuT*Eps(1,k)/TKE(1,k,i)**2.0D0
            else
                nuT = 0.0
                Cmu = 0.0
         end if

	 if (nd2.ge.15) then
	  write(5,10)  x(Mcoord), y(Mcoord), z(Mcoord), Mean(1:3, k), TKE(1,k,i), Eps(1,k)
     1     , Alpha(k), Mean2(1:3, k), TKE(2,k,i), Eps(2,k), Mean(14,k), Mean2(14,k)
     2     , tau_xy, dUdy, nuT/nu_cl, Cmu
	 else
          write(5,10)  x(Mcoord), y(Mcoord), z(Mcoord), Mean(1:3, k), TKE(1,k,i), Eps(1,k), Mean(14,k)
     1     , tau_xy, dUdy, nuT/nu_cl, Cmu
	 end if
	 if (Nrun.gt.1) write(6,13)  x(Mcoord), y(Mcoord), z(Mcoord), Tdiff(1,k), Adiff(1,k), FT(1:2, k)
          write(7,10)  x(Mcoord), y(Mcoord), z(Mcoord), Stress(1:6, k), TStress(1:6, k)
!         write(*,*) 'l. 619'
         if (ncurl.gt.0) then
         do ih2 = 0, isym - ifhalf  ! These are two components of symmetry in channel flow cases
         do ih = 1, Nhom2, isym + 1
!	  write(*,*) 'Printing Curl, i, ih, ih2 = ', i, ih, ih2
         if (abs(x((i-1)*Nhom+ih+ih2)).lt.1.0E-06)   
     1      write(11,10)  x((i-1)*Nhom+ih+ih2), y((i-1)*Nhom+ih+ih2), 
     1      z((i-1)*Nhom+ih+ih2), Curl(k, 1:3, (i-1)*Nhom+ih+ih2)
!            write(*,*) 'Into Curl: ', x((i-1)*Nhom+ih+ih2), y((i-1)*Nhom+ih+ih2),
!     1      z((i-1)*Nhom+ih+ih2), Curl(k, 1:3, (i-1)*Nhom+ih+ih2)
         end do ! ih
         end do ! ih2
         end if
	if (nFixChan.eq.0) then  ! Standart treatment
	 if (i.lt.np) then
	    goto 901
	   else
	    goto 903
	 end if
	else	! Fix channel treatment
	 if (i.eq.np) then 
            i = 2
            goto 902
	 end if
	 if (i.eq.2) goto 903
	 if (i.ne.2.and.i.ne.np) goto 901
	end if 
         write(*,*) 'l. 646'
	
!	end do  ! i, 1, np ********* i loop end
  903    continue  
	 close(5)
	 if (Nrun.gt.1) close(6)
         close(7)
	 if (ncurl.gt.0) close(11)

!	 if (np.eq.50) then    ! CHECK THE NECESSETY OF THIS !!!
!	    isym = 1
!	 else
!	    isym = 0
!	 end if
! Create a separate result file with symmetry averaged points (works for channel):
!	if (isym.eq.1) then
! Read the data file:
         open(5, file=trim(ipath)//'/inflow_'//MyChar3(k)//'_plane'//MyChar2(mzp)//'.dat')
	 read(5,*); read(5,*); read(5,*)  ! Skip the header
	 if (nd2.ge.15) then
	 do i = 1, np
	  read(5,*) Dsym(1:14, i)
 	 end do  ! i loop 
	 else
         do i = 1, np
          read(5,*) Dsym(1:8, i)
         end do ! i loop
	 end if
	 close(5)
! Write the result:
         if (isym.eq.1) 
     1     open(9, file=trim(ipath)//'/symflow_'
     2      //MyChar3(k)//'_plane'//MyChar2(mzp)//'.dat')
         open(10, file=trim(ipath)//'/TKE_RHS__'//MyChar3(k)//'_plane'//MyChar2(mzp)//'.dat')
         if (isym.eq.1) write(9,*) 
     1 '*SYMMETRIC DATA* N of points, t range, t inst, n of ave sampl: '
         if (isym.eq.1) 
     1     write(9,58) np, phtime(1:3,k), (iphtime(2, k)-iphtime(1, k))*Nhom*Nrun*2
         write(10,*) 
     1 '*TKE RHS terms* N of points, t range, t inst, n of ave sampl: '
         write(10,58) np, phtime(1:3,k), (iphtime(2, k)-iphtime(1, k))*Nhom*Nrun*2
        if (isym.eq.1) then 
        if (nd2.ge.15) then
         write(9,'(A190)') '          x:               y:              z:              U:              V:      '
     1     //'          W:             TKE:            EPS:            Alpha:          U2:     '
     2     //'      V2:          W2:       TKE2:      EPS2: '
         else
         write(9,'(A110)') 
     1'x:            y:           z:            U:            V:   '
     1     //'           W:           TKE:          EPS:'
         end if
        end if
	 if (nd2.eq.19) then
         write(10,'(A169)') 
     1 'x:            y:           z:          Prod:            Diss: '
     1     //'         DiffV:     DiffT:        DiffP:   '
     2  //'    Total:      DiffPL:    Work_1:     Work_2:    Total_W: '
	 else
         write(10,'(A134)')
     1 'x:            y:           z:          Prod:            Diss: '
     1     //'         DiffV:     DiffT:        DiffP:   '
     2     //'    Total:      DiffPL: '
	 end if
!$omp   parallel do
	 do i = 1, np
! Average the data before the print out:
	  i1 = i
	  i2 = np - i + 1
	  Csym(1:3) = Dsym(1:3, i1)
	  Csym(4:14) = 0.5D0*(Dsym(4:14,i1) + Dsym(4:14,i2))
         if (isym.eq.1) then
         if (nd2.ge.15) then
          write(9,10)  Csym(1:14)
         else
          write(9,10)  Csym(1:8)
         end if
         end if ! isym
!   Compute the diffusion derivatives:
! This will be more complicated as we will need to compute the derivatives in stream-wise and span-wise directions, which are not available from the neighboring points approach !!!
	if (i.eq.1) then
         DDiffV = (TKE(1,k,i+2)-TKE(1,k,i+1))/(Dsym(2,i+2)-Dsym(2,i+1))
         DDiffV = DDiffV - (TKE(1,k,i+1)-TKE(1,k,i))/(Dsym(2,i+1)-Dsym(2,i))
         DDiffV = 2.0E0*nu_cl*DDiffV / (Dsym(2,i+2)-Dsym(2,i))
         DDiffP = (DiffP(k,i+1)-DiffP(k,i))/(Dsym(2,i+1)-Dsym(2,i))
         DDiffPL = (DiffPL(k,i+1)-DiffPL(k,i))/(Dsym(2,i+1)-Dsym(2,i))
         DDiffT = (DiffT(k,i+1)-DiffT(k,i))/(Dsym(2,i+1)-Dsym(2,i))
	else if (i.eq.np) then
         DDiffV = (TKE(1,k,i)-TKE(1,k,i-1))/(Dsym(2,i)-Dsym(2,i-1))
         DDiffV = DDiffV - (TKE(1,k,i-1)-TKE(1,k,i-2))/(Dsym(2,i-1)-Dsym(2,i-2))
         DDiffV = 2.0E0*nu_cl*DDiffV / (Dsym(2,i)-Dsym(2,i-2))
         DDiffT = (DiffT(k,i)-DiffT(k,i-1))/(Dsym(2,i)-Dsym(2,i-1))
         DDiffP = (DiffP(k,i)-DiffP(k,i-1))/(Dsym(2,i)-Dsym(2,i-1))
         DDiffPL = (DiffPL(k,i)-DiffPL(k,i-1))/(Dsym(2,i)-Dsym(2,i-1))
	else
         DDiffV = (TKE(1,k,i+1)-TKE(1,k,i))/(Dsym(2,i+1)-Dsym(2,i))
	 DDiffV = DDiffV - (TKE(1,k,i)-TKE(1,k,i-1))/(Dsym(2,i)-Dsym(2,i-1))
	 DDiffV = 2.0E0*nu_cl*DDiffV / (Dsym(2,i+1)-Dsym(2,i-1))
         DDiffT = (DiffT(k,i+1)-DiffT(k,i-1))/(Dsym(2,i+1)-Dsym(2,i-1))
         DDiffP = (DiffP(k,i+1)-DiffP(k,i-1))/(Dsym(2,i+1)-Dsym(2,i-1))
         DDiffPL = (DiffPL(k,i+1)-DiffPL(k,i-1))/(Dsym(2,i+1)-Dsym(2,i-1))
	end if
	 if (nd2.ne.19) then
	  write(10,10) Csym(1:3), Prod(k,i), -Dsym(8,i), DDiffV, DDiffT, DDiffP
     1        , Prod(k,i) - Dsym(8,i) + DDiffV + DDiffT + DDiffP, DDiffPL
	 else
         write(10,10) Csym(1:3), Prod(k,i), -Dsym(8,i), DDiffV, DDiffT, DDiffP
     1        , Prod(k,i) - Dsym(8,i) + DDiffV + DDiffT + DDiffP, DDiffPL
     2        , WorkI(1:2, k, i), WorkI(1,k,i)+WorkI(2,k,i)
	 end if ! nd2
         end do   ! i loop
	if (isym.eq.1)  close(9)
	 close(10)

! Write the current window flow rAte: 
	if (nArea.eq.1) then
         write(8,11)  NStart+Nskip*iphtime(2,k), 
     1   frate(1:3,k), frate(2,k)+frate(3,k), frate(4,k)/frate(5,k)
     2    , frate(2,k)/frate(5,k), frate(3,k)/frate(5,k)
	 if (k.eq.1) write(*,*) ' Inflow Area is ', frate(5,k)
	end if
	end do ! k, 1, Nave

	if (nArea.eq.1)  close(8)  ! frate file is done here

        end do  ! mzp, Nzp - the end of different planes loop

          close(20)

! Deallocate the arrays:
	 deallocate(Mean)
         deallocate(Mean2)
         deallocate(TKE)
         deallocate(Dsym)
         deallocate(Prod)
         deallocate(DiffT)
         deallocate(DiffP)
         deallocate(DiffPL)
         deallocate(EPS)
         deallocate(Alpha)
         deallocate(raw)
         deallocate(iphase)
         deallocate(Ctime)
         deallocate(phtime)
         deallocate(iphtime)
	 deallocate(Tdiff)
         deallocate(Adiff)
         deallocate(x)
         deallocate(y)
         deallocate(z)
	if (nArea.eq.1) deallocate(A)
	if (nd2.eq.19) deallocate(WorkI)
	if (ncurl.gt.0) deallocate(Curl)
        call date_and_time(values=time_array_1)
        Global2 = time_array_1 (5) * 3600 + time_array_1 (6) * 60
     1       + time_array_1 (7) + 0.001 * time_array_1 (8)

	write(*,*) 'The DNS data has been processed in ', Global2-Global1, ' s'
	goto 2
1 	continue
	write(*,*)
	write(*,*) 'Program terminated. Check the error message above.'
	write(*,*)
	write(*,*)
2	continue

10	format(1x, 30E17.6)
11	format(1x, I7, 10E15.6)
12      format(1x, 3I6, 20E20.10)
13	format(1x, 5E14.6, 2F12.5)

	end program Average

	include "char_func.f"

	include "io.f"
