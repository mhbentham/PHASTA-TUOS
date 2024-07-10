! Version 3.0:
!   Has the ability to skip some points to reduce the final file size
! Version 4.0: Has the option of shifting the pressure at each time step to have the average = 0
! Version 4.1: Accepts 19-variable files (which has gradient of the level set function recorded) 
! Igor A. Bolotnov, 11/08/2010.

	program Merge1 

	implicit none

	character*1 MyChar1
	external MyChar1
        character*2 MyChar2
        external MyChar2
        character*3 MyChar3
        external MyChar3
        character*4 MyChar4
        external MyChar4
        character*5 MyChar5
        external MyChar5
        character*6 MyChar6
        external MyChar6

	character*80 ipath

	real tol, varts(1:20), Saved_dt, pressure

	integer*8 i1, i, j, k
	integer*8 Nrun, Nfiles, np, nskip, nd1, nd2 
	integer*8 reclength, lstep, jj, icount, iphase, ipres
	integer*8 istep(1:180), istart, istop
	integer*8 fout, fin, isn, iread, NewSkip, RecordSkip

! read the current case path:
	open(101, file = '../Latest/path.dat')
	 read(101, 50) ipath
   50    format(A80)
	close(101)

	write(*,*) 'Processing the case located in ', trim(ipath)

! Read the input data:
	open(1, file = trim(ipath)//'/merge.inp')

	do i = 1, 6		! skip the header
	  read(1,*)
	end do

	read(1,*) Nrun
	read(1,*)
        read(1,*) istart		! Range of timesamples in the result
	read(1,*) istop
        read(1,*)			! Number of steps to skip
        read(1,*) NewSkip
        read(1,*)
        read(1,*) Nfiles
	read(1,*)
	do i = 1, Nfiles+1
	 read(1,*) istep(i)
	end do
        read(1,*)
        read(1,*) ipres

	close(1)

        write(*,*) 'timestep range : ', istart, istop

        if (ipres.eq.1) then
          write(*,*) ' Pressure in each time step will be adjusted to result in zero averaged pressure'
          open(15, file = trim(ipath)//'/pressure_run'//MyChar2(Nrun)//'.dat')
        end if
! Read the point data:
	open(2, file = trim(ipath)//'/xyzts.dat')
	
	read(2, *) np, nskip, tol, nd1, nd2
	close(2)

	write(*,*) 'Number of points per step = ', np

	write(*,*) 'PHASTA data has the following step skip parameter (nskip): ', nskip
	write(*,*) 'The merged file may be reduced using the following step skip (NewSkip): ', NewSkip

	 if (nskip .ge. NewSkip) then
	  RecordSkip = 1
	 write(*,*) ' All the records are transferred since the PHASTA skip is >= than the requested skip'
	 end if

	 if (nskip .lt. NewSkip) then
	  RecordSkip = int(NewSkip/nskip)
	  write(*,*) ' Every ', RecordSkip, ' step will be transfered'
	 end if
	write(*,*) 'Processing run #', Nrun

	if (nd2.eq.19) then
	  reclength = 2*8+3+15*19          ! record length
	 else
          reclength = 2*8+3+15*15      
	end if
          open(20, file = trim(ipath)//'/varts_run'//MyChar2(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')

! loop over files:
	icount = 0
	do j = 1, Nfiles
	 write(*,*) 'Processing the step number: ', istep(j) 
	if (Nrun.gt.9) then
        if (istep(j).lt.1000) then
          open(3, file = trim(ipath)//'/varts.'//MyChar3(istep(j))//'.run.'//MyChar2(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
	 else if (istep(j).lt.10000) then
          open(3, file = trim(ipath)//'/varts.'//MyChar4(istep(j))//'.run.'//MyChar2(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
	 else if (istep(j).lt.100000) then
          open(3, file = trim(ipath)//'/varts.'//MyChar5(istep(j))//'.run.'//MyChar2(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
	 else
          open(3, file = trim(ipath)//'/varts.'//MyChar6(istep(j))//'.run.'//MyChar2(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
	 end if
	else
        if (istep(j).lt.1000) then
          open(3, file = trim(ipath)//'/varts.'//MyChar3(istep(j))//'.run.'//MyChar1(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
         else if (istep(j).lt.10000) then
          open(3, file = trim(ipath)//'/varts.'//MyChar4(istep(j))//'.run.'//MyChar1(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
         else if (istep(j).lt.100000) then
          open(3, file = trim(ipath)//'/varts.'//MyChar5(istep(j))//'.run.'//MyChar1(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
         else
          open(3, file = trim(ipath)//'/varts.'//MyChar6(istep(j))//'.run.'//MyChar1(Nrun)//'.dat'
     1     , status='unknown', form='formatted', recl=reclength, access='direct')
         end if
	end if

! read/write the information
         iread = 0
	 do i = 1, istep(j+1)-istep(j)   ! Loop over number of steps   !/nskip
          isn = istep(j) + i - 1
! We check if the step is dividable by Nskip:
	  fin = 0
	  if (mod(isn,NSkip).eq.0) fin = 1	  
!          if (isn.eq.105199) fin = 0
          if (isn.eq.58999) fin = 0
!          if (isn.eq.19996) fin = 0
!          if (isn.eq.19499) fin = 0
!          if (isn.eq.69899) fin = 0
!          if (isn.eq.69898) fin = 0
	  if (fin.gt.0) then
!	 write(*,*) 'start reading time step number : ', isn
	  iread = int((NSkip - 1 + i)/NSkip)
!	 write(*,*) 'iread, i, Nskip = ', iread, i, Nskip

! If required, read all the data for the sole purpose of obtaining the averaged pressure:
          if (ipres.eq.1) then
            pressure = 0.0E0
           do i1 = 1, np
	 if (nd2.eq.19) then
          read(3, '(2I8, I3, 19E15.7)', REC=np*(iread-1)+i1)
     1            lstep,jj,iphase,(varts(k), k=1, 19)
 	 else
          read(3, '(2I8, I3, 15E15.7)', REC=np*(iread-1)+i1)
     1            lstep,jj,iphase,(varts(k), k=1, 15)
	 end if
           pressure = pressure + varts(1)/real(np)
           end do    ! i1 loop
           write(15, 10) lstep, pressure
	  end if   ! ipress
	   do i1 = 1, np
         if (nd2.eq.19) then
          read(3, '(2I8, I3, 19E15.7)', REC=np*(iread-1)+i1)
     1            lstep,jj,iphase,(varts(k), k=1, 19)
         else
          read(3, '(2I8, I3, 15E15.7)', REC=np*(iread-1)+i1)
     1            lstep,jj,iphase,(varts(k), k=1, 15)
         end if
	  Saved_dt = varts(5)
	  fout = 0
	  if (lstep.ge.istart.and.lstep.lt.istop) fout = 1
!           if (i1.eq.1) write(*,*) 'istart = ', istart, fout, lstep
!	  if (i.eq.100) write(*,*) varts(1:15)
! Here we will skip some steps of dt is small enough (set up for two-phase jet cases):
	   if (fout.eq.1.and.mod(isn, RecordSkip).eq.0) then
		varts(5) = real(RecordSkip)*varts(5)
		fout = 1
	   else  
		fout = 0
	   end if
          if (i1.eq.np.and.fout.eq.1) write(*,20) i, lstep, isn, jj, iphase, fout
	  if (fout.gt.0) then
	  icount = icount + 1 
          if (nd2.eq.19) then
          if (ipres.eq.1) then
	    write(20,'(2I8, I3, 19E15.7)', REC=icount)
     1            lstep,jj,iphase, varts(1) - pressure, (varts(k), k=2, 19) 
            else
            write(20,'(2I8, I3, 19E15.7)', REC=icount)
     1            lstep,jj,iphase, (varts(k), k=1, 19)
           end if   ! ipress
	  else
          if (ipres.eq.1) then
            write(20,'(2I8, I3, 15E15.7)', REC=icount)
     1            lstep,jj,iphase, varts(1) - pressure, (varts(k), k=2, 15)
            else
            write(20,'(2I8, I3, 15E15.7)', REC=icount)
     1            lstep,jj,iphase, (varts(k), k=1, 15)
           end if   ! ipress
	  end if ! nd2.eq.19
	  end if  ! if fout

	  end do ! i1, np

	  end if ! fin
	 end do   ! i, isteps
	 close(3)
	end do  ! j, Nfiles

	close(20)

	if (ipres.eq.1) close(15)
   10   format(1x, I10, 5E15.7)
   20   format(1x, 20I10)
        write(*,*) ' Merge has been successful! '
	end program Merge1 

	character*(1) function MYCHAR1(i)
	integer i
	 MYCHAR1(1:1) = ACHAR(ICHAR('0')+i)
	end function

	character*(2) function MYCHAR2(i)
	integer i
	 MYCHAR2(1:1) = ACHAR(ICHAR('0')+INT(i/10))
	 MYCHAR2(2:2) = ACHAR(ICHAR('0')+i-10*INT(i/10))
	end function

	character*(3) function MYCHAR3(i)
	integer i
	 MYCHAR3(1:1) = ACHAR(ICHAR('0')+INT(i/100))
	 MYCHAR3(2:2) = ACHAR(ICHAR('0')+INT(i/10)-10*INT(i/100))
	 MYCHAR3(3:3) = ACHAR(ICHAR('0')+i-10*INT(i/10))
	end function

	character*(4) function MYCHAR4(i)
	integer i
	 MYCHAR4(1:1) = ACHAR(ICHAR('0')+INT(i/1000))
	 MYCHAR4(2:2) = ACHAR(ICHAR('0')+INT(i/100)-10*INT(i/1000))
	 MYCHAR4(3:3) = ACHAR(ICHAR('0')+INT(i/10)-10*INT(i/100))
	 MYCHAR4(4:4) = ACHAR(ICHAR('0')+i-10*INT(i/10))
	end function

	character*(5) function MYCHAR5(i)
	integer i
	 MYCHAR5(1:1) = ACHAR(ICHAR('0')+INT(i/10000))
	 MYCHAR5(2:2) = ACHAR(ICHAR('0')+INT(i/1000)-10*INT(i/10000))
	 MYCHAR5(3:3) = ACHAR(ICHAR('0')+INT(i/100)-10*INT(i/1000))
	 MYCHAR5(4:4) = ACHAR(ICHAR('0')+INT(i/10)-10*INT(i/100))
	 MYCHAR5(5:5) = ACHAR(ICHAR('0')+i-10*INT(i/10))
	end function

	character*(6) function MYCHAR6(i)
	integer i
	 MYCHAR6(1:1) = ACHAR(ICHAR('0')+INT(i/100000))
	 MYCHAR6(2:2) = ACHAR(ICHAR('0')+INT(i/10000)-10*INT(i/100000))
	 MYCHAR6(3:3) = ACHAR(ICHAR('0')+INT(i/1000)-10*INT(i/10000))
	 MYCHAR6(4:4) = ACHAR(ICHAR('0')+INT(i/100)-10*INT(i/1000))
	 MYCHAR6(5:5) = ACHAR(ICHAR('0')+INT(i/10)-10*INT(i/100))
	 MYCHAR6(6:6) = ACHAR(ICHAR('0')+i-10*INT(i/10))
	end function
