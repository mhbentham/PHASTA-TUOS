! This set of routines help simplify the i/o part of the code:

	subroutine get_path(ipath)
        character*80 ipath
	open(101, file = 'path.dat')
         read(101, 50) ipath
   50    format(A80)
        close(101)
        write(*,*) 'Processing the case located in ', trim(ipath)

	end subroutine

	subroutine get_input
	include "common.f"
! Read the input data:
        open(1, file = trim(ipath)//'/average.inp')

        do i = 1, 6             ! skip the header
          read(1,*)
        end do

        read(1,*) Nrun
        write(*,*) 'Number of runs to process: ', Nrun
        read(1,*)
        read(1,*) Nhom0, nRegions, nRlayers
        write(*,*) 'Number of homogeneous points (coarse region): ', Nhom0
        write(*,*) 'Number of regions: ', nRegions
        write(*,*) 'Number of layers in the dense region: ', nRlayers
        write(*,*) 'Number of homogeneous points (dense region): ', Nhom0*2**(real(nRegions-1))
        read(1,*)
        read(1,*) isym, ifhalf
        write(*,*) 'The symmetric points flag: ', isym, ifhalf
        read(1,*)
        read(1,*) Nzp   ! Number of z planes to process
        write(*,*) 'Number of planes to average: ', Nzp
        read(1,*)
        read(1,*) NStart, Ntime
        write(*,*) 'Timestep number to start, overall number: ', NStart, Ntime
        read(1,*)
        read(1,*) Averdt
        write(*,*) 'Averaging window width: ', Averdt
        read(1,*)
        read(1,*) Nphasedt
        write(*,*) 'NPHASE dt: ', Nphasedt
        read(1,*)
         if (nd2.ge.15) then
                read(1,*) nu_cl, nu_cl2
          write(*,*) 'Fluid kinematic viscosities are: ', nu_cl, nu_cl2
         else
                read(1,*) nu_cl
          write(*,*) 'Fluid kinematic viscosity is: ', nu_cl
         end if
! We have to read in two more parameters (V9.0 update):
! This is the interface half-thickness parameter (case-dependent):
        read(1,*)
        read(1,*) eps_ls
! This is domain length (used for computing flow-throughs, case-dependent):
        read(1,*)
        read(1,*) Lx
        read(1,*)
        read(1,*) Nskip   ! Update Nskip here
        write(*,*) ' Number of points to be skipped each step (Nskip): ', Nskip
        read(1,*)
        read(1,*) nArea   ! Area availability flag
        read(1,*)
        read(1,*) nFixChan   ! Flag to move the second point in the channel flow to the last place
        read(1,*)
        read(1,*) nCurl   ! Print curl(v) request
        if (nCurl.gt.0) write(*,*) ' Curl(U) will be printed '
        close(1)
        if (nd2.ge.15)
     1   write(*,*) 'Assuming the following interface half-thickness : ', eps_ls
        write(*,*) 'Characteristic stream-wise domain length: ', Lx
	end subroutine


