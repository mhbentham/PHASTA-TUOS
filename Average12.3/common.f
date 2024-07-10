! Common blocks and the corresponding variable declarations:

      character*80 ipath      	! Path to the case location
      real*8 nu_cl, nu_cl2,    	! Liquid and gas kinematic viscosity
     1		eps_ls,		! Level set interface thickness (used to distinguish the "pure" liquid and gas)
     2		Lx,		! Channel length used for the computation of "flow-throughs" for estimating convergence criteria    
     3		Averdt,		! Averaging time window width
     4		Nphasedt	! Time between the centers of two sequential time windows (does not have to be equal to Averdt)

      integer*8 	Nrun,   		! Number of runs to process
     1        	Nhom0,		! Number of homogeneous points
     2        	isym,		! Presence of symmetrical points flag
     3		ifhalf,		! Flag to process the second half of the symmetrical points
     4		Nzp,		! Number of data planes available (sets of distinct probe locations)
     5		Nstart,		! First time step number (for reference purposes only)
     5		Ntime,		! Number of timesteps to be processed
     6		nd2, nd22,	! Number of variable in the varts data files (14,15,19, more in the future)
     7		Nskip,		! Number of points to skip (saves time and more efficient if Solver CFD < 1.0)
     8		nArea,		! Availability of the xyz_area.dat data file which allows to compute the flow rates
     9		nFixChan,	! Flag for "fixing" data point sequence in some runs (should not be needed long term)
     *		nCurl,		! Flag for expensive computation of Curl vectors (not sure if finished) 
     1          nRegions,       ! number of regions in non-uniform homogeneous setup
     2          nRlayers        ! number of layers in the first (most dense) region

      common /input/ ipath, Nrun, Nhom0, isym, ifhalf, Nzp, Nstart,
     1               Ntime, Averdt, Nphasedt, nd2, nd22, nu_cl, nu_cl2,
     2               eps_ls, Lx, Nskip, nArea, nFixChan, nCurl,
     3               nRegions, nRlayers
     
