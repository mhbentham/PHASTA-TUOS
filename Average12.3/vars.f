! VARIABLES MODULE

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

        character*100 fname

        integer, parameter :: mts = 130000
        real*8, parameter :: pi = 3.14159265D0
        real*8, parameter :: rho = 1.00                 ! Consistent with channel flow data sets only !!
        real*8, parameter :: rho_g = 0.001165

        integer i, j, k, m, kk, k1, k2

        integer Nwin, its, its2
        integer np, nd1

! Added in version 12.3:
        integer MSNlhom, MSNrhom, NSpecial
        integer  NLcount(0:10)

        real*8 tol
        real*8 FT(2, mts)
        real*8 Csym(1:20)

! Allocatable arrays are here:
        real*8, allocatable :: x(:), y(:), z(:), A(:)

        real*8, allocatable :: raw(:,:,:,:), Mean(:,:), Mean2(:,:)
     1   , TKE(:,:,:), Eps(:,:), phtime(:,:), Ctime(:), Alpha(:)
     2   , Tdiff(:,:), Adiff(:,:), Stress(:,:), TStress(:,:)
     3   , Prod(:,:), DiffP(:,:), DiffT(:,:), DiffPL(:,:)   ! The latest is the diffusion term
     4   , WorkI(:,:,:), Curl(:,:,:), Dsym(:,:)

        real*8, allocatable :: frate(:,:)

        integer, allocatable :: iphase(:,:,:), iphtime(:,:)

        integer itime(1:mts), Nave, Nlow, ip1, itime1
        integer Nrhom(1:10), Nlhom(1:10) ! Ver. 12.3 addition
        real*8 raw1(1:20), DDiffV
        real*8 fluct, flprev, CTdiff, CAdiff, ttfl, DDiffT, DDiffP, DDiffPL

        integer reclength, nd3, jj
        real*8 Ttime, time1, time2, timeVF1, timeVF2
        real*8 time1p, time2p, Time_Int1, Time_Int2
        integer Nhom2, ih, ih2, ir, i1, i2, mzp, ishift, Nhom
        integer iregion  ! Ver. 12.3 addition
        integer i3, j3 ! Ver 11.8 addition
        integer Recnum, Mcoord, nskip0 ! the last one is the redundant one
        integer npnh

! code timing variables:
        integer time_array_0(8), time_array_1(8)
        real*8 start_time, finish_time, Global1, Global2

! Version 12.2: additional vars into Inflow output:
        real*8 tau_xy, Cmu, dUdy, nuT




