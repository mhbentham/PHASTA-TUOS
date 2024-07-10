!        module cf_arrays
!!----------------------------------------------------------------------
!!
!!       This module contains some important arrays and variables, which
!!       are not convinent to transfer between subroutines.
!!       Jun Fang,                                             Fall,2014
!!
!!----------------------------------------------------------------------
!        integer iyhistst, iyhisten
!        real*8, allocatable ::  cf_var(:,:)
!        real*8, allocatable ::  xcf(:), xcf_old(:)
!        real*8, allocatable ::  ycf(:), ycf_old(:)
!        real*8, allocatable ::  zcf(:), zcf_old(:)
!
!        end module
!
        module bub_track
!----------------------------------------------------------------------
!
!       This module bridges the bubble information array among different
!       subroutines
!       Jun Fang,                                           Summer,2015
!
!----------------------------------------------------------------------
        integer :: ts_hold
        real*8  XLEN, YLEN, ZLEN
        real*8  Xmid, Ymid, Zmid
        real*8  NbrhdRatio      !determine the neighborhood of a certain
                                !bubble 
        real*8  GhostRatio      !determine the ghost bubble regions 
        real*8, allocatable ::  bub_cent(:,:)
        real*8, allocatable ::  bub_info(:,:)
        real*8, allocatable ::  avg_info(:,:)

        real*8, allocatable ::  procs_dataset(:,:)
        real*8, allocatable ::  procs_coordDn(:,:)
        real*8, allocatable ::  procs_coordUp(:,:)
        real*8, allocatable ::  Shear_NodeMin(:,:)
        real*8, allocatable ::  Shear_NodeMax(:,:)

        real*8, allocatable ::  unive_dataset(:,:)
        real*8, allocatable ::  unive_coordDn(:,:)
        real*8, allocatable ::  unive_coordUp(:,:)
        real*8, allocatable ::  bubbl_coordDf(:,:)
        real*8, allocatable ::  Shear(:,:)
!       the following variable and array are for bubble tracking based
!       coalenscence control
        integer ncoalEvent
        real*8, allocatable ::  coalCenter(:,:)
!       the following array is for bubble breakup recognition capability
        real*8, allocatable ::  breakupSeeder(:,:)

        end module
