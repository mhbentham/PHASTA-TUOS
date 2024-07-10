!!*****************************************************
!
!         Turbo AMG (ramg) interface functions
!         1. V_cycle
!         2. plain
!         3. driver      
!
!*****************************************************

      !*******************************************
      !    ramg_V_cycle:
      !     One V_cycle call. 
      !     Now use SAMG, later will use own code
      !*******************************************
      recursive subroutine ramg_V_cycle(
     &               ramg_sol,ramg_rhs,
     &               ramg_res_i,ramg_res_o,
     &               ramg_mode,clevel,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper
     &           )!ramg_V_cycle
      use ramg_data
      include "common.h"
      include "mpif.h"

      !Variable Declaration
      ! arguments
      integer, intent(in)         :: ramg_mode,clevel
      real(kind=8),intent(inout),dimension(amg_nshg(clevel))
     &               :: ramg_sol,ramg_rhs
      real(kind=8),intent(inout)    ::ramg_res_i,ramg_res_o

      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)  :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)  :: lhsP
      integer,intent(in),dimension(nlwork)          :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC


      ! local
      real(kind=8)                          :: dummy_eps;
      real(kind=8) :: tlen,pos,ramg_r_i
      integer::i,j,k
      character*10                                     :: fname
      real(kind=8),dimension(amg_nshg(clevel))  :: myvF,myvE
      real(kind=8),dimension(:),allocatable :: myvC,myvCS
      real(kind=8)          :: cpusec(10)

      ! In scale
      
      if (clevel.ne.1) then
      ramg_rhs = ramg_rhs*amg_ppeDiag(clevel)%p
      endif

      !*********************
      ! SAMG mode
      !*********************
      if (ramg_mode .lt. 10) then
          ! do nothing

      !*********************
      ! finest solver
      !*********************
      else if (ramg_mode+ramg_levelx-clevel .eq. 10) then 
      ! finest level solver
          call ramg_direct_solve(amg_A_colm(clevel)%p,
     &    amg_A_rowp(clevel)%p,amg_A_lhs(clevel)%p,amg_nshg(clevel),
     &    amg_nnz(clevel),ramg_rhs,ramg_sol,
     &    ilwork,BC,iBC,iper)
      !**********************
      ! higher level smoother
      !**********************
      else if (ramg_mode+ramg_levelx-clevel .gt. 10) then
          allocate(myvC(amg_nshg(clevel+1)))
          allocate(myvCS(amg_nshg(clevel+1)))
          ! smoothing on fine
          call ramg_smoother(clevel,myvF,-ramg_rhs,-ramg_rhs,
     &              ramg_r_i,ramg_res_o,colm,rowp,lhsK,lhsP,
     &              ilwork,BC,iBC,iper,2) !dx1
          ramg_sol =  - myvF !x2
          ! restriction
          call ramg_calcAv_g(clevel,myvF,ramg_sol,colm,rowp,lhsK,lhsP,
     &               ilwork,BC,iBC,iper,1)
          myvF = myvF - ramg_rhs ! r2
          call ramg_calcIvFC(clevel,clevel+1,myvF,myvC)
          ! up one level
          myvCS = 0
          call ramg_V_cycle(myvCS,myvC,
     &                      ramg_r_i,ramg_res_o,ramg_mode,
     &                      clevel+1,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
          ! prolongation
          call ramg_calcIvCF(clevel,clevel+1,myvCS,myvF,
     &         ilwork,BC,iBC,iper) ! v2hat
          ! smoothing on fine
          ramg_sol = ramg_sol - myvF
         call ramg_calcAv_g(clevel,myvF,ramg_sol,colm,rowp,lhsK,lhsP,
     &               ilwork,BC,iBC,iper,1)
          myvE = myvF - ramg_rhs ! r3
          call ramg_smoother(clevel,myvF,myvE,myvE,
     &              ramg_r_i,ramg_res_o,colm,rowp,lhsK,lhsP,
     &              ilwork,BC,iBC,iper,1)
          ramg_sol = ramg_sol - myvF
          deallocate(myVC)
          deallocate(myVCS)
      end if
      
      ! Out scale
      if (clevel.ne.1) then
      ramg_sol = ramg_sol * amg_ppeDiag(clevel)%p
      endif

      return

      end subroutine ramg_V_cycle
      

!*****************************************************
!      ramg Cg solver
!      using ramg_V_cycle as preconditioner
!*****************************************************
      subroutine ramg_CG(
     &              ramg_sol,ramg_rhs,mrhs,run_mode,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper
     &           )!ramg_CG
      use ramg_data
      include "common.h"
      include "mpif.h"
      
      real(kind=8),intent(inout),dimension(amg_nshg(1)) :: ramg_sol
      real(kind=8),intent(in),dimension(amg_nshg(1))    :: ramg_rhs
      real(kind=8),intent(in),dimension(nshg,nflow)     :: mrhs
      integer,intent(in)                                :: run_mode

      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)  :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)  :: lhsP
      integer,intent(in),dimension(nlwork)          :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC


      ! local
      real(kind=8),dimension(amg_nshg(1))   :: cgP,cgQ,cgZ,cgR
      real(kind=8)                       :: rz,rz_0,pq,alpha,beta
      real(kind=8)               :: tmp,norm_0,norm_1,norm_e,norm_c
      real(kind=8)               :: mres_0,mres_n
      real(kind=8),dimension(nshg,4)               :: diag
      integer                    :: iterNum
 

      !call drvLesPrepDiag(diag,ilwork,iBC,BC,iper,rowp,colm,lhsK,lhsP)
      diag = ramg_flowDiag%p
     
      !cgR = ramg_rhs 
      call ramg_PPErhs(cgR,mrhs,diag,colm,rowp,lhsK,lhsP,
     &                 ilwork,BC,iBC,iper)
!      call ramg_L2_norm(amg_nshg(1),cgR,norm_0)
      call ramg_L2_norm_p(cgR,1,norm_0)
      norm_1 = norm_0
      norm_c = norm_0
      norm_e = norm_0 * ramg_eps
      write(*,*)'norm_0 ', norm_0

      cgZ = 0
      cgZ = cgR

     
!      call ramg_V_cycle(cgZ,cgR,
!     &                  mres_0,mres_n,10,1,
!     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
!      call ramg_G_cycle(cgZ,cgR,
!     &                  mres_0,mres_n,10,1,
!     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
  
!      call ramg_dot(amg_nshg(1),cgZ,cgR,rz)
      call ramg_dot_p(cgZ,cgR,1,rz)
      rz_0 = rz

      cgP = cgZ
      ramg_sol = 0

      iterNum = 1

      do while (norm_c.gt.norm_e)
!      do while (iterNum .le. 24)
        call ramg_PPEAps(cgQ,cgP,diag,colm,rowp,lhsK,lhsP,
     &                ilwork,BC,iBC,iper)
!        call ramg_calcAv(amg_A_colm(1)%p,amg_A_rowp(1)%p,
!     &                    amg_A_lhs(1)%p,amg_nshg(1),amg_nnz(1),
!     &                    cgQ,cgP)
!         call ramg_dot(amg_nshg(1),cgP,cgQ,pq)
         call ramg_dot_p(cgP,cgQ,1,pq)
         alpha = rz/pq
         ramg_sol = ramg_sol + alpha*cgP
         cgR = cgR - alpha*cgQ
!         call ramg_L2_norm(amg_nshg(1),cgR,tmp)
         call ramg_L2_norm_p(cgR,1,tmp)
         norm_c = tmp
         cgZ = 0!cgR
         cgZ = cgR
!         call ramg_V_cycle(cgZ,cgR,
!     &                  mres_0,mres_n,10,1,
!     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
!         call ramg_G_cycle(cgZ,cgR,
!     &                  mres_0,mres_n,10,1,
!     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
         tmp = rz
!         call ramg_dot(amg_nshg(1),cgZ,cgR,rz)
         call ramg_dot_p(cgZ,cgR,1,rz)
         beta = rz/tmp
         cgP = beta*cgP + cgZ
         if (myrank.eq.0) then 
            write(*,*)iterNum,norm_c/norm_0
         end if
         iterNum = iterNum + 1
      enddo

      end subroutine !ramg_CG

!*******************************************
!      ramg Interface
!      Interface for libLES
!******************************************* 
      subroutine ramg_interface(
     &           colm,rowp,lhsK,lhsP,flowDiag,
     &           mcgR,mcgZ,
     &           ilwork,BC,iBC,iper
     &           )
      use ramg_data
      include "common.h"
      include "mpif.h"
      integer,intent(in),dimension(nshg+1)             :: colm
      integer,intent(in),dimension(nnz_tot)            :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)     :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)     :: lhsP
      real(kind=8),intent(in),dimension(nshg)          :: mcgR
      real(kind=8),intent(inout),dimension(nshg)       :: mcgZ
      integer, intent(in), dimension(nlwork)           :: ilwork
      integer, intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)   :: BC
      real(kind=8),dimension(nshg,nflow)               :: dummyR
      real(kind=8),dimension(nshg)                     :: myR
      real(kind=8),intent(in),dimension(nshg,4)        :: flowDiag
      
      real(kind=8)                                     :: pos,tlen
      integer::nTest
      integer::i,j,k
      character*10                                     :: fname

      dummyR = 0

      mcgZ = 0

      if (irun_amg_prec.ne.0) then
!      if (.false.) then
        if (ramg_flag.eq.0) then
!            call ramg_dump_ir(ncorp_map,mcgR,nshg,1,'mcgR      ')
!            call ramg_ggb_setup(ilwork,BC,iBC,iper)
        end if
       myR = mcgR*amg_ppeDiag(1)%p
       call ramg_driver(colm,rowp,lhsK,lhsP,
     &           nshg,nnz_tot,nflow,
     &           dummyR,myR,
     &           nlwork,ilwork,ndofBC,BC,iBC,iper,mcgZ,1,
!     &           ramg_eps,irun_amg_prec+10)
     &           ramg_eps,1+10)
        ramg_flag = ramg_flag + 1
        mcgZ = mcgZ*amg_ppeDiag(1)%p
      else
          mcgZ = mcgR
          !call ramg_zeroOut(mcgZ,ilwork,BC,iBC,iper)
      endif
      
      !call ramg_dump(mcgR,nshg,'mcgR      ')
 
      end subroutine ramg_interface

!*******************************************      
!        ramg Driver
!        1. extract PPE / system
!        2. call coarsening
!        3. solve     
!********************************************

      !*************************************
      !    ramg_driver
      !     Input:  global matrix,# of systems
      !             control params
      !     Output: solution
      !*************************************
      subroutine ramg_driver(
     &           colm,rowp,lhsK,lhsP,
     &           nshg,nnz_tot,nflow,
     &           mrhs,pperhs,
     &           nlwork,ilwork,ndofBC,BC,iBC,iper,
     &           ramg_sol,n_sol,
     &           amg_eps,amg_mode
     &           )
      
      use ramg_data
      implicit none
      
      !***********parameters**************
      !the matrix
      integer,intent(in)                               :: nshg
      integer,intent(in)                        :: nnz_tot
      integer,intent(in)                               :: nflow
      !the matrix
      integer,intent(in),dimension(nshg+1)             :: colm
      integer,intent(in),dimension(nnz_tot)            :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)     :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)     :: lhsP
      !the forcing term, rhs
      real(kind=8),intent(in),dimension(nshg,nflow)    :: mrhs
      real(kind=8),intent(in),dimension(nshg)       :: pperhs
      ! the boundary info
      integer, intent(in)                              :: nlwork
      integer, intent(in), dimension(nlwork)           :: ilwork
      integer, intent(in)                              :: ndofBC
      integer, intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)   :: BC
      !the output solution
      integer, intent(in)                              :: n_sol
      real(kind=8),intent(inout),dimension(nshg,n_sol) :: ramg_sol
      real(kind=8),dimension(:,:),allocatable          :: amg_sol

      !the tolerance
      real(kind=8),intent(in)                          :: amg_eps
      !control run mode
      integer,intent(inout)                            :: amg_mode
      ! my AMG parameter
      !*********parameters end**************

      !*****local variables*****************
      integer                           :: mem_err,mem_err_s
      integer                           :: ppe_nnz
      real(kind=8),dimension(nshg)      :: precrhs
      real(kind=8)                      :: res_i,res_o
      integer                           :: ramg_mode,i
      logical                           :: extmtx
      real(kind=8)                      :: cpusec(10)

      call cpu_time(cpusec(1))

      extmtx = .false.
          if (amg_mode.gt.0) then
            if (amg_nnz(1).eq.0) then
            call cpu_time(cpusec(2))
            call ramg_extract_ppe(colm,rowp,lhsK,lhsP,
     &         mrhs,ilwork,BC,iBC,iper)
            call cpu_time(cpusec(3))
            !ramg_time(2) = ramg_time(2) + cpusec(3)-cpusec(2)
            endif
            allocate(amg_sol(nshg,1))
          else
              extmtx = .true.
              amg_mode = - amg_mode
c              call samgread;
              allocate(amg_sol(amg_nshg(1),1))
          end if
          amg_sol = 0

      if (amg_mode .eq. 11 ) then     ! solve PPE les Precondition ramg
         call cpu_time(cpusec(2))
         call ramg_V_cycle(amg_sol,pperhs,
     &                     res_i,res_o,10,1,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
         !call ramg_zeroOut(amg_sol,ilwork,BC,iBC,iper)
         call ramg_G_cycle(amg_sol,pperhs,
     &                     res_i,res_o,10,1,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
         !call ramg_zeroOut(amg_sol,ilwork,BC,iBC,iper)
!         call ramg_V_cycle(amg_sol,pperhs,
!     &                     res_i,res_o,10,1,
!     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
         call cpu_time(cpusec(3))
         ramg_time(4)=ramg_time(4)+cpusec(3)-cpusec(2)
      else if (amg_mode .eq. 12) then ! solve PPE CG Precondition SAMG
          ramg_mode = 2
          if (ramg_flag.eq.0) then
              ramg_mode = 1
          end if
          call cpu_time(cpusec(2))
          call ramg_V_cycle(amg_sol,pperhs,
     &                      res_i,res_o,ramg_mode,1,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
          call cpu_time(cpusec(3))
!          if (ramg_mode .eq. 1) then
!              ramg_time(3) = ramg_time(3) + cpusec(3) - cpusec(2)
!          else
!              ramg_time(4) = ramg_time(4) + cpusec(3) - cpusec(2)
!          end if
      else if (amg_mode .eq. 2) then  ! solve PPE SAMG
           ramg_mode = 0
           
           call ramg_V_cycle(amg_sol,amg_A_rhs(1)%p,
     &                     res_i,res_o,ramg_mode,1,1,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
      else if (amg_mode .eq. 1) then ! solve PPE stand alone
      else if (amg_mode .eq. 3) then ! solve PPE CG
          call ramg_prep(ilwork,BC,iBC,iper)
              call ramg_ggb_setup(ilwork,BC,iBC,iper)
          call ramg_CG(amg_sol,amg_A_rhs(1)%p,mrhs,4,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
      end if

      if (extmtx.eq..true.) then
          call ramg_deallocate(-1)
          if (amg_mode.ne.2) then
          do i=2,ramg_levelx
             call ramg_deallocate(i)
          enddo
          endif
      else
          ramg_sol = amg_sol
      end if

      deallocate(amg_sol)

      return

      end subroutine ramg_driver 
!*******************************************
!      <EOF> ramg Driver
!*******************************************      
