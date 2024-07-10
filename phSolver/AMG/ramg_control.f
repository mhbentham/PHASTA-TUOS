!*********************************************************
!      ramg control 
!      Control AMG preparation/solve process
!      will have dynamic control (based on solver vecters)
!*********************************************************

      subroutine ramg_control
      use ramg_data
      include "common.h"

      integer i

      if (iamg_init .eq. 0 ) then
          ! do init
          ramg_ppe_flag = 0
          ramg_setup_flag = 0
          ramg_iai_flag = 0
          iamg_init = 1
          ramg_flag = 0
          ramg_time = 0
          ramg_ggb_flag = 0
      else
       ramg_ppe_flag   = mod(ramg_ppe_flag+1   ,iamg_ppe_frez)
       ramg_setup_flag = mod(ramg_setup_flag+1 ,iamg_setup_frez)
       ramg_iai_flag   = mod(ramg_iai_flag+1   ,iamg_iai_frez)
       ramg_ggb_flag   = mod(ramg_ggb_flag+1   ,iamg_ppe_frez)
      end if 
      
      end subroutine !ramg_control

!*********************************************************
!     Ramg preparation
!*********************************************************

      subroutine ramg_prep(ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in),dimension(nlwork)      :: ilwork
      integer,intent(in),dimension(nshg)        :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC

      logical                  :: maxstopsign,mxs2
      integer                  :: i,p,p2

      if (ramg_setup_flag.ne.0) return
      maxstopsign = .false.
      i = 1
      do while ((i+1.le.iamg_nlevel) .and. (.not.maxstopsign))
         if (amg_nshg(i+1).ne.0) then
             call ramg_deallocate(i+1)
         end if
         call ramg_coarse_setup(i,i+1,strong_eps,iamg_interp,
     &        ramg_trunc,
     &        ilwork,BC,iBC,iper,nshg,nlwork,ndofBC)
         call ramg_calcITAI(i,i+1,maxstopsign)
         i = i+1
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         maxstopsign = (amg_nshg(i).eq.amg_nshg(i-1)) 
         IF (.true.) THEN
         !call ramg_checkcoarse(i,ilwork,BC,iBC,iper,maxstopsign)
         p = 1
         if (maxstopsign) p = 0
         call MPI_AllReduce(p,p2,1,MPI_INTEGER,MPI_SUM,
     &                MPI_COMM_WORLD,ierr)
         if (p2.eq.0) then
           maxstopsign = .true.
           !write(*,*)'mcheck stopped'
         else
           maxstopsign = .false.
         endif
         ENDIF
      enddo
      ramg_levelx = i-1
      call ramg_check_mem

      if ( (irun_amg_prec.eq.1).and.(mlsDeg.gt.0) ) then
      deallocate(amg_A_colm(1)%p)
      deallocate(amg_A_rowp(1)%p)
      deallocate(amg_A_lhs(1)%p)
      endif

      !call ramg_ggb_setup(ilwork,BC,iBC,iper)

      end subroutine !ramg_prep

      subroutine ramg_checkcoarse(level1,ilwork,BC,iBC,iper,
     &           cfstop)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in) :: level1
      integer,intent(in),dimension(nlwork)      :: ilwork
      integer,intent(in),dimension(nshg)        :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC

      logical,intent(inout)                  :: cfstop
      integer                  :: i,numtask,itkbeg,numseg,iacc
      integer                  :: j,k,p
      integer,allocatable,dimension(:) 
     &            :: subcfstop,subnnz,subcf,subcfrev,subnei
      real(kind=8) :: rhoratio
      logical      :: subneireduce

      IF (( iamg_reduce.le.1).or.(numpe.gt.1)) THEN
!      IF (.TRUE.) THEN
      if (numpe.ge.2) then 
          numtask = ilwork(1)+1
      else 
          numtask = 1
      endif
      
      allocate(subcfrev(numpe))
      subcfrev = 0
      allocate(subnei(numtask))
      subnei = 0
 
      subnei(1) = myrank+1
      subcfrev(subnei(1))=1
      itkbeg = 1
      do i=2,numtask
         subnei(i)=ilwork(itkbeg+3)+1
         subcfrev(subnei(i)) = i
         itkbeg = itkbeg + 4 + 2*ilwork(itkbeg+4)
      enddo
   
      ELSE !reduced
          numtask = rmapmax-1
          allocate(subcfrev(numtask))
          allocate(subnei(numtask))
          do i=1,numtask
             subcfrev(i) = i
          enddo
      ENDIF

      allocate(subcfstop(numtask))
      allocate(subcf(numtask))
      allocate(subnnz(numtask))
      subcfstop = 0
      subcf = 0
      subnnz = 0


      do i = 1,amg_nnz(level1)
         k = amg_A_rowp(level1)%p(i)
         p = iabs(amg_paramap(level1)%p(k))
!         if (iamg_reduce.gt.1) then
!            p = 1 ! for reduced only
!         endif
         p = subcfrev(p)
         subnnz(p) = subnnz(p) + 1
      enddo
      do i = 1,amg_nshg(level1)
         p = iabs(amg_paramap(level1)%p(i))
!         if (iamg_reduce.gt.1) then
!            p = 1 ! for reduced only
!         endif
         p = subcfrev(p)
         subcf(p) = subcf(p) + 1
      enddo

      do i=1,numtask
         IF ((iamg_reduce.le.1).or.(numpe.gt.1)) THEN
!         IF (.TRUE.) THEN
         p = subnei(i)
         subneireduce = (p.ne.(myrank+1))
         ELSE !reduced
             subneireduce = (i.gt.iamg_reduce)
         END IF
         if ((subcf(i).lt.30).and.subneireduce) then
             subcfstop(i) = 1
         endif
         if (.not.subneireduce) then
             if (subcf(i).lt.200) then
                 subcfstop(i) = 1
             endif
             if ((subnnz(i)/(subcf(i)**2)).gt.0.6) then
                 subcfstop(i) = 1
             endif
         end if
      enddo
  
      cfstop = .true.
      do i=1,numtask
         if (subcfstop(i).eq.0) then
             cfstop = .false.
             exit
         endif
      enddo
    

      IF ((iamg_reduce.le.1).or.(numpe.gt.1)) THEN
      ! if interior is coarsened, boundary should be fixed whatever
      if (subcfstop(1).eq.1) then
          cfstop = .true.
      endif
      ENDIF

      if (.not.cfstop) then
      do i=1,amg_nshg(level1)
         p = iabs(amg_paramap(level1)%p(i))
!         if (iamg_reduce.gt.1) then
!             p = 1 ! reduced only
!         endif
         p = subcfrev(p)
         if (subcfstop(p).eq.1) then
         amg_paramap(level1)%p(i) = -iabs(amg_paramap(level1)%p(i))
         endif
      enddo
      endif

      deallocate(subcfrev)
      deallocate(subnei)
      deallocate(subcfstop)
      deallocate(subcf)
      deallocate(subnnz)
     
      end subroutine ! ramg_checkcoarse

!*********************************************************
!      <EOF> RAMG Control
!*********************************************************
