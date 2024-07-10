      module spebc
      
      integer  nfin,  nelint, npin, npint, nfint
      real*8   sang, deltaint
      real*8   xnrml,  ynrml,  znrml,  aR, aI
      
      

      real*8, allocatable :: xyn(:), xynin(:)
      real*8, allocatable :: xcyl(:,:),xintl(:,:,:),xsinfin(:,:,:)
      real*8, allocatable :: xtsX(:), xtsY(:), xtsC(:)

      integer, allocatable :: ien2D(:,:), nen1(:), elcnfin(:,:)
      integer, allocatable :: nrint(:), imax(:) 
  
      end module

c-----------------------------------------------------------------------
c allocate the spebc arrays
c-----------------------------------------------------------------------
      subroutine setSPEBC(numnp,nsd)
      
      use spebc

      allocate (xyn(numnp))
      allocate (xynin(numnp))
      allocate (xcyl(numnp,nsd))
      allocate (nen1(numnp))

c      allocate (elcnpin(numnp))
c      allocate (xsi(numnp,nsd))
      

      
      return
      end

