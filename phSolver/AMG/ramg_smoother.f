!************************************************************
!      ramg_calcPPEv
!      calculate u = PPE*v ( parallelly using lesLib calls)
!************************************************************
      subroutine ramg_calcPPEv(colm,rowp,lhsK,lhsP,
     &                   ilwork,BC,iBC,iper,
     &                         u,v)
      use ramg_data
      include "common.h"
      include "mpif.h"
      
      integer,intent(in),dimension(nshg+1) :: colm
      integer,intent(in),dimension(nnz_tot) :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot) :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot) :: lhsP
      integer,intent(in),dimension(nlwork)  :: ilwork
      integer,intent(in),dimension(nshg)    :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC
      
      real(kind=8),intent(in),dimension(nshg) :: v
      real(kind=8),intent(inout),dimension(nshg) :: u

      real(kind=8),dimension(nshg,3) :: utmp
      real(kind=8),dimension(nshg,4) :: utmp4
      integer                         :: i,j

      call commOut(v,ilwork,1,iper,iBC,BC)
      call fLesSparseApG(colm,rowp,lhsP,v,utmp,nshg,nnz_tot)
      call commIn(utmp,ilwork,3,iper,iBC,BC)
!      call commOut(utmp,ilwork,3,iper,iBC,BC)

      do i=1,3
         utmp(:,i) = utmp(:,i)*ramg_flowDiag%p(:,i)**2
      enddo

      do i=1,3
         utmp4(:,i) = -utmp(:,i)
      enddo

      utmp4(:,4) = v

      call commOut(utmp4,ilwork,4,iper,iBC,BC)
      call fLesSparseApNGtC(colm,rowp,lhsP,utmp4,u,nshg,nnz_tot)
      call commIn(u,ilwork,1,iper,iBC,BC)
!      call commOut(u,ilwork,1,iper,iBC,BC)
!     There is a slight modify at lesSparse.f:
!     row(20*nNodes) ==> row(nnz_tot)
      
      end subroutine ! ramg_calcPPEv

      
!************************************************************
!      ramg_PPE_jacobi
!      calculate global u=(D^-1)*2(-PPE*v + Dr)
!
!      Ap product is global
!      
!************************************************************
      subroutine ramg_PPE_jacobi(colm,rowp,lhsK,lhsP,
     &                   ilwork,BC,iBC,iper,
     &                         r,u,v,resin,resout)
      use ramg_data
      include "common.h"
      include "mpif.h"
      
      integer,intent(in),dimension(nshg+1) :: colm
      integer,intent(in),dimension(nnz_tot) :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot) :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot) :: lhsP
      integer,intent(in),dimension(nlwork)  :: ilwork
      integer,intent(in),dimension(nshg)    :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC
      
      real(kind=8),intent(in),dimension(nshg) :: v,r
      real(kind=8),intent(inout),dimension(nshg) :: u
      real(kind=8),intent(inout)  :: resin,resout

      real(kind=8),dimension(nshg) :: utmp,vtmp
      
      integer                         :: i,j

      call ramg_calcPPEv(colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper,
     &                   utmp,v)! Pv
      vtmp = r-utmp
      call ramg_L2_norm(nshg,vtmp,resin)
      do i = 1,nshg
         u(i) = ramg_flowDiag%p(i,5)*(r(i)-utmp(i))+v(i)
      enddo
      call ramg_calcPPEv(colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper,
     &                   utmp,u)
      vtmp = r-utmp
      call ramg_L2_norm(nshg,vtmp,resout)
      resin = sqrt(resin)
      resout = sqrt(resout)

      end subroutine ! ramg_PPE_jacobi
 

!************************************************************
!      ramg_jacobi
!      calculate u = D^-1 * [ -( L+U ) v + r ]
!
!      and also the residule in and out resout = | r-Au |
!                                       resin = | r-Av | 
!      
!************************************************************
      subroutine ramg_jacobi(Acolm,Arowp,Alhs,Anshg,Annz_tot,!A
     &                      r,u,v,resin,resout)
     
      include "common.h"
      integer,intent(in)                   :: Anshg,Annz_tot
      integer,intent(in),dimension(Anshg+1) :: Acolm
      integer,intent(in),dimension(Annz_tot) :: Arowp
      real(kind=8),intent(in),dimension(Annz_tot) :: Alhs
      real(kind=8),intent(in),dimension(Anshg) :: v,r
      real(kind=8),intent(inout),dimension(Anshg) :: u
      real(kind=8),intent(inout)              :: resin,resout
      real(kind=8)                        :: tmp,tmpin,tmpout
      
      integer                              :: i,j,k,p

      real(kind=8)                     :: damp_jacobi

      u = 0
      resin = 0
      resout = 0
      damp_jacobi = 1.0/ramg_relax
      do i=1,Anshg
         tmpin = 0
         do k=Acolm(i),Acolm(i+1)-1
            tmpin = tmpin + u(Arowp(k))*Alhs(k)
         enddo
         tmpin = tmpin-r(i)
         resin = resin + tmpin*tmpin
      enddo
      do i=1,Anshg
         do k = Acolm(i)+1,Acolm(i+1)-1
            j = Arowp(k)
            tmp = Alhs(k)*v(j)
            u(i) = u(i) - tmp
         enddo
         u(i) = damp_jacobi*u(i) + r(i)
         u(i) = u(i)/Alhs(Acolm(i))
      enddo
      do i=1,Anshg
         tmpout = 0
         do k=Acolm(i),Acolm(i+1)-1
            tmpout = tmpout + u(Arowp(k))*Alhs(k)
         enddo
         tmpout = tmpout-r(i)
         resout = resout + tmpout*tmpout
      enddo
      resin = sqrt(resin)
      resout = sqrt(resout)
       
      end subroutine ! ramg_jacobi

!************************************************************
!      ramg_gauss
!      
!      forward and backward, defined in fwdbck
!      Gauss-Seidel smoothing, u = gaussseidel(M,r)
!
!      and also the residule in and out resout = | r-Au |
!                                       resin = | r-Av | 
!      
!************************************************************
      subroutine ramg_gauss(acolm,arowp,alhs,anshg,annz_tot,!A
     &                      r,u,v,resin,resout,fwdbck,clevel,
     &                      ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in)                   :: anshg,annz_tot
      integer,intent(in),dimension(anshg+1) :: acolm
      integer,intent(in),dimension(annz_tot) :: arowp
      real(kind=8),intent(in),dimension(annz_tot) :: alhs
      real(kind=8),intent(in),dimension(anshg) :: v,r
      real(kind=8),intent(inout),dimension(anshg) :: u
      real(kind=8),intent(inout)              :: resin,resout
      integer,intent(in)                      :: fwdbck!1=fwd,2=bck
      integer,intent(in)                      :: clevel
      integer,intent(in),dimension(nlwork)        :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC

     
      real(kind=8)                        :: tmp,tmpin,tmpout,tmp2
      real(kind=8),dimension(nshg)        :: gs_tmp
      integer,dimension(nshg)             :: gs_flag
      integer                             :: istr,iend,iint
      integer                              :: i,j,k,p,ki,kj,kp
!      integer            :: itag,iacc,iother,numseg,isgbeg,itkbeg
!      integer            :: numtask

      u = 0!v
      resin = 0
      resout = 0
      if (fwdbck.eq.1) then
          istr = 1
          iend = anshg
          iint = 1
      else
          istr = anshg
          iend = 1
          iint = -1
      end if

      !if (anshg.eq.nshg) then  ! If we do level-1 smoothing
          gs_tmp = 0
      !end if
      
!      do i=1,anshg
!         tmpin = 0
!         do k=acolm(i),acolm(i+1)-1
!            tmpin = tmpin + u(arowp(k))*alhs(k)
!         enddo
!         tmpin = tmpin-r(i)
!         resin = resin + tmpin*tmpin
!      enddo

      !write(*,*)'mcheck:ggb:start_gauss',resin

      do i=istr,iend,iint
         tmp = 0
         if (clevel.ne.0) then
             ki = CF_map(clevel)%p(i)
         else
             ki = i
         end if
         !ki = i
         do k = acolm(ki)+1,acolm(ki+1)-1
            j = arowp(k)
            if (clevel.ne.0) then
                kj = CF_revmap(clevel)%p(j)
            else
                kj = j
            end if
            tmp = tmp + alhs(k)*u(j)
         enddo
         gs_tmp(ki) = tmp
         u(ki) = r(ki) - tmp
         u(ki) = u(ki)/alhs(acolm(ki))
      enddo

      !write(*,*)'mcheck:ggb:end gauss 1'
      
      do i=istr,iend,iint
         u(i) = r(i) - gs_tmp(i)
         u(i) = u(i)/alhs(acolm(i))
      enddo
!      do i=1,anshg
!         tmpout = 0
!         do k=acolm(i),acolm(i+1)-1
!            tmpout = tmpout + u(arowp(k))*alhs(k)
!         enddo
!         tmpout = tmpout-r(i)
!         resout = resout + tmpout*tmpout
!      enddo
      resin = 0!sqrt(resin)
      resout = 0!sqrt(resout)
      
      end subroutine ! ramg_gauss

!**********************************************
!      ramg dense direct solver routines.
!      
!      ramg_nr_ludcmp
!      LU Decomposition from Numerical Recipies
!      
!      ramg_nr_lubksb
!      Back substitution after LU
!      
!      ramg_sparse2dense
!      Sparse to dense form
!      
!      ramg_direct_LU
!      Outside package for direct solve sparse
!        matrix by LU
!**********************************************
      subroutine ramg_nr_ludcmp(a,n,np,indx,d)
      integer          ::  n,np,dnr_NMAX
      integer,dimension(n)   :: indx
      real(kind=8)           :: d
      real(kind=8),dimension(np,np)    :: a
      parameter (dnr_NMAX=500)
      integer i,imax,j,k
      real(kind=8)           ::  aamax,dum,summ
      real(kind=8),dimension(dnr_NMAX) :: vv
      d = 1.
      do i=1,n
         aamax = 0.
         do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax = abs(a(i,j))
         enddo 
         if (aamax.eq.0) pause 'singular matrix in LU decomp'
         vv(i) = 1./aamax
      enddo 
      do j=1,n
         do i=1,j-1
            summ = a(i,j)
            do k=1,i-1
               summ = summ - a(i,k)*a(k,j)
            enddo
            a(i,j) = summ
         enddo
         aamax = 0.
         do i=j,n
            summ = a(i,j)
            do k = 1,j-1
               summ = summ - a(i,k)*a(k,j)
            enddo
            a(i,j) = summ
            dum = vv(i)*abs(summ)
            if (dum.ge.aamax) then
                imax = i
                aamax = dum
            end if
         enddo
         if (j.ne.imax) then
             do k=1,n
                dum = a(imax,k)
                a(imax,k) = a(j,k)
                a(j,k) = dum
             enddo 
             d = -d
             vv(imax) = vv(j)
         end if
         indx(j) = imax
         if (a(j,j).eq.0) then
             a(j,j) = 1.0e-20
             write(*,*) 'tiny tiny singular'
         end if
         if (j.ne.n) then
             dum = 1./a(j,j)
             do i=j+1,n
                a(i,j) = a(i,j)*dum
             enddo
         end if
      enddo
      end subroutine ! ramg_nr_ludcmp

      subroutine ramg_nr_lubksb(a,n,np,indx,b)
      integer                        :: n,np
      integer,dimension(n)           :: indx
      real(kind=8),dimension(np,np)  :: a
      real(kind=8),dimension(n)      :: b
      integer :: i,ii,j,ll
      real(kind=8) :: summ
      ii = 0
      do i=1,n
         ll = indx(i)
         summ = b(ll)
         b(ll) = b(i)
         if (ii.ne.0) then
             do j=ii,i-1
                summ = summ-a(i,j)*b(j)
             enddo
         else if (summ.ne.0) then
             ii = 1
         end if
         b(i) = summ
      enddo
      do i=n,1,-1
         summ = b(i)
         do j=i+1,n
            summ = summ - a(i,j)*b(j)
         enddo
         b(i)=summ/a(i,i)
      enddo
      end subroutine ! ramg_nr_lubksb

      subroutine ramg_direct_LU(Acolm,Arowp,Alhs,Anshg,Annz_tot,
     &                          Arhs,Asol)

      integer,intent(in)                   :: Anshg,Annz_tot
      integer,intent(in),dimension(Anshg+1) :: Acolm
      integer,intent(in),dimension(Annz_tot) :: Arowp
      real(kind=8),intent(in),dimension(Annz_tot) :: Alhs
      real(kind=8),intent(in),dimension(Anshg) :: Arhs
      real(kind=8),intent(inout),dimension(Anshg) :: Asol

      real(kind=8),dimension(Anshg,Anshg)      :: mtxA
      integer,dimension(Anshg)             :: indx
      real(kind=8)                     :: d

      integer :: i,j

      call ramg_sparse2dense(Acolm,Arowp,Alhs,Anshg,Annz_tot,mtxA)

      Asol = Arhs
      
      call ramg_nr_ludcmp(mtxA,Anshg,Anshg,indx,d)
      call ramg_nr_lubksb(mtxA,Anshg,Anshg,indx,Asol)
      
      end subroutine ! ramg_direct_LU

      subroutine ramg_sparse2dense(Acolm,Arowp,Alhs,Anshg,Annz_tot,
     &                             mtxA)

      integer,intent(in)                   :: Anshg,Annz_tot
      integer,intent(in),dimension(Anshg+1) :: Acolm
      integer,intent(in),dimension(Annz_tot) :: Arowp
      real(kind=8),intent(in),dimension(Annz_tot) :: Alhs
      real(kind=8),intent(inout),dimension(Anshg,Anshg) :: mtxA

      integer                       :: i,j,k,ip,jp,kp

      mtxA = 0
      do i=1,Anshg
         do j= Acolm(i),Acolm(i+1)-1
            k = Arowp(j)
            mtxA(i,k) = Alhs(j)
         enddo
      enddo

      end subroutine !ramg_sparse2dense
      
 
!**********************************************
!      ramg_direct_solve
!**********************************************
      subroutine ramg_direct_solve(Acolm,Arowp,Alhs,Anshg,Annz_tot,
     &             Arhs,Asol,ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      integer,intent(in)                   :: Anshg,Annz_tot
      integer,intent(in),dimension(Anshg+1) :: Acolm
      integer,intent(in),dimension(Annz_tot) :: Arowp
      real(kind=8),intent(in),dimension(Annz_tot) :: Alhs
      real(kind=8),intent(in),dimension(Anshg) :: Arhs
      real(kind=8),intent(inout),dimension(Anshg) :: Asol
      integer,intent(in),dimension(nlwork)        :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC

      real(kind=8)              :: mres_i,mres_n,mres_0
      real(kind=8),dimension(Anshg) :: myV
      integer                      :: titer

!      Asol=Arhs
!      return

      !write(*,*)'mcheck direct solve begin,',myrank
      if (iamg_c_solver.eq.0) then
          Asol=Arhs
      else if ( iamg_c_solver .eq. 1) then ! direct solve using G-S
      call ramg_gauss(Acolm,Arowp,Alhs,Anshg,Annz_tot,Arhs,
     &                 myV,Asol,mres_0,mres_n,1,0,
     &                 ilwork,BC,iBC,iper)
      titer = 1
      !write(*,*)myrank,'iter:',iter,'res:',mres_0,mres_n
      do while (mres_n > mres_0*(1e-3) )
         Asol = myV
         call ramg_gauss(Acolm,Arowp,Alhs,Anshg,Annz_tot,Arhs,
     &                    myV,Asol,mres_i,mres_n,1,0,
     &                    ilwork,BC,iBC,iper)
         titer = titer+1
         if (titer.gt.99) then
             !write(*,*)'C_Solve exceed max iter',mres_n/mres_0
             exit
         end if
         !write(*,*)myrank,'iter:',iter,'res:',mres_i,mres_n
      enddo
      Asol = myV
      else ! direct solve using dense-matrix LU-decomposition
          call ramg_direct_LU(Acolm,Arowp,Alhs,Anshg,Annz_tot,
     &                        Arhs,Asol)
      endif
      !call ramg_commOut(Asol,ramg_levelx,1,iper,iBC,BC)
      !write(*,*)'mcheck direct solve end,',myrank

      end subroutine ! ramg_direct_solve


      subroutine ramg_smoother(level,xnew,xold,xrhs,resin,resout,
     &                         colm,rowp,lhsK,lhsP,
     &                         ilwork,BC,iBC,iper,fwdbck)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      
      integer level
      real(kind=8),dimension(amg_nshg(level)),intent(in) ::xrhs
      real(kind=8),dimension(amg_nshg(level)),intent(inout) 
     &                            :: xnew,xold
      real(kind=8),intent(inout) :: resin,resout
      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp

      real(kind=8),intent(in),dimension(9,nnz_tot)    :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)    :: lhsP

      integer,intent(in),dimension(nlwork)   :: ilwork
      integer,intent(in),dimension(nshg)     :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC
      integer fwdbck

      integer i,k
      logical :: jacobi

      jacobi = (ramg_relax.gt.0)

      k = 1
      do i=1,k
          if ( (mlsDeg.gt.0) .and. (level.ge.1) ) then
              call ramg_mls_apply(xnew,xold,xrhs,level,
     &                  colm,rowp,lhsK,lhsP,
     &                  ilwork,BC,iBC,iper,
     &                  ramg_res_i,ramg_r_i,fwdbck)
          else if (jacobi) then
             call ramg_jacobi(amg_A_colm(level)%p,
     &                 amg_A_rowp(level)%p,amg_A_lhs(level)%p,
     &                 amg_nshg(level),amg_nnz(level),
     &                 xrhs,xnew,xold,resin,resout)
          else
          call ramg_gauss(amg_A_colm(level)%p,amg_A_rowp(level)%p,
     &                    amg_A_lhs(level)%p,amg_nshg(level),
     &                    amg_nnz(level),xrhs,xnew,xold,
     &                    resin,resout,fwdbck,level,
     &                    ilwork,BC,iBC,iper)
          endif
          xold = xnew
      enddo

      end subroutine !ramg_smoother
