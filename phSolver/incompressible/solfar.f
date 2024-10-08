      subroutine SolFlow(y,          ac,         u,
     &                   banma,
     &                   yold,       acold,      uold,
     &                   x,          iBC,
     &                   BC,         res,             
     &                   nPermDims,  nTmpDims,  aperm,
     &                   atemp,      iper,       
     &                   ilwork,     shp,        shgl, 
     &                   shpb,       shglb,      rowp,     
     &                   colm,       lhsK,       lhsP, 
     &                   solinc,     rerr,       tcorecp,
     &                   GradV,      elemvol_global,
     &                   avgxcoordf, avgycoordf, avgzcoordf,
     &                   svLS_lhs, svLS_ls, svLS_nFaces)
c
c!----------------------------------------------------------------------
c
c! This is the 2nd interface routine to the Farzin's linear equation 
c! solver library that uses the CGP and GMRES methods.
c
c! input:
c!  y      (nshg,ndof)           : Y-variables at n+alpha_f
c!  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c!  yold   (nshg,ndof)           : Y-variables at beginning of step
c!  acold   (nshg,ndof)          : Primvar. accel. at beginning of step
c!  x      (numnp,nsd)            : node coordinates
c!  iBC    (nshg)                : BC codes
c!  BC     (nshg,ndofBC)         : BC constraint parameters
c!  iper   (nshg)                : periodic nodal information
c
c! output:
c!  res    (nshg,nflow)           : preconditioned residual
c!  y      (nshg,ndof)           : Y-variables at n+alpha_f
c!  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c
c
c! The followings are preliminary steps required to use Farzin''s
c! solver library.  New way of writing has to be used such as
c
c!          |  K     G | | du |    | Rmom  |
c!          |          | |    | =  |       |
c!          | G^t    C | | dp |    | Rcon  |
c
c!          |     E    | | dT | =  | Rtemp |
c
c!     where
c
c!      xKebe : K_ab = dRmom_a/du_b    xTe : E_ab = dRtemp_a/dT_b 
c
c!              G_ab = dRmom_a/dp_b
c!      xGoC  :
c!              C_ab = dRcon_a/dp_b       
c
c!              resf = Rmon Rcon       rest = Rtemp
c
c  
c! Zdenek Johan,  Winter 1991.  (Fortran 90)
c! Juin Kim, Summer 1998. (Incompressible flow solver interface)
c! Alberto Figueroa.  CMM-FSI
c----------------------------------------------------------------------
c
      use pointer_data
c!#ifdef AMG      
c!      use ramg_data
c!#endif     
        
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      include "svLS.h"  !MB, include svLS

      TYPE(svLS_lhsType) svLS_lhs   ! MB
      TYPE(svLS_lsType) svLS_ls     ! MB
c     
      real*8    y(nshg,ndof),             ac(nshg,ndof),
     &          yold(nshg,ndof),          acold(nshg,ndof),
     &          u(nshg,nsd),              uold(nshg,nsd),
     &          x(numnp,nsd),             BC(nshg,ndofBC),
     &          res(nshg,nflow),          tmpres(nshg,nflow),
     &          flowDiag(nshg,4),         
     &          aperm(nshg,nPermDims),    atemp(nshg,nTmpDims),
     &          sclrDiag(nshg,1),         
     &          lhsK(9,nnz_tot),          lhsP(4,nnz_tot),
!MR CHANGE
     &          GradV(nshg,nsdsq)
!MR CHANGE
      dimension banma(nshg,1)
c
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &          shpb(MAXTOP,maxsh,MAXQPT),
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
      integer   usr(100),                 eqnType,temp,
     &          rowp(nshg*nnz),           colm(nshg+1),
     &          iBC(nshg),                ilwork(nlwork),
     &          iper(nshg)
c
      real*8    yAlpha(nshg,ndof),        acAlpha(nshg,ndof),
     &          uAlpha(nshg,nsd),
     &          lesP(nshg,4),             lesQ(nshg,4),
     &          solinc(nshg,ndof),        CGsol(nshg)

!MR CHANGE
      real*8     tcorecp(2)
!MR CHANGE END
      
      real*8    rerr(nshg,numerr),            rtmp(nshg,4),rtemp
      
      real*8    msum(4),mval(4),cpusec(10)

! BEGIN MB ----------------------------------------------------------------
      integer dof, svLS_nfaces, i, j, k, lesId
      integer, allocatable :: incL(:)
      real*8, allocatable :: faceRes(:), Res4(:,:), Val4(:,:)
! END MB ------------------------------------------------------------------
c!.... Matt Talley's Bubble Coal Control
      real*8 avgxcoordf(coalest), avgycoordf(coalest), avgzcoordf(coalest)

c     
c.... *******************>> Element Data Formation <<******************
c
c
c.... set the parameters for flux and surface tension calculations
c
c
        temp = npro


        idflx = 0 
        if(idiff >= 1 )  idflx= (nflow-1) * nsd  
        if (isurf == 1)  idflx=nflow*nsd
        if (idiff >= 1)  idflx = idflx + 10     ! used in the dissipation timeseries computation  
c.... compute solution at n+alpha
c
c           if (myrank .eq. master) write(*,*) 'idflx = ', idflx

      call itrYAlpha( uold,    yold,    acold,
     &                u,       y,       ac,  
     &                uAlpha,  yAlpha,  acAlpha)
c           if (myrank .eq. master) write(*,*) 'itrYAlpha is done'

c
c.... form the LHS matrices, the residual vector (at alpha)
c
!MR CHANGE
      impistat=1
      impistat2=1
      telmcp1 = TMRC()
!MR CHANGE END
      !if (myrank.eq.master) write (*,*) 'calling ElmGMR'
      call ElmGMR (uAlpha,    yAlpha,     acAlpha,
     &             x,         banma,
     &             shp,       shgl,       iBC,       
     &             BC,        shpb,       shglb,
     &             res,       iper,       ilwork,   
     &             rowp,      colm,       lhsK,      
     &             lhsP,      rerr,       GradV,
     &             elemvol_global,
     &             avgxcoordf,  avgycoordf,  avgzcoordf)
      !if (myrank.eq.master) write (*,*) 'called ElmGMR'

      telmcp2 = TMRC()
      impistat=0
      impistat2=0
!MR CHANGE END

c           if (myrank .eq. master) write(*,*) 'ElmGMR is done'

            tmpres(:,:) = res(:,:)
            iblk = 1
! BEGIN MB ---------------------------------------------------------
!###################################################################
!  Calling svLs to solve
! The following block of code has been implemented just as the svLS
! is implemented in the developBoiling_IB verison for the code
! this will likely lead to mistakes if unedited with syncIO format of
! the current code -- MB, 01 Aug 2024

      IF (svLSFlag.eq.1) THEN

      incL = 1
      dof = 4
      if (.not.allocated(Res4)) then
            allocate(Res4(dof,nshg), Val4(dof*dof, nnz_tot))
      end if

      do i=1, nshg
            Res4(1:dof,i) = res(i,1:dof)
      end do

      do i=1, nnz_tot
            Val4(1:3,i) = lhsK(1:3,i)
            Val4(5:7,i) = lhsK(4:6,i)
            Val4(9:11,i) = lhsK(7:9,i)
            Val4(13:15,i) = lhsP(1:3,i)
            Val4(16,i) = lhsP(4,i)
      end do

      do i=1, nshg
            do j=colm(i), colm(i+1) - 1
                  k =rowp(j)
                  do l=colm(k), colm(k+1) - 1
                        if (rowp(l).eq.i) then
                              Val4(4:12:4,l) = -lhsP(1:3,j)
                              exit
                        end if
                  end do
            end do
      end do

      do j=1, 0
            i = iper(j)

            do k1=colm(i), colm(i+1) - 1
                  if (rowp(k1).eq.i) then
                        i1=k1
                        exit
                  end if
            end do
            
            do k1=colm(j), colm(j+1) - 1
                  if (rowp(k1).eq.j) then
                        j1=k1
                        exit
                  end if 
            end do

            if (i.ne.j) then        ! periodic condition on node

                  do k1=colm(j), colm(j+1) - 1
                        write(*,*)'k1 = ', k1, ': rowp(k)=', rowp(k1)
                        if (rowp(k1).eq.j) then
                              Val4(1:16, k1) = zero
                              Val4(1:16:5, k1) = one
                        else
                              Val4(1:16, k1) = zero
                        end if
                        if (rowp(k1).eq.i) then
                              Val4(1:16:5, k1) = -one
                              write(*,*) ' component (j,i) is modified'
                        end if
                  end do

            end if ! end periodic node condition
          
          
 101   FORMAT(1x, 'VAL4(1,6,11,16) @', 2I7, 4F12.7)
 102   FORMAT(1x, 'Res4(1:4) @', I7, 4F12.7)
 103   FORMAT(1x, 'val4 row # ', I7, 104F12.7)
      
      end do            ! nshg loop
      !if (myrank.eq.master) write(*,*) 'calling svLS_SOLVE'
      !if (myrank.eq.master) write(*,*) 'dof = ', dof
      !if (myrank.eq.master) write(*,*) 'Res4 = ', Res4
      !if (myrank.eq.master) write(*,*) 'Val4 = ', Val4
      call svLS_SOLVE(svLS_lhs, svLS_ls, dof, Res4, Val4)
      !if(myrank.eq.master)write(*,*)'svLS_SOLVE done'

      DO i=1, nshg
            solinc(i,1:dof) = Res4(1:dof,i)
      END DO

      ELSE  ! MB, use lesSolve instead
           

! ###############################################################################
! END MB ------------------------------------------------------------------------

c.... lesSolve : main matrix solver
c
      lesId   = numeqns(1)
      eqnType = 1

!MR CHANGE
      impistat=1
      impistat2=1
      tlescp1 = TMRC()
!MR CHANGE END
c
c.... setup the linear algebra solver
c
      rtmp = res(:,1:4)
      call usrNew ( usr,        eqnType,          aperm,
     &              atemp,      rtmp,             solinc,          
     &              flowDiag,   sclrDiag,         lesP,   
     &              lesQ,       iBC,              BC,
     &              iper,       ilwork,           numpe,
     &              nshg,       nshl,             nPermDims,  
     &              nTmpDims,   rowp,             colm,     
     &              lhsK,       lhsP,             rdtmp,      
     &              nnz_tot,    CGsol )
c
c.... solve linear system
c

c#ifdef AMG     
c      ramg_verb=iamg_verb
c      ! Initial Time Profiling
c      call cpu_time(cpusec(1))
c      !call MPI_Barrier(MPI_COMM_WORLD,ierr)
c      if (irun_amg_prec.ne.0) then
c          call ramg_control
c          !call MPI_Barrier(MPI_COMM_WORLD,ierr)
c          call cpu_time(cpusec(4))
c          if (irun_amg_sa.eq.0) then
c         call ramg_extract_ppe(colm,rowp,lhsK,lhsP,
c     &               tmpres,ilwork,BC,iBC,iper)
c          endif
c          call cpu_time(cpusec(5))
c          if (irun_amg_sa.eq.0) then
c          call ramg_prep(ilwork,BC,iBC,iper)
c          call ramg_init_ilwork(ilwork,BC,iBC,iper)
c!          call ramg_check_paracomm(ilwork,BC,iBC,iper)
c          !deallocate(levelbaseflag)
c          !deallocate(reducecolm)
c          !deallocate(reducerowp)
c          !deallocate(reducelhs)
c          endif
c          if (mlsDeg.gt.0) then
c          call ramg_mls_setup(colm,rowp,lhsK,lhsP,
c     &                        ilwork,BC,iBC,iper)
c          endif
c          if (maxnev.gt.0) then
c        call ramg_ggb_setup(colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
c          endif
c      end if
c
c      call cpu_time(cpusec(6))
c      call cpu_time(cpusec(2))
c#endif
      call myfLesSolve ( lesId, usr )
c#ifdef AMG      
c      call cpu_time(cpusec(3))
c
c      ramg_time(1) = ramg_time(1) + cpusec(3)-cpusec(1)
c      ramg_time(7) = ramg_time(7) + cpusec(6)-cpusec(1)
c
c      ! ramg_time: 1 : local total
c      !            2 : max total
c      !            3 : average total
c      !            4 : local VG-cycle
c      !            5 : max VG-cycle
c      !            6 : average VG-cycle
c      !            7 : local setup time
c      !            8 : max setup time
c      !            9 : average setup time
c      !           10 :
c      !           11 : Ap-product level 1
c      !           12 : Ap-product level >1
c      !           13 : Prolongation/Restriction
c      
c      if (numpe.gt.1) then
c          call MPI_AllReduce(ramg_time(1),ramg_time(2),1,
c     &    MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
c          call MPI_AllReduce(ramg_time(1),ramg_time(3),1,
c     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
c          ramg_time(3)=ramg_time(3)/numpe;
c          call MPI_AllReduce(ramg_time(4),ramg_time(5),1,
c     &    MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
c          call MPI_AllReduce(ramg_time(4),ramg_time(6),1,
c     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
c          ramg_time(6)=ramg_time(6)/numpe;
c          call MPI_AllReduce(ramg_time(7),ramg_time(8),1,
c     &    MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
c          call MPI_AllReduce(ramg_time(7),ramg_time(9),1,
c     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
c          ramg_time(9)=ramg_time(9)/numpe;
c      else
c          ramg_time(2) = ramg_time(1)
c          ramg_time(3) = ramg_time(1)
c          ramg_time(5) = ramg_time(4)
c          ramg_time(6) = ramg_time(4)
c          ramg_time(8) = ramg_time(7)
c          ramg_time(9) = ramg_time(7)
c      endif
c      
c      ! Time profiling output
c      if ((iamg_verb.gt.1).and.(myrank.eq.0)) then
c      write(*,*)
c      write(*,7101)' == ACUSIM DONE, Total Time :',ramg_time(2),
c     &                                             ramg_time(3)
c      write(*,7101)' == AMG == Setup :',ramg_time(8),ramg_time(9)
c      write(*,7101)' == AMG == VG-cycle :',ramg_time(5),ramg_time(6)
c      write(*,7101)' == AMG == Ap(=1)(>1) :',
c     &   ramg_time(11),ramg_time(12)
c      write(*,7101)' == AMG == Ap(Int/Prol) :',
c     &   ramg_time(13),ramg_time(13)
c7101  format(T1,A,T40,f9.2,' sec,  ',f9.2,' sec')
c      write(*,*)
c      endif
c
c#endif     
      
      ! End Time profiling output
      
      call getSol ( usr, solinc )

      if (numpe > 1) then
         call commu ( solinc, ilwork, nflow, 'out')
      endif

      END IF      ! MB, svLS / lesLIB choice condition 
!MR CHANGE 
      tlescp2 = TMRC()
      impistat=0
      impistat2=0

      tcorecp(1) = tcorecp(1) + telmcp2-telmcp1 ! elem. formation
      tcorecp(2) = tcorecp(2) + tlescp2-tlescp1 ! linear alg. solution
!MR CHANGE END
      call rstatic (res, y, solinc) ! output flow stats
c     
c.... end
c     
      return
      end

      subroutine SolSclr(y,          ac,         u,
     &                   yold,       acold,      uold,
     &                   x,          iBC,
     &                   BC,         nPermDimsS,  nTmpDimsS,  
     &                   apermS,     atempS,     iper,       
     &                   ilwork,     shp,        shgl, 
     &                   shpb,       shglb,      rowp,     
     &                   colm,       lhsS,       solinc,
     &                   cfl,
     &                   svLS_lhs_sc, svLS_sc, svLS_nFaces)
c
c----------------------------------------------------------------------
c
c This is the 2nd interface routine to the Farzin''s linear equation 
c solver library.
c
c input:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c  yold   (nshg,ndof)           : Y-variables at beginning of step
c  acold  (nshg,ndof)           : Primvar. accel. variable at begng step
c  x      (numnp,nsd)            : node coordinates
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  iper   (nshg)                : periodic nodal information
c
c output:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c
c
c The followings are preliminary steps required to use Farzin''s
c solver library.  New way of writing has to be used such as
c
c          |     E    | | dS | =  | RScal |
c
c----------------------------------------------------------------------
c
      use pointer_data
        
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      INCLUDE "svLS.h"

      TYPE(svLS_lhsType) svLS_lhs_sc
      TYPE(svLS_lsType) svLS_sc
c     
      real*8    y(nshg,ndof),             ac(nshg,ndof),
     &          yold(nshg,ndof),          acold(nshg,ndof),
     &          u(nshg,nsd),              uold(nshg,nsd),
     &          x(numnp,nsd),             BC(nshg,ndofBC),
     &          res(nshg,1),
     &          flowDiag(nshg,4),
     &          sclrDiag(nshg,1),         lhsS(nnz_tot),
     &          apermS(nshg,nPermDimsS),  atempS(nshg,nTmpDimsS)

c
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &          shpb(MAXTOP,maxsh,MAXQPT),
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
      integer   usr(100),                 eqnType,
     &          rowp(nshg*nnz),           colm(nshg+1),
     &          iBC(nshg),                ilwork(nlwork),
     &          iper(nshg)
c
      real*8    yAlpha(nshg,ndof),        acAlpha(nshg,ndof),
     &          uAlpha(nshg,nsd),
     &          lesP(nshg,1),             lesQ(nshg,1),
     &          solinc(nshg,1),           cfl(nshg)

      !REAL*8 sumtime, timekeeper
      INTEGER dof, svLS_nFaces, i, j, k, l
      INTEGER, ALLOCATABLE :: incL(:)
      REAL*8, ALLOCATABLE :: faceRes(:), Res4(:,:), Val4(:,:)
      
c     
c.... *******************>> Element Data Formation <<******************
c
c.... compute solution at n+alpha
c
      call itrYAlpha( uold,    yold,    acold, 
     &                u,       y,       ac,  
     &                uAlpha,  yAlpha,  acAlpha)
c
c.... form the LHS matrices, the residual vector (at alpha)
c
      call ElmGMRSclr(yAlpha,acAlpha,    x,
     &             shp,       shgl,       iBC,       
     &             BC,        shpb,       shglb,
     &             res,       iper,       ilwork,   
     &             rowp,      colm,       lhsS,
     &             cfl )


! ################################################################
! MB, call the svLS solver
      IF (svLSFlag .EQ. 1) THEN 

            lesId   = numeqns(1+nsolt+isclr)
      !        if (myrank.eq.master) write(*,*) 'lesID = ', lesId
            incL = 1
            dof = 1   ! Should be equal to 1 ? (was 4 for N.-S.)
            IF (.NOT.ALLOCATED(Res4)) THEN
                  ALLOCATE (Res4(dof,nshg), Val4(dof*dof,nnz_tot))
            END IF

            DO i=1, nshg
                  Res4(1:dof,i) = res(i,1)
            END DO

            DO i=1, nnz_tot
                  Val4(1,i)   = lhsS(i)
            END DO

      !      if (lesId.eq.2) then   ! Temperature
      !         CALL svLS_SOLVE(svLS_lhsT, svLS_sc, dof, Res4, Val4)
      !      else   ! Level set 
                  CALL svLS_SOLVE(svLS_lhs_sc, svLS_sc, dof, Res4, Val4)
      !      end if

            DO i=1, nshg
                  solinc(i,1:dof) = Res4(1:dof,i)
            END DO
      ELSE  ! not using svLS
     

! ################################################################

c
c.... lesSolve : main matrix solver
c
      lesId   = numeqns(1+nsolt+isclr)
      eqnType = 2
c
c.... setup the linear algebra solver
c
      call usrNew ( usr,        eqnType,          apermS,
     &              atempS,     res,              solinc,          
     &              flowDiag,   sclrDiag,         lesP,   
     &              lesQ,       iBC,              BC,
     &              iper,       ilwork,           numpe,
     &              nshg,       nshl,             nPermDimsS,  
     &              nTmpDimsS,  rowp,             colm,     
     &              rlhsK,      rlhsP,            lhsS,      
     &              nnz_tot )
c
c.... solve linear system
c
      call myfLesSolve ( lesId, usr )
      call getSol ( usr, solinc )

      if (numpe > 1) then
         call commu ( solinc, ilwork, 1, 'out')
      endif

      END IF      ! MB, end svLS / lesLIB choice condition      
      
      nsolsc=5+isclr
      call rstaticSclr (res, y, solinc, nsolsc) ! output scalar stats
c     
c.... end
c     
      return
      end





