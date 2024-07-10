!-------------------------------------------------------------------------------
!
!	This file contains several subrountines which are related to bubble boi-
!	ling phenomenon including bubble evaporation capability and condensation
!	capability.
!	                                                                        
!	Mengnan Li,                                                Spring, 2016
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!	Name                    Description
!       dVolume   :             Volume change due to boiling or condensation
!       numshell  :             numbers of element in bubble info collection   
!                               shell
!                               1. shell outside the bubble
!                               2. shell inside  the bubble
!       bubble_vol:             volume of each bubble
!       B_factor  :             Factor controls bubble growth rate from
!                               analytical solution
!       Tempb     :             temperature in each element
!       gytemp    :             temperature gradient in each element
!
!
!
!-------------------------------------------------------------------------------

        subroutine Bubvolgenera(dVolume, shell_num,bml,sclr_ls,elemvol_local) 
!-------------------------------------------------------------------------------
!
!	This subroutine is calculated the volume generation term during boiling
!	and condensation process.(caleld by e3res.f)
!
!-------------------------------------------------------------------------------
        use spat_var_eps ! for spatially varying epsilon_ls	
!        use bub_track
!        use bubboil_info       
        include "common.h"

        dimension dVolume(npro), elem_shell_num(4), bubvol(4), R(4)
        dimension R_1(npro), B_factor(4), shell_num(npro),sclr_ls(npro)
        dimension bml(npro,nshl,1), bubdVolume(npro,4),elemvol_local(ibksiz)
        dimension R_0(npro), shell_num_old(npro), shell_num_new(npro)

        integer i, j
	real*8 Rho_l,     Rho_v,     cp_l,       k_l,            elem_shell_num
	real*8 R,         a_l,       B_factor,   dVolume,        epsilon_lsd_tmp
	real*8 R_1,       bubvol,    sclr_ls,    elemvol_local,  R_0
        real*8 source

        Rho_l = datmat(1,1,1)  !958
        Rho_v = datmat(1,1,2) !0.579
        cp_l = datmat(1,3,1)!1.22
        k_l = datmat(1,4,1)!0.679
        dVolume(:) = 1.0e-10
        shell_num(:)=1.0e-10

        do i=1, i_num_bubbles
        if (lstep+1 .eq. int(1+inistep))then
          elem_shell_num(i) = 2000
          bubvol(i)=(4.0/3.0)*pi*((2E-4)**3.0)
          bubble_tempG(i)=1.0E-12
!          avg_info(1,19)=1.0E-15
        else
          elem_shell_num(i) = numshell(i,2)
          bubvol(i)=bubble_vol(i)
        endif  
        R(i) =(((3.0E0/4.0E0)*bubvol(i))/pi)**(1.0E0/3.0E0)
!        if(myrank.eq.master)then       
!         write(*,*)'source', source
!        endif
        if(bubboil.eq.0.and.bubgrow.eq.1)then
!         if(bubgrow.eq.1)then  ! for test only
          a_l = k_l/(Rho_l*cp_l*1.0e3) !1.0e3
          B_factor(i) = ((12.0E0*a_l/pi)**(0.5E0))* 
     &               ((delt_T(i)*cp_l*1.0e3*Rho_l)/(h_fg*Rho_v))
!        if (myrank.eq.master)write(*,*)'delt_T', delt_T(1,i)
        do j = 1, npro
          bubdVolume(j,i) = 2.0E0*pi*R(i)*(B_factor(i)**(2.0E0))
!     &        +source/npro
!        if(myrank.eq.master)write(*,*)'i,dVolume',i,dVolume(j,i)
        enddo
        endif

        if(bubboil.eq.1)then  ! loop over elements (at the local level)
        do j = 1, npro 
          epsilon_lsd_tmp = 
     &      elem_local_size(lcblk(1,iblk)+j-1)
!          epsilon_lsd_tmp = epsilon_lsd*
!     &         (elemvol_local(j)**(1.0/3.0))
!          epsilon_lsd_tmp = 
!     &         (elemvol_local(j)**(1.0/3.0))
!          R_1(i)=R(i)+epsilonBT*1.0   ! sqrt(3.0) for structure mesh
                                                    ! sqrt(2.0) for parasolid mesh
!        write(*,*)'elemvol_local',elemvol_local(j)   
        bubdVolume(j,i) = bubble_tempG(i)*k_l*(1.0/Rho_v-1.0/Rho_l)
     &                    /h_fg
        bubdVolume(j,i) = bubdVolume(j,i)*4.0d0*pi*(R(i)**(2.0))
!         bubdVolume(j,i) = bubdVolume(j,i)
!        if (lstep+1 .lt. int(40+inistep))then
!        bubdVolume(j,i) = 0.0d0
!        endif
!        bubdVolume(j,i) = bubdVolume(j,i)/(epsilonBT)
!        bubdVolume(j,i) =  bubdVolume(j,i)/sqrt(epsilonBT)
!        if(myrank.eq.master) write(*,*)'R(i)',R(i)
!         bubdVolume(j,i) = bubble_tempG(i)*k_l/(h_fg*Rho_v)/epsilonBT   
        enddo
        endif

        enddo ! i_num_bubbles


!       Assembly the source matrix
        do i = 1, i_num_bubbles
           do j = 1, npro
          if ((sclr_ls(j).GT.-2.0E0*epsilonBT).and.
     &    (sclr_ls(j).LT.-1.0E0*epsilonBT)) then
              do n = 1, nshl
                if(INT(bml(j,n,1)).eq.i) then
                 dVolume(j) = bubdVolume(j,i)
                 shell_num(j) = elem_shell_num(i)
                endif
              enddo
           else
                 dVolume(j) = 0.0d0
           endif
!           if(isnan(dVolume(j))) dVolume(j) = 0.0d0
         
!        if(dVolume(j).ne.0.0)then
!         write(*,*)'dVolume',dVolume(j)
!        endif
          enddo
!          do k = 1, npro
!            if (lstep+1 .eq. int(1+inistep))then
!            shell_num_old(k)=1.75E+4
!            endif
!            shell_num(k) = 0.8*shell_num_new(k)+ shell_num_old(k)
!            shell_num_old(k)= shell_num(k)
!          enddo
!            if (lstep+1 .eq. int(3+inistep))then
!            write(*,*)shell_num_new(10),shell_num_old(10)
!            endif
        enddo
        return
        end
!-------------------------------------------------------------------------------
c===============================================================================
c===============================================================================

        subroutine Bubheatflux(yl, shpfun, shg, elemvol_local,Tempb, gytemp,bml)
!-------------------------------------------------------------------------------
!
!	This subroutine is used to calculated total heat flux flowing into bubb-
!	-le due to temperature difference(called in e3ivar.f)
!
!-------------------------------------------------------------------------------
        use  spat_var_eps ! for spatially varying epsilon_ls
!        use  bubboil_info
!        use  bub_track
        include "common.h"
        integer    i,  n
        dimension  yl(npro,nshl,ndof),       shpfun(npro,nshl),
     &             shg(npro,nshl,nsd),       Temp(npro),
     &             elemvol_local(ibksiz),    
     &             gyti(npro,nsd),
     &             Sclr(npro),               Tempb(npro),
     &             gytemp(npro,5),           R_0(npro),
     &             R_1(npro)
        dimension  bml(npro,nshl,1)
	real*8     elemvol_local,            elemvol            
        real*8     epsilon_lsd_tmp
	real*8     epsilon_ls_tmp,           gytemp,              
     &             R_0,                      R_1,
     &             num,
     &             local_volume,             bml,
     &             RR
        Sclr = zero
        Tempb = zero
        gyti = zero
        gytemp(:,:) = zero
!         write(*,*) 'bubble_vol',bubble_vol(1)
        do i = 1, npro
          do n = 1, nshl
            Sclr(i) = Sclr(i) + shpfun(i,n) * yl(i,n,6) !scalar
c     
c!     .... compute the global gradient of Scalar variable
c     
            gyti(i,1) = gyti(i,1) + shg(i,n,1) * yl(i,n,6) 
            gyti(i,2) = gyti(i,2) + shg(i,n,2) * yl(i,n,6)
            gyti(i,3) = gyti(i,3) + shg(i,n,3) * yl(i,n,6)
c     
           enddo
         enddo
c
c!  .... compute the global gradient of Temperature outside bubble 2 epsilon
c
       
        do i=1, npro
         do k=1, i_num_bubbles
          if (lstep+1 .eq. int(1+inistep))then
            bubble_vol(k)=(4.0/3.0)*pi*((2.0E-4)**3)
          endif
             do n = 1, nshl
                if(INT(bml(i,n,1)).eq.k) then
                R_0(i)=(((3.0E0/4.0E0)*bubble_vol(k))/pi)**(1.0E0/3.0E0)
!                R_1(i)=R_0(i)+epsilon_ls_tmp*2.5/sqrt(3.0)  
!                R_1(i)=R_0(i)+epsilonBT*2.5
                endif
              enddo
            RR = maxval(R_0(:))
!          if(myrank.eq.master)write(*,*)'RR',RR
          enddo
         enddo ! i_num_bubbles
         elemvol = 0.0d0
         do i=1, npro
           if (Sclr(i).le. 1.0E0*epsilonBT
     &       .and. Sclr(i).ge. 0.0E0*epsilonBT) then
           elemvol=elemvol+elemvol_local(i)
           endif
         enddo
!           R_1=0.0
!           R_0=(((3.0E0/4.0E0)*bubble_vol(k))/pi)**(1.0E0/3.0E0)
!           R_1=R_0+epsilon_lsd_tmp*2.5/sqrt(2.0)
!         if(myrank.eq.master) write(*,*)'elemvol',elemvol
!          write(*,*) 'bubble_vol,R_0',bubble_vol(1),R_0
!      .... R_0 is zero in the begining of the run due to the sequence issue
!        if (lstep+1 .gt. int(10+inistep))then
        do i = 1, npro
              do n = 1, nshl
                   Tempb(i) = Tempb(i) + shpfun(i,n) * yl(i,n,5)
                   R_1(i)=RR + shpfun(i,n)*yl(i,n,6)
              enddo
!           if (Sclr(i).le. 3.0E0*epsilon_ls_tmp
!     &       .and. Sclr(i).ge. 2.0E0*epsilon_ls_tmp) then
           if (Sclr(i).le. 1.0E0*epsilonBT
     &       .and. Sclr(i).ge. 0.0E0*epsilonBT) then

              do n = 1, nshl
!                   Tempb(i) = Tempb(i) + shpfun(i,n) * yl(i,n,5) !Temperature
                   gytemp(i,1) = gytemp(i,1) + shg(i,n,1) * yl(i,n,5)
                   gytemp(i,2) = gytemp(i,2) + shg(i,n,2) * yl(i,n,5)
                   gytemp(i,3) = gytemp(i,3) + shg(i,n,3) * yl(i,n,5)
             enddo ! for nshl
           gytemp(i,4) = (gytemp(i,1)*gyti(i,1)+gytemp(i,2)*gyti(i,2)
     &   +gytemp(i,3)*gyti(i,3))/sqrt(gyti(i,1)**2+gyti(i,2)**2+gyti(i,3)**2)
!           gytemp(i,5)=gytemp(i,4)*((R_1(i)/R_0(i))**2.0) ! *elemvol_local(i)
!            gytemp(i,5)=gytemp(i,4)*(1.375d0*1.375d0)
!            gytemp(i,5)=gytemp(i,4)
!            if(myramk.eq.master) write(*,*) 'g4',gytemp(i,4)
!            if(myrank.eq.master) write(*,*) R_1(i)/R_0(i)
!            write(*,*)'gytemp(i,5)', gytemp(i,5)
             gytemp(i,5)=gytemp(i,4)*elemvol_local(i)*((R_1(i)/RR)**2.0)/elemvol
!              gytemp(i,5)=gytemp(i,4)*elemvol_local(i)*((R_1(i)/RR)**2.0)
!             gytemp(i,5)=gytemp(i,4)*((R_1(i)/R_0(i))**2.0)
           else
            gytemp(i,5)=0.0d0
            Tempb(i)=0.0
           endif
!          if(myrank.eq.master)then
!          endif
        enddo
!        else
!          gytemp(:,5)=0.0d0
!        endif
!        enddo
        return
        end

!-------------------------------------------------------------------------------
c===============================================================================
c===============================================================================

!        subroutine Bubibc(y,bml)
!-------------------------------------------------------------------------------
!
!	This subroutine is used to set up the boudary condition inside the bubb-
!	le(called in itrbc.f)
!
!-------------------------------------------------------------------------------
!        include "common.h"
!        dimension y(nshg,ndof)
!        dimension bml(npro,nshl,1)
!        real*8  y, bml
        
!        if(isclr.eq.0) then
!        do i=1, npro
!             do n = 1, nshl
!                if(INT(bml(i,n,1)).eq.1) then
!                  do j=1, nshg
!                  if (y(j,6).LE.0.0) then
!                  y(j,id) = 373.15 ! Kelvin
!                  endif
!                  enddo
!                else
!                  do j=1, nshg
!                  if (y(j,6).LE.0.0) then
!                  y(j,id) = 370.65 ! Kelvin
!                  endif
!                  enddo
!                endif
!             enddo
!        enddo
!        endif

!        return
!        end

