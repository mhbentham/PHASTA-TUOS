c-----------------------------------------------------------------------
c
c    Initialize the predictor multicorrector (set up parameters)
c
c    Modified by Alberto Figueroa.  Winter 2004
c-----------------------------------------------------------------------
      subroutine itrSetup ( y,  acold ) 
      
      include "common.h"
      
      real*8     y(nshg,ndof),     acold(nshg,ndof)
      
c
c  Define the Hulbert parameters
c  second order if between (and including) 0 and 1 otherwise backward Euler
c
      if( rhoinf(itseq).lt.0.or.rhoinf(itseq).gt.1) then ! backward Euler
         almi   = one
         alfi   = one
         gami   = one
         ipred  = 1
      else           !second order family
         almi   = (three-rhoinf(itseq))/(one+rhoinf(itseq))/two
         alfi   = one/(one+rhoinf(itseq))
         gami   = pt5+almi-alfi
         if(ideformwall.eq.1) then
            betai=1.0
         else
            betai=0.0
         endif
      endif
c     
c.... set the jacobian type
c     
      Jactyp=0
      if (impl(itseq) .eq. 3) then
         Jactyp = 1
         impl(itseq) = 2
      endif
c     
c.... same_Dy predictor special case
c     
      if (ipred.eq.4 .and. itseq .eq. 1 ) then
         y=y-(one-alfi)*Delt(1)*acold
         if ( rhoinf(itseq) .eq. 0.0 ) then
            ipred = 3
         endif
      endif
c
c.... set the global time increment and the CFL data
c
      Dtgl   = one / Delt(itseq)  ! caution: inverse of time step
      CFLfld = CFLfl(itseq)
      CFLsld = CFLsl(itseq)
      
      return
      end


c-----------------------------------------------------------------------
c
c    Predict solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrPredict (yold,  y,  acold,   ac,   uold,   u)
      
      include "common.h"
      
      real*8        y(nshg,ndof),               ac(nshg,ndof),
     &              u(nshg,nsd),                yold(nshg,ndof),
     &              acold(nshg,ndof),           uold(nshg,nsd)

c Predict Flow Variables
      
      if ( ipred.eq.1) then     ! yn+1_pred=yn
         fct = (gami-one)/gami
         y(:,1:5)  = yold(:,1:5)
         ac(:,1:5) = acold(:,1:5) * fct
         if(ideformwall.eq.1) 
     &          u(:,1:3) = uold(:,1:3) + Delt(itseq)*yold(:,1:3) + 
     &              pt5*((gami-two*betai)/gami)*
     &              Delt(itseq)*Delt(itseq)*acold(:,1:3)
      endif

c     
      if ( ipred.eq.2) then     ! an+1_pred=0
         y(:,1:5)  = yold(:,1:5) + (one - gami)/Dtgl * acold(:,1:5)
         ac(:,1:5) = 0.0
      endif

c     
      if(ipred.eq.3 ) then      ! an+1_pred=an
         y(:,1:5)  = yold(:,1:5)+alfi*Delt(itseq)*acold(:,1:5)
         ac(:,1:5) = acold(:,1:5)
      endif
c             do i=1,nshg
c             if(isnan(y(i,5)))write(*,*)'y'
c             if(isnan(yold(i,5)))write(*,*)'yold'
c             if(isnan(alfi))write(*,*)'alfi'
c             enddo

c     
      if ( ipred.eq.4 ) then    ! protect from DC=4 rho=0, same dV
         fct1 = alfi/(one-alfi)
         fct2 = one-almi/gami
         fct3 = almi/gami/alfi*Dtgl
         y(:,1:5)    = yold(:,1:5)+fct1*(yold(:,1:5)-y(:,1:5))
         ac(:,1:5)   = acold(:,1:5)*fct2+(y(:,1:5)-yold(:,1:5))*fct3
      endif

!             do i=1,nshg
!             if(isnan(y(i,5)))write(*,*)'y'
!             if(isnan(yold(i,5)))write(*,*)'yold'
!             if(isnan(alfi))write(*,*)'alfi'
!             enddo


c
c Predict Level Set Variables
c
      if (ndof.gt.5 .and. ilset.eq.2) then
c Scalar 1
         if (iSolvLSSclr1.eq.1) then
            y(:,6)  = yold(:,6)
            ac(:,6) = acold(:,6)
            ac(:,6) = zero
         elseif ( ipred.eq.1) then     ! yn+1_pred=yn
            fct = (gami-one)/gami
            y(:,6)  = yold(:,6)
            ac(:,6) = acold(:,6) * fct
         elseif ( ipred.eq.2) then     ! an+1_pred=0
            y(:,6)  = yold(:,6) + (one - gami)/Dtgl * acold(:,6)
            ac(:,6) = 0.0
         elseif (ipred.eq.3 ) then     ! an+1_pred=an
            y(:,6)  = yold(:,6)+alfi*Delt(itseq)*acold(:,6)
            ac(:,6) = acold(:,6)
         elseif ( ipred.eq.4 ) then    ! protect from DC=4 rho=0, same dV
            fct1 = alfi/(one-alfi)
            fct2 = one-almi/gami
            fct3 = almi/gami/alfi*Dtgl
            y(:,6)    = yold(:,6)+fct1*(yold(:,6)-y(:,6))
            ac(:,6)   = acold(:,6)*fct2+(y(:,6)-yold(:,6))*fct3
         endif
c Scalar 2
         if ((iSolvLSSclr2.eq.1) .or. (iSolvLSSclr2.eq.2)) then
            y(:,7)  = yold(:,7)
            ac(:,7) = acold(:,7)
         elseif ( ipred.eq.1) then     ! yn+1_pred=yn
            fct = (gami-one)/gami
            y(:,7)  = yold(:,7)
            ac(:,7) = acold(:,7) * fct
         elseif ( ipred.eq.2) then     ! an+1_pred=0
            y(:,7)  = yold(:,7) + (one - gami)/Dtgl * acold(:,7)
            ac(:,7) = 0.0
         elseif (ipred.eq.3 ) then     ! an+1_pred=an
            y(:,7)  = yold(:,7)+alfi*Delt(itseq)*acold(:,7)
            ac(:,7) = acold(:,7)
         elseif ( ipred.eq.4 ) then    ! protect from DC=4 rho=0, same dV
            fct1 = alfi/(one-alfi)
            fct2 = one-almi/gami
            fct3 = almi/gami/alfi*Dtgl
            y(:,7)    = yold(:,7)+fct1*(yold(:,7)-y(:,7))
            ac(:,7)   = acold(:,7)*fct2+(y(:,7)-yold(:,7))*fct3
         endif
      endif              
c     

      return
      end

c-----------------------------------------------------------------------
c
c    Correct solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrCorrect ( y,     ac,   u,   solinc )
      
      include "common.h"
      
      real*8      y(nshg,ndof),               ac(nshg,ndof),  
     &            u(nshg,nsd),                solinc(nshg,4)
      
      fct1 = gami*Delt(itseq)
      fct2 = gami*alfi*Delt(itseq)
      
      y(:,1:3)  = y(:,1:3)  + fct1 * solinc(:,1:3)
      y(:,4  )  = y(:,4  )  + fct2 * solinc(:,4  )
      ac(:,1:3) = ac(:,1:3) + solinc(:,1:3)
      if(ideformwall.eq.1) 
     &   u(:,1:3)  = u(:,1:3)  + 
     &            Delt(itseq)*Delt(itseq)*betai*solinc(:,1:3)
c     
      return
      end

c-----------------------------------------------------------------------
c
c    Correct solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrCorrectSclr ( y,     ac,   solinc )
      
      include "common.h"
      
      real*8      y(nshg,ndof),               ac(nshg,ndof),  
     &            solinc(nshg)
      
      is=5+isclr
      if (((isclr.eq.1) .and.(iSolvLSSclr1.eq.1)) .or.
     &    ((isclr.eq.2) .and.(iSolvLSSclr2.eq.1)) .or.
     &    ((isclr.eq.2) .and.(iSolvLSSclr2.eq.2))) then
        fct1 = Delt(itseq)
      y(:,is)  = y(:,is)  + fct1 * solinc(:)
      ac(:,is) = ac(:,is) + solinc(:)
      else
        fct1 = gami*Delt(itseq)
        y(:,is)  = y(:,is)  + fct1 * solinc(:)
        ac(:,is) = ac(:,is) + solinc(:)
      endif
c     
      return
      end

c-----------------------------------------------------------------------
c
c    Correct solution at time n+1, protecting against negative values
c
c-----------------------------------------------------------------------
      subroutine itrCorrectSclrPos ( y,     ac,   solinc )
      
      include "common.h"
      
      real*8      y(nshg,ndof),               ac(nshg,ndof),  
     &            solinc(nshg),               posinc(nshg)
      
      fct1 = gami*Delt(itseq)
      if(fct1.eq.0) then
         call itrCorrectSclr( y, ac, solinc)
      else
         turbUpdFct = 1.0
         updFct = turbUpdFct
         updFct2 = 1.0
c         if(topHdTscFlag.and.PRM_TCS_PARTIAL_UPDATE)then
c            updFct2 = max(MTH_SQRT_EPS,min(topHdUpdFct,1.0))
c         endif
         fctCr = fct1*updFct*updFct2
         c0 = 0.01
         c1 = -(1.0-c0)/fctCr
         is=5+isclr
c         if(any(updFct*updFct2*solinc(:).lt.c1*y(:,is))) then
c            write(*,*) 'INSTEAD OF GETTING NEGATIVE FIELD'
c            write(*,*) 'BROUGHT FIELD DOWN TO 1 PERCENT'
c            write(*,*) 'FOR SCALAR NUMBER ',isclr
c            write(*,*) '(SEE itrCorrectSclr in itrPC.f)'
c         endif
         posinc = max(updFct*updFct2*solinc,c1*y(:,is))
         y(:,is)  = y(:,is)  + fct1 * posinc(:)
         ac(:,is) = ac(:,is) + posinc(:)
      endif
c     
      return
      end


c-----------------------------------------------------------------------
c
c    Compute solution and acceleration at n+alpha
c
c-----------------------------------------------------------------------
      subroutine itrYAlpha ( uold,        yold,        acold,        
     &                       u,           y,           ac,
     &                       uAlpha,      yAlpha,      acAlpha )

c      use readarrays       !reads in uold and acold     
      include "common.h"
      
      real*8        yold(nshg,ndof),            acold(nshg,ndof),
     &              y(nshg,ndof),               ac(nshg,ndof),
     &              yAlpha(nshg,ndof),          acAlpha(nshg,ndof),
     &              u(nshg,nsd),                uold(nshg,nsd),
     &              uAlpha(nshg,nsd)



      acAlpha(:,4) = zero  !pressure acceleration is never used but....

      yAlpha (:,1:3) = yold(:,1:3) 
     &                  + alfi * (y(:,1:3) - yold(:,1:3))

      acAlpha(:,1:3) = acold(:,1:3)
     &                  + almi * (ac(:,1:3) - acold(:,1:3))

      yAlpha (:,4  ) = y(:,4)

      if(ideformwall.eq.1) uAlpha (:,1:3) = uold(:,1:3) 
     &                  + alfi * (u(:,1:3) - uold(:,1:3))
      
      if(ndof.ge.5) then
c
c  Now take care of temperature, turbulence, what have you
c
      

         yAlpha (:,5:ndof  ) = yold(:,5:ndof) 
     &                       + alfi * (y(:,5:ndof) - yold(:,5:ndof))
         acAlpha(:,5:ndof  ) = acold(:,5:ndof) 
     &                       + almi * (ac(:,5:ndof) - acold(:,5:ndof))
     
       endif
      return
      end

c-----------------------------------------------------------------------
c
c    Update solution at end of time step
c
c-----------------------------------------------------------------------
      subroutine itrUpdate( yold,          acold,        uold,
     &                      y,             ac,           u )

c      use readarrays            !reads in uold and acold

      include "common.h"
      
      real*8        yold(nshg,ndof),            acold(nshg,ndof),
     &              y(nshg,ndof),               ac(nshg,ndof),
     &              u(nshg,nsd),                uold(nshg,nsd)

            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
!	Allow for redistance of the level set function, AD 5/8/00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      if we are re-distancing the levelset function the we need 
!      to update the levelset variable with the corrected value
!      stored in the "next" scalar position.
            
        if(iLSet .eq.2 .and. isclr .eq. 2) then 
     
           y(:,6) = y(:,7) ! will need to also fix the yold ac etc.?????
     	   ac(:,7)=zero	
        endif
        

      
      yold  = y
      acold = ac
      if(ideformwall.eq.1)  uold  = u
      
      return
      end


c-----------------------------------------------------------------------
c
c    Update solution at end of time step
c
c-----------------------------------------------------------------------
      subroutine itrUpdateDist( yold,          acold,
     &                          y,             ac )

      include "common.h"
      
      real*8        yold(nshg,ndof),            acold(nshg,ndof),
     &              y(nshg,ndof),               ac(nshg,ndof)

c update the old value with the current solutions at n+1
c
      ac(:,7)    = zero
      yold(:,7)  =  y(:,7)
      acold(:,7) = ac(:,7)
      
      return
      end
