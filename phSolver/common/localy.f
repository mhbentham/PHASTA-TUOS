        subroutine localy (global, rlocal, ientmp, n, code)
c
c----------------------------------------------------------------------
c
c This subroutine performs a vector gather/scatter operation.
c
c input:
c  global (nshg,n)              : global array
c  rlocal (npro,nshl,n)         : local array
c  ien    (npro,nshl)           : nodal connectivity
c  n                            : number of d.o.f.''s to be copied
c  code                         : the transfer code
c                                  .eq. 'gather  ', from global to local
c                                  .eq. 'scatter ', add  local to global 
c                                  .eq. 'globaliz', from local to global
c
c
c Zdenek Johan, Winter 1992.
c----------------------------------------------------------------------
c
        include "common.h"

        dimension global(nshg,n),           rlocal(npro,nshl,n),
     &            ien(npro,nshl),           ientmp(npro,nshl)
c
        character*8 code
        
c
c.... cubic basis has negatives in ien
c
        if (ipord > 2) then
           ien = abs(ientmp)
        else
           ien = ientmp
        endif
c
c.... ------------------------>  'localization  '  <--------------------
c
        if (code .eq. 'gather  ') then
c
c.... set timer
c
c          call timer ('Gather  ')
c
c.... gather the data to the current block 
c

CAD      rlocal = yl={P, u, v, w, T, scalar1, ...}
CAD	 global = y = {u, v, w, P, T, scalar1, ...}

CAD      Put u,v,w in the slots 2,3,4 of yl 

          do j = 1, 3
            do i = 1, nshl
              rlocal(:,i,j+1) = global(ien(:,i),j)
            enddo
          enddo

CAD      Put Pressure in the first slot of yl

          do i = 1, nshl
             rlocal(:,i,1) = global(ien(:,i),4)
          enddo

CAD      Fill in the remaining slots with T, and additional scalars
          
          if(n.gt.4) then
             do j = 5, n
                do i = 1, nshl
                   rlocal(:,i,j) = global(ien(:,i),j)
                enddo
             enddo
          endif
c
c.... transfer count
c
          gbytes = gbytes + n*nshl*npro
c
c.... return
c
c          call timer ('Back    ')
          return
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        if (code .eq. 'scatter ') then
           write(*,*) 'do not use localy here'
        endif
c
c.... ------------------------->  'globalizing '  <----------------------
c
        if (code .eq. 'globaliz') then
           write(*,*) 'do not use localy here'
        endif
c
c.... --------------------------->  error  <---------------------------
c
        call error ('local   ', code, 0)
c
c.... end
c
        end
        subroutine localyotwn (global, rlocal, ientmp, n, code)
c
c----------------------------------------------------------------------
c
c This subroutine performs a vector gather/scatter operation.
c The local array is assigned values from the Off-The-Wall Nodes (otwn)
c
c input:
c  global (nshg,n)              : global array
c  rlocal (npro,nshl,n)         : local array
c  ien    (npro,nshl)           : nodal connectivity
c  n                            : number of d.o.f.'s to be copied
c  code                         : the transfer code
c                                  .eq. 'gather  ', from global to local
c                                  .eq. 'scatter ', add  local to global 
c                                  .eq. 'globaliz', from local to global
c
c
c Zdenek Johan, Winter 1992.
c----------------------------------------------------------------------
c
        use turbSA		! need otwn(nshg), otwn_vec(nshg,3)
        include "common.h"

        dimension global(nshg,n),           rlocal(npro,nshl,n),
     &            ien(npro,nshl),           ientmp(npro,nshl)
c
        character*8 code
	real*8 rho_gh(nshg)
        
c
c.... cubic basis has negatives in ien
c
        if (ipord > 2) then
           ien = abs(ientmp)
        else
           ien = ientmp
        endif
c
c.... ------------------------>  'localization  '  <--------------------
c
        if (code .eq. 'gather  ') then
c
c.... gather the data to the current block 
c

CAD      rlocal = ylotwn={P, u, v, w, T, scalar1, ...}
CAD	 global = y =    {u, v, w, P, T, scalar1, ...}

CAD      Put u,v,w in the slots 2,3,4 of ylotwn 

          do j = 1, 3
            do i = 1, nshl
              rlocal(:,i,j+1) = global(otwn(ien(:,i)),j)
            enddo
          enddo

CAD      Put Pressure in the first slot of ylotwn
c For convective pressure boundary condition then we must
c correct for the hydrostatic head that may exist between the 
c wall and off-the-wall nodes.  Note localyotwn is only called
c for any(btest(iBCB(:,1)).
c
c First, compute the pressure correction terms (rho_gh = 
c density*gravity).  otwn_vec(nshg,:) gives the distance vector
c to the otw node while datmat(:,5,1) gives the body force vector.
c
          rho_gh = zero
	  do j = 1, 3
	    rho_gh(:) = rho_gh(:) + otwn_vec(:,j)*datmat(j,5,1)
	  enddo
c
c now multiply by density to get rho*g*height
c
	  do i = 1, nshg
	    if (iLSet .eq. 0) then
	      rho_gh(i) = rho_gh(i) * datmat(1,1,1)
	    else
	      isc = abs(iRANS)+6
	      if (global(i,isc) .ge. 0.0) then
	        rho_gh(i) = rho_gh(i) * datmat(1,1,1)
	      else
	        rho_gh(i) = rho_gh(i) * datmat(1,1,2)
	      endif
	    endif
	  enddo  
c
          do i = 1, nshl
             rlocal(:,i,1) = global(otwn(ien(:,i)),4) - 
     &                       rho_gh(ien(:,i))
          enddo

CAD      Fill in the remaining slots with T, and additional scalars
          
          if(n.gt.4) then
             do j = 5, n
                do i = 1, nshl
                   rlocal(:,i,j) = global(otwn(ien(:,i)),j)
                enddo
             enddo
          endif
c
c.... transfer count
c
          gbytes = gbytes + n*nshl*npro
c
c.... return
c
c          call timer ('Back    ')
          return
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        if (code .eq. 'scatter ') then
           write(*,*) 'do not use localyotwn here'
        endif
c
c.... ------------------------->  'globalizing '  <----------------------
c
        if (code .eq. 'globaliz') then
           write(*,*) 'do not use localyotwn here'
        endif
c
c.... --------------------------->  error  <---------------------------
c
        call error ('localyotwn   ', code, 0)
c
c.... end
c
        end
