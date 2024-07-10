      subroutine tnanq (u, n, arrname)

      include "common.h"

      dimension   u(nshg,n),xmax(1000),xmin(1000)
      dimension   ucopy(127)
      character*8 arrname

      nnanq = 0
      DO j = 1,n
        xmax(j) = -1.d32
        xmin(j) = 1.d32
	DO i = 1,nshg
	  IF (u(i,j) .ne. u(i,j)) nnanq = nnanq + 1
          xmax(j)=max(xmax(j),u(i,j))
          xmin(j)=min(xmin(j),u(i,j))
	ENDDO
      ENDDO
      ucopy(1:127)=u(1:127,1)

      IF (nnanq .ne. 0) call error('tnanq   ',arrname,nnanq)
c     return
      if (myrank .eq. master) then
      write(*,*) arrname
      jskip=n/5
      if(jskip.lt.1) then
        write(*,110) (xmax(i),i=1,n)
        write(*,110) (xmin(i),i=1,n)
      else
      do i=1,jskip
        write(*,110) (xmax(i+j*jskip),j=0,4)
      enddo
      if(jskip*5.lt.n) write(*,110) xmax(jskip*5+1:n) !leftovers
      do i=1,jskip
        write(*,110) (xmin(i+j*jskip),j=0,4)
      enddo
      if(jskip*5.lt.n) write(*,110) xmin(jskip*5+1:n) !leftovers
      endif
      endif
110   format(9(e14.7))

      return
      end
