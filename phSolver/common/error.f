        subroutine error (routin, variab, num)
c
c----------------------------------------------------------------------
c
c This utility routine prints out the error and stops the program.
c
c input:
c  routin       : name of the routine where the error occurred
c  variab       : an 8-character error message
c  num          : any integer number associated with the error
c
c Farzin Shakib, Summer 1985.
c----------------------------------------------------------------------
c
        include "common.h"
        include "mpif.h"
c
        character*8 routin, variab
c
        data ierchk /0/
c
c.... check for redundant error
c
        if (ierchk .eq. 1) stop
        ierchk = 1
c
c.... open file
c
        open (unit=ierror, file=ferror, status='unknown')
c
c.... print the error
c
	if (myrank.eq.master) then
        write (*,1000) title, routin, variab, num
        if (num .ne. 0) write (ierror,1000) title, routin, variab, num
        if (num .eq. 0) write (ierror,1000) title, routin, variab
	end if
c
c.... halt the process
c
        close (ierror)

	if (myrank.eq.master) then
        WRITE(6,'(A,G14.6)') 'Life: ',death - birth
	end if
        if (numpe > 1) then
           call MPI_ABORT(MPI_COMM_WORLD)
        endif
        
 
1000    format(' ',a80,//,
     &         ' ****** Error occurred in routine <',a8,'>',/,
     &          '  Error code :',a8,:,' : ',i8,//)
        end
