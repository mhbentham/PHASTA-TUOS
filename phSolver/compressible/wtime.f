      subroutine wtime

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      character*5 cname

      dimension timeg(100,numpe)

      if (numpe > 1) then
      CALL MPI_GATHER (ttim, 100, MPI_DOUBLE_PRECISION, timeg, 100, 
     &                 MPI_DOUBLE_PRECISION, master, 
     &                 MPI_COMM_WORLD, ierr)
      endif
      if (myrank .eq. master) then

      ftime = trim(ftime) // cname(lstep)

      open (itime, file=ftime, form='formatted')

      write (itime,101) 'Tcore cpu/wc: ',ttim(1), ttim(2)
      write (itime,100) 'local       : ',ttim(3) ,ttim(3) /ttim(2)
      write (itime,100) 'global      : ',ttim(4) ,ttim(4) /ttim(2)
      write (itime,100) 'communic.   : ',ttim(5) ,ttim(5) /ttim(2)
      write (itime,100)
      write (itime,100) 'e3q         : ',ttim(33),ttim(33)/ttim(2)
      write (itime,100) '  e3qvar    : ',ttim(34),ttim(34)/ttim(2)
      write (itime,100) 'e3 outside  : ',ttim(31),ttim(31)/ttim(2)
      write (itime,100) 'e3          : ',ttim(6) ,ttim(6) /ttim(2)
      write (itime,100) 'e3ivar out  : ',ttim(8), ttim(8) /ttim(2)
      write (itime,100) 'e3ivar      : ',ttim(20),ttim(20)/ttim(2)
      write (itime,100) '  dui       : ',ttim(10),ttim(10)/ttim(2)
      write (itime,100) '  prim vars : ',ttim(11),ttim(11)/ttim(2)
      write (itime,100) '  acc       : ',ttim(12),ttim(12)/ttim(2)
      write (itime,100) '  therm prop: ',ttim(27),ttim(27)/ttim(2)
      write (itime,100) '  dxidx     : ',ttim(26),ttim(26)/ttim(2)
      write (itime,100) '  shape grad: ',ttim(13),ttim(13)/ttim(2)
      write (itime,100) '  grad Y    : ',ttim(7), ttim(7) /ttim(2)
      write (itime,100) '  div (qi)  : ',ttim(32),ttim(32)/ttim(2)
      write (itime,100) 'e3mtrx out  : ',ttim(9) ,ttim(9) /ttim(2)
      write (itime,100) 'e3mtrx      : ',ttim(21),ttim(21)/ttim(2)
      write (itime,100) 'e3conv out  : ',ttim(14),ttim(14)/ttim(2)
      write (itime,100) 'e3conv      : ',ttim(22),ttim(22)/ttim(2)
      write (itime,100) 'e3visc out  : ',ttim(15),ttim(15)/ttim(2)
      write (itime,100) 'e3visc      : ',ttim(23),ttim(23)/ttim(2)
      write (itime,100) 'e3LS  out   : ',ttim(16),ttim(16)/ttim(2)
      write (itime,100) 'e3LS        : ',ttim(24),ttim(24)/ttim(2)
      write (itime,100) '  e3tau     : ',ttim(25),ttim(25)/ttim(2)
      write (itime,100) 'e3juel out  : ',ttim(17),ttim(17)/ttim(2)
      write (itime,100) 'e3juel      : ',ttim(28),ttim(28)/ttim(2)
      write (itime,100) 'e3bdg  out  : ',ttim(18),ttim(18)/ttim(2)
      write (itime,100) 'e3bdg       : ',ttim(30),ttim(30)/ttim(2)
      write (itime,100) 'e3assm out  : ',ttim(19),ttim(19)/ttim(2)
      write (itime,100) 'e3assm      : ',ttim(29),ttim(29)/ttim(2)
      write (itime,100)
      write (itime,100) 'e3b         : ',ttim(40),ttim(40)/ttim(2)
      write (itime,100) 'e3nref      : ',ttim(42),ttim(42)/ttim(2)
      write (itime,100)
      ttim(50) = ttim(51) + ttim(52) + ttim(53)
      write (itime,100) 'bc          : ',ttim(50),ttim(50)/ttim(2)
      write (itime,100) '  bc3res    : ',ttim(51),ttim(51)/ttim(2)
      write (itime,100) '  bcbdg     : ',ttim(52),ttim(52)/ttim(2)
      write (itime,100) '  itrbc     : ',ttim(53),ttim(53)/ttim(2)
      write (itime,100)
      write (itime,100) 'vy          : ',ttim(55),ttim(55)/ttim(2)
      write (itime,100) 'vq          : ',ttim(56),ttim(56)/ttim(2)
      write (itime,100)
      write (itime,100) 'getdmc      : ',ttim(60),ttim(60)/ttim(2)
      write (itime,100) 'sumgat-wait : ',ttim(61),ttim(61)/ttim(2)
      write (itime,100) 'sumgat-work : ',ttim(62),ttim(62)/ttim(2)
      write (itime,100) 'gerbar (cur): ',ttim(64),ttim(64)/ttim(2)
      write (itime,100) 'sponge      : ',ttim(66),ttim(66)/ttim(2)
      write (itime,100) 'rstat-wait  : ',ttim(67),ttim(61)/ttim(2)
      write (itime,100) 'rstat-work  : ',ttim(68),ttim(62)/ttim(2)
      write (itime,100)
      write (itime,100) 'barrier 1   : ',ttim(71),ttim(71)/ttim(2)
      write (itime,100) 'barrier 2   : ',ttim(72),ttim(72)/ttim(2)
      write (itime,100)
      write (itime,100) 'blocs elmmfg: ',ttim(80),ttim(80)/ttim(2)
      write (itime,100) 'blocs itrres: ',ttim(81),ttim(81)/ttim(2)

      write (itime,100) 
      write (itime,100) 'Temps min/max:'
      write (itime,100) 
      write (itime,101) '1 :', MINVAL(timeg(1 ,:)), MAXVAL(timeg(1 ,:))
      write (itime,101) '2 :', MINVAL(timeg(2 ,:)), MAXVAL(timeg(2 ,:))
      write (itime,101) '3 :', MINVAL(timeg(3 ,:)), MAXVAL(timeg(3 ,:))
      write (itime,101) '4 :', MINVAL(timeg(4 ,:)), MAXVAL(timeg(4 ,:))
      write (itime,101) '5 :', MINVAL(timeg(5 ,:)), MAXVAL(timeg(5 ,:))
      write (itime,101) '6 :', MINVAL(timeg(6 ,:)), MAXVAL(timeg(6 ,:))
      write (itime,101) '20:', MINVAL(timeg(20,:)), MAXVAL(timeg(20,:))
      write (itime,101) '21:', MINVAL(timeg(21,:)), MAXVAL(timeg(21,:))
      write (itime,101) '22:', MINVAL(timeg(22,:)), MAXVAL(timeg(22,:))
      write (itime,101) '23:', MINVAL(timeg(23,:)), MAXVAL(timeg(23,:))
      write (itime,101) '24:', MINVAL(timeg(24,:)), MAXVAL(timeg(24,:))
      write (itime,101) '25:', MINVAL(timeg(25,:)), MAXVAL(timeg(25,:))
      write (itime,101) '26:', MINVAL(timeg(26,:)), MAXVAL(timeg(26,:))
      write (itime,101) '27:', MINVAL(timeg(27,:)), MAXVAL(timeg(27,:))
      write (itime,101) '28:', MINVAL(timeg(28,:)), MAXVAL(timeg(28,:))
      write (itime,101) '29:', MINVAL(timeg(29,:)), MAXVAL(timeg(29,:))
      write (itime,101) '30:', MINVAL(timeg(30,:)), MAXVAL(timeg(30,:))
      write (itime,101) '40:', MINVAL(timeg(40,:)), MAXVAL(timeg(40,:))
      write (itime,101) '42:', MINVAL(timeg(42,:)), MAXVAL(timeg(42,:))
      write (itime,101) '50:', MINVAL(timeg(50,:)), MAXVAL(timeg(50,:))
      write (itime,101) '55:', MINVAL(timeg(55,:)), MAXVAL(timeg(55,:))
      write (itime,101) '60:', MINVAL(timeg(60,:)), MAXVAL(timeg(60,:))
      write (itime,101) '61:', MINVAL(timeg(61,:)), MAXVAL(timeg(61,:))
      write (itime,101) '62:', MINVAL(timeg(62,:)), MAXVAL(timeg(62,:))
      write (itime,101) '64:', MINVAL(timeg(64,:)), MAXVAL(timeg(64,:))
      write (itime,101) '66:', MINVAL(timeg(66,:)), MAXVAL(timeg(66,:))
      write (itime,101) '67:', MINVAL(timeg(67,:)), MAXVAL(timeg(67,:))
      write (itime,101) '68:', MINVAL(timeg(68,:)), MAXVAL(timeg(68,:))
      write (itime,101) '71:', MINVAL(timeg(71,:)), MAXVAL(timeg(71,:))
      write (itime,101) '72:', MINVAL(timeg(72,:)), MAXVAL(timeg(72,:))
      write (itime,101) '80:', MINVAL(timeg(80,:)), MAXVAL(timeg(80,:))
      write (itime,101) '81:', MINVAL(timeg(81,:)), MAXVAL(timeg(81,:))

      write (itime,100) 
      write (itime,100) 'Temps detail:'
      write (itime,100) 
      write (itime,100) 'Tcore-cpu: ',timeg(1 ,:)
      write (itime,100) 'Tcore-wc : ',timeg(2 ,:)
      write (itime,100) 'local    : ',timeg(3 ,:)
      write (itime,100) 'global   : ',timeg(4 ,:)
      write (itime,100) 'commu    : ',timeg(5 ,:)
      write (itime,100) 'e3       : ',timeg(6 ,:)
      write (itime,100) 'e3ivar   : ',timeg(20,:)
      write (itime,100) 'e3mtrx   : ',timeg(21,:)
      write (itime,100) 'e3conv   : ',timeg(22,:)
      write (itime,100) 'e3visc   : ',timeg(23,:)
      write (itime,100) 'e3LS     : ',timeg(24,:)
      write (itime,100) 'e3tau    : ',timeg(25,:)
      write (itime,100) 'e3eig1   : ',timeg(26,:)
      write (itime,100) 'e3eig2   : ',timeg(27,:)
      write (itime,100) 'e3juel   : ',timeg(28,:)
      write (itime,100) 'e3assm   : ',timeg(29,:)
      write (itime,100) 'e3bdg    : ',timeg(30,:)
      write (itime,100) 'e3b      : ',timeg(40,:)
      write (itime,100) 'e3nref   : ',timeg(42,:)
      write (itime,100) 'bc       : ',timeg(50,:)
      write (itime,100) 'vy       : ',timeg(55,:)
      write (itime,100) 'getdmc   : ',timeg(60,:)
      write (itime,100) 'sumgat-wa: ',timeg(61,:)
      write (itime,100) 'sumgat-wo: ',timeg(62,:)
      write (itime,100) 'gerbar   : ',timeg(64,:)
      write (itime,100) 'sponge   : ',timeg(66,:)
      write (itime,100) 'rstat-wai: ',timeg(67,:)
      write (itime,100) 'rstat-wor: ',timeg(68,:)
      write (itime,100) 'barrier 1: ',timeg(71,:)
      write (itime,100) 'barrier 2: ',timeg(72,:)
      write (itime,100) 'bl-elmmfg: ',timeg(80,:)
      write (itime,100) 'bl-itrres: ',timeg(81,:)

      close (itime)

      endif

100   format (1X,A,G14.6,F9.5)
101   format (1X,A,G14.6,G14.6)

      return
      end