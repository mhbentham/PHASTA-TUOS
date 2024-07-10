        subroutine get_bcredist(x,  y, iBCredist, BCredist,
     &                          primvert, primvertval) 
c
c----------------------------------------------------------------------
c
c 
c----------------------------------------------------------------------
c
        use pvsQbi  ! brings in NABI
        use stats   !  
        use pointer_data  ! brings in the pointers for the blocked arrays
        use local_mass
        use redist_freeze
c
        include "common.h"
        include "mpif.h"

c
        dimension x(numnp,nsd)               
        dimension y(nshg,ndof)
        real*8  BCredist(nshg), primvertval(nshg)
        integer iBCredist(nshg), primvert(nshg)
c
c JMR Debug
c
c      do iedge = 1, nedges_tet
c        ivtx1 = edge_vtx_tet(iedge,1)
c        ivtx2 = edge_vtx_tet(iedge,2)
c        write(*,1000) iedge, ivtx1, ivtx2
c      enddo
c 1000   format (" edge: ",i6," has vertices ",i4," and ",i4)
c
c
c initialize array
c
        do inode = 1, nshg
          BCredist(inode) = zero
          iBCredist(inode) = 0	
          primvert(inode) = 0
          primvertval(inode) = zero
        enddo
c
c loop over interior element blocks
c
        
           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lelCat = lcblk(2,iblk)
              lcsyst = lcblk(3,iblk)
              iorder = lcblk(4,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              nsymdl = lcblk(9,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              ngauss = nint(lcsyst)

c
c compute levelset for Scalar 2 for primary nodes.
c Primary nodes are the nodes of elements of which the interace intersects.
c
               call e3redistfreeze (x,y,mien(iblk)%p,
     &                              iBCredist,BCredist,
     &                              primvert, primvertval)

           enddo
c
        return
        end






c..***********************************************************
c....This routine identifies those elements that are interected
c    by the interface and then computes the level set for the
c    element nodes.  The BCredist and iBCredist arrays are
c    set to enforce those computed scalar values.
c
c    Note that iBCredist and BCredist compliment the BCs set
c    in the standard IBC and BC arrays
c*************************************************************
        subroutine e3redistfreeze (x,y,ien,iBCredist,BCredist,
     &             primvert, primvertval)
c                                                                      
c----------------------------------------------------------------------
c
c input: 
c  x     (numnp,nsd)		: nodal coordinates at current step
c  y      (nshg,ndof)		: nodal solution
c  ien (npro,nshl)
c output:
c  iBCredist(nshg)		: BC code for Scalar 2(0=no BC st, 1=BC set) 
c  BCredist(nshg)		: BC value for scalar 2
c  primvert(nshg)		: flag if primary vertex
c  primvertval(nshg)		: primary vertex value
c  
c
c----------------------------------------------------------------------
c
        include "common.h"
c....Passed arrays
        dimension shp(nshl,ngauss),     shgl(nsd,nshl,ngauss),
     &            x(numnp,nsd),         y(nshg,ndof),
     &            ien(npro,nshl)
        integer  iBCredist(nshg), primvert(nshg)
        real*8   BCredist(nshg), primvertval(nshg)

c
c local arrays
c
        real*8 dxidx(npro,nsd,nsd),WdetJ(npro)

        real*8 tmp(npro), disttmp(npro)
        real*8 dxdxi(npro,nsd,nsd)

        real*8 sgn(npro,nshl),  shape(npro,nshl),
     &         shdrv(npro,nsd,nshl)
        real*8 xl(npro,nenl,nsd), yl(npro,nshl,ndofl),
     &         BCredistl(npro,nenl), primvertvall(npro,nenl)
        integer iBCredistl(npro,nenl), primvertl(npro,nenl)
        real*8 scalar(npro,nshl), xpts(4,nsd), normvec(nsd)

        integer ipro, numpts

c
c Set sign
c
        do i=1,nshl
          where ( ien(:,i) < 0 )
            sgn(:,i) = -one
          elsewhere
            sgn(:,i) = one
          endwhere
        enddo
c
c*************Localizing coordinates and solution*************
c	
        call localx (x,xl,ien,nsd,'gather  ')
        call localy (y,yl,ien,ndofl,'gather  ')
        call localBC (iBCredist,BCredist,iBCredistl,BCredistl,
     &                ien,1,'gather  ')
        call localBC (primvert, primvertval,primvertl,primvertvall,
     &                ien,1,'gather  ')

        ib = 5+isclr
        scalar(:,:) = yl(:,:,ib)
c
c Loop over elements
c
        do ipro = 1, npro
c          write(*,*) "local elem ", ipro," node 1 ",xl(ipro,1,:)
c          write(*,*) "local elem ", ipro," node 2 ",xl(ipro,2,:)
c          write(*,*) "local elem ", ipro," node 3 ",xl(ipro,3,:)
c          write(*,*) "local elem ", ipro," node 4 ",xl(ipro,4,:)
c          write(*,*) "local elem ", ipro," scalar: ",scalar(ipro,:)
          i_primary_flag = 0
          if((maxval(scalar(ipro,:)).ge.0.0) .and. 
     &       (minval(scalar(ipro,:)).le.0.0)) then    !true = primary element
c             write(*,*) "is primary : local elem ", ipro
             i_primary_flag = 1
c             write(*,*) "Element is primary: ", ipro
             primvertl(ipro,:) = 1
             call calc_interface_pts(xl(ipro,:,:), scalar(ipro,:), xpts, 
     &                               numpts)
             call calc_normal(xpts, numpts, normvec)
             call compute_dist(xpts, normvec, xl(ipro,:,:), 
     &                         scalar(ipro,:), numpts, 
     &                         iBCredistl(ipro,:),
     &                         BCredistl(ipro,:))

c             write(*,*) "  numpts = ", numpts
c             write(*,*) "  xpts 1 = ", xpts(1,:)
c             write(*,*) "  xpts 2 = ", xpts(2,:)
c             write(*,*) "  xpts 3 = ", xpts(3,:)
c             write(*,*) "  xpts 4 = ", xpts(4,:)
c             write(*,*) "  iBCredistl = ", iBCredistl(ipro,:)
c             write(*,*) "  BCredistl = ", BCredistl(ipro,:)
             where (iBCredistl(ipro,:).eq.1)
               primvertl(ipro,:) = 2
               primvertvall(ipro,:) =  BCredistl(ipro,:)
             endwhere

          endif
        enddo

cc JMR Debug
cc  write out primvertl and primvertvall arrays
c        do  ipro = 1, npro
c          write(*,2000) ipro, primvertl(ipro,:), primvertvall(ipro,:)
c        enddo
c 2000 format ("local elem ",i6," primvert: ",4(1x,i2),
c     &" primvertvall: ",4(1x,e10.4))
c        do  ipro = 1, npro
c          write(*,2002) ipro, iBCredistl(ipro,:), BCredistl(ipro,:)
c        enddo
c 2002 format ("local elem ",i6," iBCredistl: ",4(1x,i2),
c     &" BCredistl: ",4(1x,e10.4))


c scatter back to global level
        call localBC (iBCredist, BCredist, iBCredistl, BCredistl,
     &                ien, 1, 'scatter ')
        call localBC (primvert, primvertval,primvertl,primvertvall,
     &                ien, 1, 'scatter ')
c        do inode = 1, nshg
c           write(*,2004) inode, iBCredist(inode), BCredist(inode)
c        enddo
c 2004 format ("node ",i6," iBCredist: ",i6," BCredist: ",e10.4)
cc
c        do inode = 1, nshg
c           write(*,2006) inode, primvert(inode), primvertval(inode)
c        enddo
c 2006 format ("node ",i6," primvert: ",i6," primvertval: ",e10.4)
c
c
       return
       end


      subroutine calc_interface_pts(xl, scalar, xpts, numpts)
c
      use redist_freeze !access to BC arrays if freezing value of primary LS vertices
      include "common.h"
c
      dimension xl(nshl,nsd),   scalar(nshl),
     &          xpts(4,nsd)

      real*8 rij(nsd)
      real*8 factor

      real*8 xptstmp(4,nsd)
      integer ivtx1, ivtx2, iord(4), int_order(4)

c initialize
      numpts = 0
      do iedge = 1, nedges_tet
        ivtx1 = edge_vtx_tet(iedge,1)
        ivtx2 = edge_vtx_tet(iedge,2)
        if ((scalar(ivtx1) * scalar(ivtx2)).lt.0.0) then  !interface crosses this edge
          numpts = numpts+1
          rij(:) = xl(ivtx2,:) - xl(ivtx1,:)
          factor = (scalar(ivtx1)-0)/(scalar(ivtx1)-scalar(ivtx2))
          xpts(numpts,:) = xl(ivtx1,:) + factor*rij(:)
          iord(numpts) = iedge
        endif
      enddo
      if (numpts.gt.4) then
         write(*,*) "ERROR: numpts = ",numpts," > 4"
      endif
c correct the order for 4-pt intersections
      if (numpts .eq. 4) then
c        int_order(:) = 
c   &       int_order_tet_4(iord(1),iord(2),iord(3),iord(4),:) 
        int_order(1) = 
     &         int_order_tet_4(iord(1),iord(2),iord(3),iord(4),1)
        int_order(2) =
     &         int_order_tet_4(iord(1),iord(2),iord(3),iord(4),2)
        int_order(3) =
     &         int_order_tet_4(iord(1),iord(2),iord(3),iord(4),3)
        int_order(4) =
     &         int_order_tet_4(iord(1),iord(2),iord(3),iord(4),4)
        do index=1,4
          xptstmp(index,:) = xpts(int_order(index),:)
        enddo
        xpts(:,:) = xptstmp(:,:)
      endif
c
      return
      end


c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine calc_normal(xpts, numpts, normvec)
c
      include "common.h"
c
      real*8  normvec(nsd),   xpts(4,nsd)

      real*8 v1(nsd), v2(nsd)
      real*8 temp
c
      v1(:) = xpts(2,:) - xpts(1,:)
      v2(:) = xpts(3,:) - xpts(1,:)
c 
      call cross_product(v1, v2, normvec)
      temp = one / 
     &       sqrt (normvec(1)**2 + normvec(2)**2 + normvec(3)**2)
      normvec(:) = normvec(:) * temp
c
      return 
      end


      subroutine cross_product(v1, v2, v3)
c
c This routine computes the cross product v1 x v2
c and returns the result in v3
c
      include "common.h"
c
      real*8 v1(nsd), v2(nsd), v3(nsd)
c
      v3(1) = v1(2) * v2(3) - v2(2) * v1(3)
      v3(2) = v2(1) * v1(3) - v1(1) * v2(3)
      v3(3) = v1(1) * v2(2) - v2(1) * v1(2)
c
      return
      end


      subroutine dot_product(v1, v2, v3)
c
c This routine computes the dot product: v1 dot v2
c and returns the result in v3
c
      include "common.h"
c
      real*8 v1(nsd), v2(nsd)
      real*8 v3
c
      v3 = v1(1)*v2(1) +
     &     v1(2)*v2(2) +
     &     v1(3)*v2(3)
c
      return
      end


      subroutine compute_dist(xpts, normvec, xl, scalar,
     &                        numpts, iBCredistl, BCredistl)
c
c  This routine computes the normal ditance from a vertice to the intersection plane.
c  If the intersection point lies within the element then the scalar 2 value is updated
c  and enforced at that value during successive redistancing steps.
c
      include "common.h"
c

      real*8 xpts(4,nsd),         xl(nshl,nsd),  
     &        scalar(nshl),       BCredistl(nshl) ,
     &        normvec(nsd)
      integer iBCredistl(nshl)
c
      real*8 vec(nsd), xintpt(nsd)
      real*8 dist
      integer iflag_inplane
c
      do inode=1,nshl
        vec(:) = xl(inode,:) - xpts(1,:)
        call dot_product(vec,normvec,dist)
        sign = scalar(inode) / abs(scalar(inode))
        xintpt(:) = xl(inode,:) - dist*normvec(:)
c test if intersection point is within intersection plane boudned by element
        iflag_inplane = 0
        call is_in_intplane(xintpt,xpts,numpts,iflag_inplane)
        if (iflag_inplane .ne. 0) then
          dist = abs(dist) * sign
          if ((iBCredistl(inode).eq.1) .and. 
     &        (abs(dist).lt.abs(BCredistl(inode)))) then
            BCredistl(inode) = dist
c            write(*,*) "scalar kept at phi = ", scalar(inode)
          else
            BCredistl(inode) = dist
c            write(*,*) "scalar changed from ",scalar(inode),"to ",dist
          endif
          iBCredistl(inode) = 1
        endif
      enddo
c
      return
      end
               


      subroutine is_in_intplane(xintpt,xpts,numpts,iflag)
c
c This routine determines if the point, xintpt, is within 
c the portion of the plane encompassed by the points, xpts.
c It is assumed that xintpt is within the same plane defined
c by xpts.
c
      include "common.h"

      dimension xintpt(nsd)
      dimension xpts(4,nsd)
      integer numpts, iflag
c
      integer ivec4(4,2), ivec3(3,2)
      data ((ivec4(i,j),i=1,4),j=1,2) / 2, 3, 4, 1,   3, 4, 1, 2 /
      data ((ivec3(i,j),i=1,3),j=1,2) / 2, 3, 1,      3, 1, 2 /
      real*8 vtx1(nsd), vtx2(nsd), vtx_offline(nsd)
      integer isameside
c
c
c Test that xintpt is in same plane as xpts - fail if not
c (must have same normal)
c
c  WILL DO THIS LATER

c Test if point is in plane by testing that the point
c is on the "inside" side of each edge.  The "inside" 
c side of the edge is the same side of the edge as the
c remaining points of the area.  This is doen by computing 
c the normal for 
c    1. the edge vector and the vector from an edge vertex 
c       to the point.
c    2. the edge vector and the vector from an edge vertex
c       and one of the other plane points in xpts.
c If the normal point in the same direction then the point,
c xintpt, lies on the same side of the edge as the other
c point in xpts.  If this happens for all edges, then xintpt
c must be contained within the area defined by xpts.
c
c Note that this method only works if the points in xpts
c are in an order that walks round the area.  This was done 
c in calc_interface_pts(). 
c
c So loop over the edges in xpts
c
      iflag = 1
      if (numpts .eq. 3) then
        do ipt = 1, numpts
          vtx1 = xpts(ipt,:)
          vtx2 = xpts(ivec3(ipt,1),:)
          vtx_offline = xpts(ivec3(ipt,2),:) 
          call onsameside(xintpt, vtx_offline, vtx1, vtx2, isameside)
          if (isameside .eq. 0) then
            iflag = 0
            exit 
          endif
        enddo
      else if (numpts .eq. 4) then
        do ipt = 1, numpts
          vtx1 = xpts(ipt,:)
          vtx2 = xpts(ivec4(ipt,1),:)
          vtx_offline = xpts(ivec4(ipt,2),:) 
          call onsameside(xintpt, vtx_offline, vtx1, vtx2, isameside)
          if (isameside .eq. 0) then
            iflag = 0
            exit
          endif
        enddo
      else
        write(*,*) "***ERROR*** - numpts not equal to 3 or 4"
      endif
c
      return
      end

      

      subroutine onsameside (p1, p2, vtxa, vtxb, iflag)
c
c function to test if point p1 is on the same side of line a-b
c as point p2.
c Returns iflag = 1 if true, iflag=0 otherwise.
c
      include "common.h"

      real*8 p1(nsd), p2(nsd), vtxa(nsd), vtxb(nsd)
      real*8 cp1(nsd), cp2(nsd)
      real*8 dp
c
      call cross_product (vtxb-vtxa, p1-vtxa, cp1)
      call cross_product (vtxb-vtxa, p2-vtxa, cp2)
      call dot_product (cp1, cp2, dp)
      if (dp .ge. 0) then
        iflag = 1
      else
        iflag = 0
      endif
 
      return
      end

