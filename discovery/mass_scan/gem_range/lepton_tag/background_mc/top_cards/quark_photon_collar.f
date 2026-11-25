      logical function dummy_cuts(p)
c     Require at least one final photon within 0.05 of a final jet
c     parton. The run card supplies the lower separation of 0.02.
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      integer i,j
      double precision p(0:3,nexternal),r2
      external r2
      logical is_a_j(nexternal),is_a_a(nexternal),is_a_l(nexternal)
      logical is_a_b(nexternal),is_a_nu(nexternal),is_heavy(nexternal)
      logical is_a_onium(nexternal),do_cuts(nexternal)
      common/to_specisa/is_a_j,is_a_a,is_a_l,is_a_b,is_a_nu,is_heavy,
     &                  is_a_onium,do_cuts
      dummy_cuts=.false.
      do i=nincoming+1,nexternal
         if(is_a_a(i))then
            do j=nincoming+1,nexternal
               if(is_a_j(j).or.is_a_b(j))then
                  if(r2(p(0,i),p(0,j)).lt.0.05d0**2)then
                     dummy_cuts=.true.
                     return
                  endif
               endif
            enddo
         endif
      enddo
      return
      end
