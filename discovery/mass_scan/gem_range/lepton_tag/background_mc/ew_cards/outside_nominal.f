      logical function dummy_cuts(p)
c     Complement of the three nominal EW run-card acceptances, evaluated
c     with the same MG5 pT, rapidity and DeltaR functions as cuts.f.
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      double precision p(0:3,nexternal),pt,rap,r2
      double precision pair(0:3),m2,pmax
      external pt,rap,r2
      logical is_a_j(nexternal),is_a_a(nexternal),is_a_l(nexternal)
      logical is_a_b(nexternal),is_a_nu(nexternal),is_heavy(nexternal)
      logical is_a_onium(nexternal),do_cuts(nexternal),inside
      common/to_specisa/is_a_j,is_a_a,is_a_l,is_a_b,is_a_nu,is_heavy,
     &                  is_a_onium,do_cuts
      integer i,k,ngam,nlep,igam(2),ilep(2)
      ngam=0
      nlep=0
      inside=.true.
      do i=nincoming+1,nexternal
         if(is_a_a(i))then
            ngam=ngam+1
            if(ngam.gt.2)stop 'unexpected photon multiplicity'
            igam(ngam)=i
            if(pt(p(0,i)).lt.14d0.or.abs(rap(p(0,i))).gt.3d0)
     &           inside=.false.
         else if(is_a_l(i))then
            nlep=nlep+1
            if(nlep.gt.2)stop 'unexpected lepton multiplicity'
            ilep(nlep)=i
         endif
      enddo
      if(ngam.eq.2)then
         do k=0,3
            pair(k)=p(k,igam(1))+p(k,igam(2))
         enddo
         m2=pair(0)**2-pair(1)**2-pair(2)**2-pair(3)**2
         if(m2.lt.60d0**2.or.m2.gt.170d0**2) inside=.false.
         pmax=0d0
         do i=1,nlep
            pmax=max(pmax,pt(p(0,ilep(i))))
         enddo
         if(pmax.lt.13d0) inside=.false.
         if(nlep.eq.1)then
            if(abs(rap(p(0,ilep(1)))).gt.3.5d0)inside=.false.
         endif
      else if(ngam.eq.1.and.nlep.eq.2)then
         do k=1,nlep
            i=ilep(k)
            if(pt(p(0,i)).lt.14d0.or.abs(rap(p(0,i))).gt.3d0
     &          .or.r2(p(0,i),p(0,igam(1))).lt.0.2d0**2)inside=.false.
         enddo
      else
         stop 'unknown EW control channel'
      endif
      dummy_cuts=.not.inside
      return
      end
