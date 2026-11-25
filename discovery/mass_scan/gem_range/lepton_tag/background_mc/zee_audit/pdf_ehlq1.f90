! Bind the proton densities implemented in PYTHIA 5.6 to the event audit.
subroutine ehlq1(n, ids, xs, qs, xfs) bind(C)
  use iso_c_binding
  implicit none
  integer(c_int), value :: n
  integer(c_int), intent(in) :: ids(n)
  real(c_double), intent(in) :: xs(n), qs(n)
  real(c_double), intent(out) :: xfs(n)
  integer :: mstp(200), msti(200), i
  real :: parp(200), pari(200), proton(-25:25)
  common /pypars/ mstp, parp, msti, pari
  mstp(51) = 1
  mstp(52) = 1
  do i = 1, n
     call pystfu(2212, real(xs(i)), real(qs(i)**2), proton)
     xfs(i) = dble(proton(ids(i)))
  end do
end subroutine
