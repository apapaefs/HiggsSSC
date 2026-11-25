! PYTHIA 5.6 / JETSET 7.3. The native subprocess-19 cross section already
! includes the branching fractions for the enabled Z decays.
program vintage
  implicit none
  integer :: msel,msub(200),kfin(2,-40:40),mstp(200),msti(200)
  integer :: mdcy(500,3),mdme(2000,2),kfdp(2000,5),mrlu(6),mstu(200),mstj(200)
  real :: ckin(200),parp(200),pari(200),brat(2000),rrlu(100),paru(200),parj(200)
  common/pysubs/msel,msub,kfin,ckin
  common/pypars/mstp,parp,msti,pari
  common/ludat1/mstu,paru,mstj,parj
  common/ludat3/mdcy,mdme,brat,kfdp
  common/ludatr/mrlu,rrlu
  integer :: seed,n,flavours,i,ev
  character(32) :: argument
  call get_command_argument(1,argument)
  read(argument,*) seed
  call get_command_argument(2,argument)
  read(argument,*) n
  call get_command_argument(3,argument)
  read(argument,*) flavours
  if (seed<=0.or.n<=0.or.(flavours/=1.and.flavours/=2)) stop 1
  mrlu(1)=seed
  msel=0
  msub(19)=1
  ckin(3)=10.
  mstp(51)=1
  mstj(41)=2
  if (command_argument_count()>3) then
    call get_command_argument(4,argument)
    select case (trim(argument))
    case ('nominal')
      continue
    case ('no_isr')
      mstp(61)=0
    case ('no_timelike_qed')
      mstj(41)=1
    case ('no_mpi')
      mstp(81)=0
    case ('isr_half_scale')
      parp(64)=0.25
    case ('isr_double_scale')
      parp(64)=4.
    case default
      error stop 'Unknown PYTHIA radiation or MPI setting'
    end select
    write(*,'(A,3I4)') 'SWITCHES ',mstp(61),mstj(41),mstp(81)
    write(*,'(A,I4,F8.2)') 'ISR_SCALE ',mstp(64),parp(64)
  end if
  do i=mdcy(23,2),mdcy(23,2)+mdcy(23,3)-1
    mdme(i,1)=0
    if(abs(kfdp(i,1))==11.and.abs(kfdp(i,2))==11) mdme(i,1)=1
    if(flavours==2.and.abs(kfdp(i,1))==13.and.abs(kfdp(i,2))==13) mdme(i,1)=1
  end do
  call pyinit('CMS','p','p',40000.)
  do ev=1,n
    call pyevnt()
    if (flavours==1) call measure(seed,ev)
  end do
  write(*,'(A,3I12,ES20.10)') 'RATE ',seed,n,flavours,dble(pari(1))*1.d9
  if (flavours==1) call finish()
end program
