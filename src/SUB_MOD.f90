MODULE SUB_MOD
USE Interface_MOD
USE gammaf90
USE gasdevf90
USE mGf90
IMPLICIT NONE

! declare all variables used within the model itself
public
integer           :: error   = 0
! number of samples to be output for the ensemble   
integer           :: EnsLen  = 2 

! Total number of iterations
integer           :: nruns   = 10

! number of runs for initial 'Burn-in'
integer           :: BurnInt = 100 
integer           :: jrun
integer           :: AMacc, DRacc, startrun = 0, length

CONTAINS
!-----------------------------------------------------------------------
function CalSSQE(Npars) result(SS)
implicit none
! Used to call the model but NOT write output (for most calls in the chain) 
real, intent(IN)  :: Npars(NPar)   ! Normalized parameters
real              :: Ymod(ANobs), OBS(ANobs)
real              :: Apars(NPar)  !Absolute parameters
real              :: SS(NDTYPE*Nstn)
  ! Convert normalized parameters to real parameters:
  Apars = Apv_(Npars)

  ! Calculate SSqE based on parameter input (a real model run):
  call CALCYSSE(Apars,Ymod,OBS,SS)
end function CalSSQE
!-----------------------------------------------------------------------
! Used to call the model with best parameters AND write output to the "bestout" file
subroutine model(subp, SS)
implicit none

! Current parameters
real,    intent(IN)  :: subp(NPar)

! Calculated SSQE
real,    intent(OUT) :: SS(NDTYPE * Nstn)

real    :: Ymod(ANobs), OBS(ANobs) 
integer :: don

! assign Apv at each step
Apv = Apv_(subp)
savefile = .TRUE.
call CalcYSSE(Apv,Ymod,OBS,SS)
savefile = .FALSE.

if (taskid .eq. 0) then
   ! Open bestout for writing
  open(bofint, file=bofn,status='replace',action='write')
  write(bofint, 3000) '  DOY    ',               &
                      'Depth    ',               &
                      'Data_type',               &
                      'Mod_Value', 'OBS_Value'
   ! Write the best simulation results into the file 
  do don = 1, ANobs
     write(bofint, 3100)  &
     OBS_DOY(don),OBS_Depth(don),OBS_Label(don),Ymod(don),OBS(don)
  enddo
  close(bofint)
endif
      
3000   format(1x,     3(A10),  3x, A10,3x, A10)
3100   format(1x, 2(F6.1, 2x), 5x, A5, 2(1x, 1pe12.3))
end subroutine model
!-----------------------------------------------------------------------
subroutine modelensout(efint, runnum, subp, SS)
  ! Used to call the model with the current parameter values 
  ! AND write output to the file for the Ensemble of model output 
  ! (only once every "outputint" runs)  
integer, intent(IN)  :: efint, runnum
real,    intent(IN)  :: subp(Np2Vary)
real,    intent(OUT) :: SS(NDTYPE  * Nstn)
real                 :: Ymod(ANobs), OBS(ANobs)
integer              :: don

if (taskid .eq. 0) then 
  Apv = Apv_(subp)

  ! Save model outputs
  savefile = .TRUE.
  call CalcYSSE(Apv,Ymod,OBS,SS)
  savefile = .FALSE.
  
   ! Open ensout for writing
  open(efint, file=eofn,status='old',&
       action='write',position='append')

  do don = 1, ANobs
     write(efint, 3200) runnum, OBS_DOY(don),OBS_Depth(don),&
                                OBS_Label(don), Ymod(don)
  enddo

  close(efint)
endif
3200   format(I6,1x, 2(F6.1, 2x), A10, 1x, 1pe12.4)
end subroutine modelensout
!-----------------------------------------------------------------------
subroutine modelnooutput(subp, SS)
implicit none
! Used to call the model but NOT write output (for most calls in the chain) 
real, intent(IN)  :: subp(Np2Vary)
real, intent(OUT) :: SS(NDTYPE * Nstn)
real              :: Ymod(ANobs), OBS(ANobs)
  ! The subp is the normalized value! 
  Apv      = Apv_(subp)
  savefile = .FALSE.
  call CalcYSSE(Apv,Ymod,OBS,SS)
end subroutine modelnooutput
!-----------------------------------------------------------------------
subroutine write_bestsigma
implicit none
integer :: i,k

if (taskid .eq. 0) then
   ! Write into best sigma file:
   open(bsfint, file=bsfn, action='write',status='replace' )
   write(bsfint,1200) 'At Run #        ',jrun
   write(bsfint,1210) 1D2*real(AMacc)/real(jrun-startrun),       &
                      1D2*real(DRacc)/real(jrun-startrun)
   write(bsfint,1220) 'Best LogL =    ', BestLogLike
   
   ! To calculate weighted SSqE(SS/sigma**2, Eq. 20), just to compare to previous non-adaptive MCMC runs
   TwtSSE = 0
   do k = 1, NDTYPE * Nstn
      TwtSSE = TwtSSE + BestSSqE(k)/(sigmabest(k)**2)
   enddo
   write(bsfint,1220) 'Best weighted SSqE = ', TwtSSE
   write(bsfint,1300)  NewLogLike, CurrLogLike, BestLogLike
   write(bsfint,1400) (SigmaLabel(i), sigmabest(i), sigmamean(i), i = 1, NDTYPE * Nstn)
   close(bsfint)
endif
1400 format('                  sigmabest              sigmamean ',/, &
           <NDTYPE*Nstn>(a15,2(1x,1pe20.3),/) )
1200 format(a15,1x,i16)
1210 format('** % 1st Accept. = ',1x,1f8.2,                         &
            '     ** % 2nd Accept. = ',1f8.2)

1220 format(a25,1x,1pe11.3)
1300 format('*** LogL:    New ',1pe11.3,'     Curr ',1pe11.3,  &
           '     Best ',1pe11.3)
1330 format('***        with  Subpcurr                Subpguess ',   &
           '    -->       Subpbest ','               Subpmean ',/,   &
           100(a15,4(1x,1pe20.3),/) )
end subroutine 
!-----------------------------------------------------------------------
subroutine write_bestpar
implicit none
integer :: i

Apvcurr = Apv_(subpcurr)
Apvguess= Apv_(subpguess)
Apvbest = Apv_(subpbest)
Apvmean = Apv_(subpmean)

if (taskid .eq. 0) then
   ! Write into best parameter file:
   open(bpfint, file=bpfn, action='write',status='replace'  )
   write(bpfint,1200) 'At Run #        ',jrun
   write(bpfint,1210) 1D2*real(AMacc)/real(jrun-startrun),         &
                      1D2*real(DRacc)/real(jrun-startrun)
   
   write(bpfint,1220) 'Best logL =    ', BestLogLike
   write(bpfint,1300)  NewLogLike, CurrLogLike, BestLogLike
   write(bpfint,1330) (ParamLabel(i),   &
     Apvcurr(i),                        &
     Apvguess(i),                       &
     Apvbest(i),                        &
     Apvmean(i), i = 1, Np2Vary)
   close(bpfint)
endif
1200 format(a15,1x,i16)
1210 format('** % 1st Accept. = ',1x,1f8.2,                         &
            '     ** % 2nd Accept. = ',1f8.2)
1220 format(a25,1x,1pe11.3)
1300 format('*** LogL:    New ',1pe11.3,'     Curr ',1pe11.3,  &
           '     Best ',1pe11.3)
1330 format('***        with  Subpcurr                Subpguess ',   &
           '    -->       Subpbest ','               Subpmean ',/,   &
           100(a15,4(1x,1pe20.3),/) )

end subroutine write_bestpar
!-----------------------------------------------------------------------
function Invmat(C) result(InvC)
implicit none
real, intent(in)  ::  C(NPar*(NPar+1)/2)
real              ::  InvC(NPar,NPar)

real              ::  work(NPar), U(NPar*(NPar+1)/2)
!!$ Invert the symmetric covariance matrix 
  call syminv(C, NPar, U, work, nullty, error )

!!$ Unpack the compressed inverse matrix U into the full (2-D) inverse covariance matrix
  call unpack(NPar,U,InvC)
end function
!-----------------------------------------------------------------------
function newsigma(SS) result(sigma_)
implicit none
integer             :: i,j,oi
real, intent(in)    :: SS(NDTYPE * Nstn)
real                :: sigma_(NDTYPE * Nstn)
!!$  sigma is randomly generated every time based on CurrSSqE
    ! Eq. 22 in Laine 2008
Do j = 1, Nstn
    do i = 1, NDTYPE
       oi = i+(j-1)*NDTYPE
       sigma_(oi) = sqrt( 1d0/                                          &
       gammadev( (n0d2(oi) +NDPTS(i,j)/2d0),                            &
                 1d0/(n0d2(oi)*S02(oi) + SS(oi)/2d0) ) )
    enddo
Enddo
end function newsigma
!-----------------------------------------------------------------------
logical function check_bounds(pars) result(NOTOK)
implicit none
real,    intent(in)  :: pars(NPar)
integer              :: k

NOTOK = .FALSE.
do k = 1, NPar  ! Loop over the Parameters to be varied

   if ((pars(k) .lt. MinValue(k)) .or. (pars(k).gt.MaxValue(k))) then
      NOTOK = .true.
      exit
   endif
        
enddo ! loop over Parameters 
return
end function
!-----------------------------------------------------------------------
function goodPar(Rchol,Pars) result(y)
implicit none
! Rchol is the lower triangluar  Cholesky factor of
! the Proposal Covariance Matrix, 
! which is calculated only every outputint steps 
! (it remains constant for each set of outputint steps). 
real,  intent(in)   :: Rchol(NPar*(NPar+1)/2)

! The desired mean of the new params
real,  intent(in)   :: Pars(NPar)
real  :: y(NPar)

! Absolute parameters
logical             :: ISBADPAR

! Pars is a one dimensional array holding proposed values for the subset of parameters being varied. 
    ISBADPAR = .true.
! If any parameter gets outside it's range, start over to generate a 
! different parameter set

    DO WHILE (ISBADPAR) 
    ! generate new set of parameters
       y = multiGauss(Rchol,Pars,NPar)
      
! Restrict all the parameters within the bounds: 
       cffpar  = Apv_(y)

       ISBADPAR = check_bounds(cffpar)
    ENDDO
End function
!-----------------------------------------------------------------------
function NewPAR(Pcvm_, OldPar) result(y)
implicit none
real, intent(in)   :: OldPar(NPar)
real, intent(in)   ::  Pcvm_(NPar*(NPar+1)/2)
real  :: y(NPar)

! Calculate the Cholesky factor for Pcvm, which is called Rchol here. 
call cholesky(Pcvm_,NPar,NPar*(NPar+1)/2,  &
                    Rchol,nullty,error)
y = goodPar(Rchol, OldPar)
end function
!-----------------------------------------------------------------------
subroutine write_ensemble
implicit none
integer i

cffpar = Apv_(subpcurr)
! Write Ensemble file(s) (Run #, LogLike and Parameters)
write(epfint,1850) jrun, CurrLogLike, (cffpar(i), i = 1, Np2Vary)
write(esfint,1850) jrun, CurrLogLike,      &
     (sigma(i),  i = 1, NDTYPE*Nstn),      &
     (CurrSSqE(i),i= 1, NDTYPE*Nstn)
1850 format(i9,1x,100(1pe12.3,2x))
end subroutine
!-----------------------------------------------------------------------
!!!
!!! second stage DR acceptance probability
!!! assumes Gaussian proposals, global variable iC = inv(C) is the inverse of
!!! first stage proposal of the point that was rejected
!!! 
subroutine MCMC_DR_alpha13(C, &
                           par1,      logLike1, &
                           par2,      logLike2, &
                           par3, ss3, logLike3, alpha13)

implicit none
real, intent(in) :: C(NPar*(NPar+1)/2)      ! the first stage proposal

! The SSqE of  new2 params
real, intent(in) ::  ss3(NDTYPE* Nstn)  

! Three sets of parameters
real, intent(in) :: par1(NPar), par2(NPar), par3(NPar)

! The probability of accepting the first move (already rejected)
real, intent(in) :: logLike1, logLike2

! The probability of the second move given the current position 
! and the first move (already rejected)
! logLike3 is the loglikelihood of the second move
real, intent(out):: alpha13, logLike3
real             :: alpha12
real             :: l2, q1, alpha32
real             :: iC(NPar, NPar)         ! the uncompacted form of inverse of C  
real             :: cff1(NPar), cff2(NPar),dpar32(NPar),dpar12(NPar)
!integer          :: m,q
! Calculate iC:
iC      = Invmat(C)

!  write(6,*) ' The Inverse of the Covariance Matrix: '
!
!  do m = 1, NPar
!     write(6,1100) (IC(m,q), q = 1, NPar)
!  end do
!
!1100 format(<NPar>(1pe9.1,1x))
! The probability of accepting the first move
alpha12 = min(1d0, exp(logLike2-logLike1))

! Calculate the logLike of the second move
logLike3 = CalcLogLike(ss3, sigma, par3)

! Calculate the probability (alpha32) of accepting y2 given y3
if (alpha12 == 0d0) then
   alpha32 = 0d0
else
! the loglikelihood from y3 to y2
!! oldpar is par3 and oldSS is ss3!!!
   alpha32 = min(1d0, exp(logLike2 - logLike3))
endif

! l2 = log(posterior(par3) - posterior(par1)  )
l2 = logLike3-logLike1

dpar32 = par3-par2
dpar12 = par1-par2

call matmuls(iC, dpar32, cff1) 
call matmuls(iC, dpar12, cff2) 

!  q1 = -0.5d0*( &
!       sum(  matmuls(iC,(par3-par2)) * (par3-par2) ) - &
!       sum(  matmuls(iC,(par1-par2)) * (par1-par2) ) )

 q1 = -0.5d0*(sum(cff1*dpar32) - sum(cff2 * dpar12))


  ! alpha2(x,y1,y2)
alpha13 = min(1d0, exp(l2+q1)*(1d0 - alpha32)/(1d0 - alpha12))
end subroutine MCMC_DR_alpha13
!-----------------------------------------------------------------------
subroutine MCMC_adapt
! Major MCMC subroutine
implicit none
integer:: Acc, Dcc
real   :: SS3(NDTYPE*Nstn), alpha13, DR_p
real   :: NewLogLike2
real   :: PCVM2(NPar*(NPar+1)/2)

! Data stored in each subprocess
real, allocatable :: par_(:), SSqE_(:), LogLike_(:), sigma_(:)

! Data stored in the root
real, allocatable :: par0(:), SSqE0(:), LogLike0(:), sigma0(:)
real, allocatable :: par1(:,:),SSqE1(:,:), sigma1(:,:)

! Output also includes:
! 1) Params, CVM, PCVM, SSqE, sigma, loglikehood at each step
! 2) Best of the above quantities so far
integer :: i,j, k, intv,row,col

integer :: stat(MPI_STATUS_SIZE)   ! required variable for receive routines

integer, parameter  :: tag2 = 2

! Downscaling factor for Delayed Rejection MCMC
real,    parameter  :: DRScale = 1d-2

! The number of times of adaptation (i.e. number of ensembles)
if ( (nruns-startrun) < EnsLen) then
   EnsLen  = 1
endif

intv = (nruns-startrun)/EnsLen

! At each process, construct a matrix to save the accepted params:
allocate(par_(intv * NPar)) ! Allocate the params for each process
allocate(LogLike_(intv))
allocate(sigma_(intv * NDTYPE * Nstn))
allocate( SSqE_(intv * NDTYPE * Nstn))

allocate(LogLike0(intv*numtasks))
LogLike0(:) = 0d0
   
if (taskid == 0) then
   allocate(par1(intv*numtasks, NPar))
   par1(:,:) = 0d0
   
   allocate(sigma1(intv*numtasks, NDTYPE*Nstn))
   sigma1(:,:) = 0d0
   
   allocate(SSqE1(intv*numtasks,  NDTYPE*Nstn))
   SSqE1(:,:) = 0d0
endif

! Initialize par_:
do i = 1, intv
  do j = 1, NPar
   par_((i-1)*NPar + j) = subpcurr(j)
  enddo
  
  LogLike_(i) = CurrLogLike

  do j = 1, NDTYPE * Nstn
     sigma_( (i-1)*NDTYPE*Nstn + j ) = sigma(j)
      SSqE_( (i-1)*NDTYPE*Nstn + j ) = CurrSSqE(j)
  enddo
enddo

MeanLogLike = 0d0
DO i = 1, EnsLen

 ! Synchronize all processes:
 if (MPIRUN==1) call MPI_BARRIER (MPI_COMM_WORLD,ierr)

 if (taskid .eq. 0) then
   if (i > 1) then
     ! write output to Ensemble files for simulated values
       call modelensout(eofint,jrun,subpcurr, dumE)
   endif
!   %%%%%%%end of block of outputting code
   do j = 1, numtasks-1
     ! Send the Pcvm (only) to child processes:
     call MPI_SEND(Pcvm, NPar*(NPar+1)/2, MPI_REAL8, j, tag2, MPI_COMM_WORLD, ierr)
   enddo
 else
   ! Receive Pcvm from root (child processes):
   call MPI_RECV(Pcvm, NPar*(NPar+1)/2, MPI_REAL8, 0, tag2, MPI_COMM_WORLD, stat, ierr)
 endif

 Acc = 0  !Acceptance counts for the 1st proposal
 Dcc = 0  !Acceptance counts for the 2nd proposal

 ! Start subcycle:
 subcycle: DO j = 1, intv
   
    call modelnooutput(subppro, SSqE)

    !   New log-likelihood
    NewLogLike = CalcLogLike(SSqE,     sigma, subppro)

    !   Old log-likelihood
    CurrLogLike = CalcLogLike(CurrSSqE, sigma,subpcurr)

    call random_number(DR_p)
    IF (log(DR_p) < (NewLogLike - CurrLogLike)) THEN
         CurrLogLike = NewLogLike
         CurrSSqE    = SSqE
         subpcurr    = subppro
         Acc         = Acc + 1
    ELSE
    ! DR:
    ! Propose a second move (Y2) based on scaled proposal covariance matrix
    ! and current position
         Pcvm2      = DRScale*Pcvm
         subppro2   = NewPAR(Pcvm2, subpcurr)
       !      Calcuate the new SSqE
         SS3        = CalSSQE(subppro2)

       ! Calculate the acceptance probability of the second move
         call MCMC_DR_alpha13(Pcvm,      &
                subpcurr,  CurrLogLike,  &
                subppro,    NewLogLike,  &
                subppro2,      SS3,  NewLogLike2, alpha13)

         ! Judge whether the second move should be accepted or not
         call random_number(DR_p)
         if (DR_p < alpha13) then
             !accept the second move (Y1)
            CurrLogLike = NewLogLike2
            CurrSSqE    = SS3
            subpcurr    = subppro2
            Dcc         = Dcc + 1
         endif
    ENDIF

    ! save to the par_ (relative params!!)
     do col = 1, NPar
        par_((j-1)*NPar + col) = subpcurr(col)
     enddo

    ! save to LogLike_
    LogLike_(j) = CurrLogLike  
    
    ! save to SSqE and sigma:
    do col = 1, NDTYPE * Nstn
       SSqE_((j-1)*NDTYPE*Nstn + col) = CurrSSqE(col)
      sigma_((j-1)*NDTYPE*Nstn + col) =    sigma(col)
    enddo

    ! Propose newpar based on old par
    subppro = NewPAR(Pcvm, subpcurr)

    ! Update sigma:
    sigma = newsigma(CurrSSqE)

   ENDDO subcycle ! END of loop of each cycle 
   
   write(6,400) taskid, dble(Acc)/dble(intv)
   write(6,401) taskid, dble(Dcc)/dble(intv)

400 format('For taskid', 1x, I0, ', Acceptance rate for the 1st proposal is', F12.3)
401 format('For taskid', 1x, I0, ', Acceptance rate for the 2nd proposal is', F12.3)

   ! Synchronize all processes:
   call MPI_BARRIER (MPI_COMM_WORLD,ierr)

   ! MPI root: receive matrix of accepted params from child processes:
   ! Collect the data (par_, LogLike_, sigma, SSqE_, acceptance ratios) from all tasks to the root:
   ! The dimension of output data cannot exceed 30*50:
   if (.not. ALLOCATED(par0))   allocate(par0(intv*NPar*numtasks))
   par0(:)   = 0d0
   if (.not. ALLOCATED(sigma0)) allocate(sigma0(intv*NDTYPE*Nstn*numtasks))
   sigma0(:) = 0d0
   if (.not. ALLOCATED(SSqE0))  allocate(SSqE0(intv*NDTYPE*Nstn*numtasks))
   SSqE0(:)  = 0d0

   call MPI_Gather(par_, NPar*intv,         MPI_REAL8,     &
                   par0, NPar*intv,         MPI_REAL8,     &
                    0,  MPI_COMM_WORLD, ierr)

   call MPI_Gather(LogLike_, intv,          MPI_REAL8,     &
                   LogLike0, intv,          MPI_REAL8,     &
                    0,  MPI_COMM_WORLD, ierr)

   call MPI_Gather(sigma_, NDTYPE*intv*Nstn, MPI_REAL8,     &
                   sigma0, NDTYPE*intv*NStn, MPI_REAL8,     &
                    0,  MPI_COMM_WORLD, ierr)

   call MPI_Gather(SSqE_,  NDTYPE*intv*Nstn, MPI_REAL8,     &
                   SSqE0,  NDTYPE*intv*NStn, MPI_REAL8,     &
                    0,  MPI_COMM_WORLD, ierr)

   IF (taskid == 0) THEN
      DO k = 1, intv*numtasks
         ! Convert par0 into a matrix (par1):
         do col = 1, NPar
            par1(k, col) = par0((k-1)*NPar + col) 
         enddo

         ! Convert sigma0 and SSqE0 into a matrix:
         do col = 1, Nstn * NDTYPE
            sigma1(k, col) = sigma0((k-1)*Nstn*NDTYPE + col) 
             SSqE1(k, col) =  SSqE0((k-1)*Nstn*NDTYPE + col) 
         enddo
      ENDDO

      ! Open enspar for writing
      open(epfint,file=epfn,status='old',&
           action='write',position='append')

      ! Open enssig for writing
      open(esfint,file=esfn,status='old',&
           action='write',position='append')

      do k = 1, intv*numtasks
         cffpar = Apv_(par1(k,:))
         ! Write Ensemble file(s) (Run #, LogLike and Parameters)
         write(epfint,1850) LogLike0(k), (cffpar(col), col = 1,NPar)
         write(esfint,1850) LogLike0(k),                   &
              (sigma1(k,col),  col = 1, NDTYPE*Nstn),      &
              ( SSqE1(k,col),  col = 1, NDTYPE*Nstn)
      enddo

      close(epfint)
      close(esfint)

      ! MPI root: find the best LogL
      col  = maxloc(LogLike0, 1)
      DR_p = maxval(LogLike0)

      !	this is to keep track of the 'best' run which 
      !	is not strictly speaking the purpose of the
      !	assimilation but an interesting output
      if(DR_p > BestLogLike ) then
         BestLogLike = DR_p     
         BestSSqE    =  SSqE1(col,:)
         sigmabest   = sigma1(col,:)
         subpbest    =   par1(col,:)  
      endif

! MPI root: Calculate new CVM based on collected accepted params:
      call cov(intv*numtasks, par1, Cvm)  

!   Calculate new PCVM based on CVM
!   Set the Proposal covariance matrix, by scaling the Parameter Covariance matrix
      Pcvm = Cvm*Spcvm/NPar

!   Add a small term to the diagonal entries, so that the matrix will not be singular. 
      do k = 1, NPar
        Pcvm(k*(k+1)/2)=Pcvm(k*(k+1)/2) + CvEpsilon*Spcvm/NPar
      enddo

      write(6, *) 'Proposal Covariance matrix for ',i+1,' ensemble run ='
      DO row = 1, NPar
         write(6,3000) (Pcvm(row*(row-1)/2+col), col = 1, row )
      End do

      ! Release memory
      if (allocated(par0))   deallocate(par0)
      if (allocated(sigma0)) deallocate(sigma0)
      if (allocated(SSqE0))  deallocate(SSqE0)

      ! Update jrun:
      jrun = jrun + intv * numtasks

     endif  !==> end of taskid ==0
ENDDO
1850 format(100(1pe12.3,2x))
3000 format(5x,<NPar>(1pe8.1,1x))
End subroutine MCMC_adapt
!-----------------------------------------------------------------------
END MODULE sub_mod
