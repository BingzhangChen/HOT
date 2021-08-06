MODULE PARAM_MOD
USE MPI
USE BIOCOM_MOD
IMPLICIT NONE

! Number of stations to run
integer, parameter      :: Nstn  = 1

! Station name:
character(4), parameter :: Stn(Nstn) = ['HOT']

integer, parameter :: bot_bound  = Neumann ! Option of bottom boundary condition

! Current biological MODEL!
integer, parameter :: Model_ID   = GeiderDroop

! Number of phytoplankton groups
integer, parameter :: NPHY       = 1
logical, parameter :: DO_IRON    = .FALSE.  ! Whether involve iron or not
integer, parameter :: N_fer      = 33    ! Number of vertical layers for Iron

!! End of definding external forcing variables

! Indices for state variables
integer, parameter :: iNO3   = 1
integer, private   :: i
integer, parameter :: iPHY(NPHY) = (/ (i+iNO3, i=1,NPHY) /)
integer, parameter :: iZOO   = iPHY(NPHY) + 1
integer, parameter :: iDET   = iZOO+1
integer, parameter :: iPHYC(NPHY) = (/ (i+iDET,        i=1,NPHY) /) 
integer, parameter ::  iCHL(NPHY) = (/ (i+iPHYC(NPHY), i=1,NPHY) /) 
integer, parameter :: NVAR   = iCHL(NPHY) ! Total number of biological tracers

! Useless, but needed for compilation
integer            :: iVNO3, iVPHY,iCOVNP, iVZOO, iCOVPZ, iCOVNZ  

! Define tracer matrix:
real               :: Vars(NVAR,nlev) = 0d0

! Define the number of sinking tracers:
integer, parameter :: NVsinkterms =  1  ! DET

! Define the index of sinking tracers in the Vars matrix: 
integer, parameter :: Windex(NVsinkterms) = [iDET]

! Indices for additional output variables on top of state variables
integer, parameter :: oCHLt       = NVAR       + 1
integer, parameter :: oLno3(NPHY) = (/ (i+oCHLt,      i=1,NPHY) /)
integer, parameter :: oSI(NPHY)   = (/ (i+oLno3(NPHY),i=1,NPHY) /)
integer, parameter :: oQN(NPHY)   = (/ (i+oSI(NPHY),  i=1,NPHY) /)
integer, parameter :: othe(NPHY)  = (/ (i+oQN(NPHY),  i=1,NPHY) /)
integer, parameter :: oNPP        = othe(NPHY) + 1
integer, parameter :: oPON        = oNPP       + 1
integer, parameter :: omuNet(NPHY)= (/ (i+oPON,        i=1,NPHY) /)
integer, parameter :: oGraz(NPHY) = (/ (i+omuNet(NPHY),i=1,NPHY) /)

! The diffusion output order must be consistent with the tracer order!
integer, parameter :: oD_VARS(NVAR) = [(oGraz(NPHY) + i, i = 1, NVAR)]
integer, parameter :: Nout   = oD_VARS(NVAR)

! Initialize Varout matrix:
real               :: Varout(Nout, nlev) = 0d0

! Initialize labels for Varout:
character(LEN=10)  :: Labelout(Nout+ow)  = 'Unknown'

! Define parameter indices:
integer, parameter :: imu00            = 1
integer, parameter :: iaI0             = imu00 +1
integer, parameter :: iKN              = iaI0  +1
integer, parameter :: iKp              = iKN   +1
integer, parameter :: imz              = iKp   +1
integer, parameter :: iwDET            = imz   +1  ! Sinking rate of detritus
integer, parameter :: NPar             = iwDET     ! Total number of parameters
real               :: params(NPar)     = 0d0  ! Define parameters
character(LEN=8)   :: ParamLabel(NPar) = 'Unknown' !Define parameter labels

!  Maxima and minima  of parameter values 
real               :: MaxValue(NPar)   = 0d0, MinValue(NPar) = 0d0
real,    parameter :: Topt(NPHY)       = [(-3. + dble(i), i = 1, NPHY)]
real,    parameter :: alphaG           = 1.1  !Kill-the-winner coefficient

integer, parameter :: grazing_formulation = 3

CONTAINS

SUBROUTINE choose_model
IMPLICIT NONE

character(len=10), parameter :: format_string = "(A3,I0)"
IF(TASKID==0) WRITE(6,*) 'Geider-Droop phytoplankton model selected!'

! Assign label names for Varout:
Labelout(oTemp     ) = 'Temp'
Labelout(oPAR      ) = 'PAR '
Labelout(oAks      ) = 'Aks '
Labelout(oDust     ) = 'Dust'
Labelout(ow        ) = 'w   '
Labelout(iNO3  + ow) = 'NO3 '
Labelout(iZOO  + ow) = 'ZOO '
Labelout(iDET  + ow) = 'DET '
Labelout(oCHLt + ow) = 'CHL'
Labelout(oNPP  + ow) = 'NPP'
Labelout(oPON  + ow) = 'PON'

DO i = 1, NPHY
   write(Labelout(iPHY(i)  +ow), format_string) 'PHY', i
   write(Labelout(iPHYC(i) +ow), format_string) 'PHC', i
   write(Labelout(omuNet(i)+ow), format_string) 'mu',  i
   write(Labelout(ograz(i) +ow), format_string) 'Gra', i
   write(Labelout(oLno3(i) +ow), format_string) 'LNO', i
   write(Labelout(oSI(i)   +ow), format_string) 'SI' , i
   write(Labelout(oQN(i)   +ow), format_string) 'QN' , i
   write(Labelout(othe(i)  +ow), format_string) 'The' , i
   write(Labelout(iCHL(i)  +ow), format_string) 'CHL', i
ENDDO

DO i = 1, NVAR
   Labelout(oD_VARS(i) + ow) = 'D_'//trim(Labelout(i+ow))
ENDDO

! Write out Varout labels to check:
if (taskid == 0) then
   do i = 1, Nout+ow
      write(6,*) 'Labelout(',i,') = ',trim(Labelout(i))
   enddo
endif

! Initialize parameters
if (taskid==0) write(6,'(I2,1x,A30)') NPar,'parameters to be estimated.'

! Initialize parameters
ParamLabel(imu00)  = 'mu0'
ParamLabel(iaI0)   = 'aI0'
ParamLabel(iKN)    = 'KN'
ParamLabel(iKp)    = 'Kp'
ParamLabel(imz)    = 'mz'
ParamLabel(iwDET)  = 'wDET'

! KN:
! Fennel et al. (2006): 0.007~1.5
! Chai et al. (2002): 0.05~1
! Franks (2009): 0.005~3

MaxValue(imu00)  =  5.1
MinValue(imu00)  =  0.01
  params(imu00)  =  2.4

MaxValue(iKN)  =  3.0
MinValue(iKN)  =  0.05
  params(iKN)  =  0.2

!Half-saturation constant of zooplankton grazing
MaxValue(iKp)  =  3.0
MinValue(iKp)  =  0.01
  params(iKp)  =  0.2

!Zooplankton mortality coefficient
MaxValue(imz)  =  0.25
MinValue(imz)  =  0.01
  params(imz)  =  0.1

!Chl-specific initial C assimilation rate (W m-2)-1 (gChl/molC)-1 d-1
MaxValue(iaI0) = 0.6
MinValue(iaI0) = 0.01
  params(iaI0) = 0.2

MaxValue(iwDET)= 50.
MinValue(iwDET)= 0.01
  params(iwDET)= 0.5

!Log transform:
do i = 1, NPar
   if (MaxValue(i) .le. 0d0) &
   stop "Parameter MaxValues uninitialized!"
   MaxValue(i) = log(MaxValue(i))
   MinValue(i) = log(MinValue(i))
   params(i)   = log(params(i))
enddo
END SUBROUTINE CHOOSE_MODEL
END MODULE PARAM_MOD
