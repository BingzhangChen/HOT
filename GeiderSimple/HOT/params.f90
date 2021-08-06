MODULE PARAM_MOD
USE MPI
USE BIOCOM_MOD
IMPLICIT NONE

! Number of stations to run
integer, parameter      :: Nstn  = 1

! Station name:
character(4), parameter :: Stn(Nstn) = ['HOT']

! Station latitude (Convert latitude to radians)
real,    parameter :: Stn_lat(Nstn)= [22.75*pi/180.]

! timestep 
real               :: dtsec  = 60D0 ! time step in seconds

integer            :: nsave  = 1 ! save the output per day

integer, parameter :: bot_bound  = Dirichlet  ! Option of bottom boundary condition

! Current biological MODEL!
integer, parameter :: Model_ID   = GeiderSimple

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
integer, parameter :: iCHL(NPHY) = (/ (i+iDET, i=1,NPHY) /) 
integer, parameter :: NVAR   = iCHL(NPHY) ! Total number of biological tracers

! Useless, but needed for compilation
integer            :: iVNO3, iVPHY,iCOVNP, iVZOO, iCOVPZ, iCOVNZ, iPHYC(NPHY)  

! Define tracer matrix:
real               :: Vars(NVAR,nlev) = 0d0

! Define the number of sinking tracers:
integer, parameter :: NVsinkterms =  1  ! DET

! Define the index of sinking tracers in the Vars matrix: 
integer, parameter :: Windex(NVsinkterms) = [iDET]

! Indices for additional output variables on top of state variables
integer, parameter :: oLno3(NPHY) = (/ (i+NVAR,       i=1,NPHY) /)
integer, parameter :: oSI(NPHY)   = (/ (i+oLno3(NPHY),i=1,NPHY) /)
integer, parameter :: othe(NPHY)  = (/ (i+oSI(NPHY),  i=1,NPHY) /)
integer, parameter :: oNPP        = othe(NPHY) + 1
integer, parameter :: oNPPr       = oNPP       + 1   !Incubation derived NPP
integer, parameter :: oPON        = oNPPr      + 1
integer, parameter :: oCHLt       = oPON       + 1
integer, parameter :: omuNet(NPHY)= (/ (i+oCHLt,       i=1,NPHY) /)  !Nitrogen-based growth rate
integer, parameter :: omuChl(NPHY)= (/ (i+omuNet(NPHY),i=1,NPHY) /)  !Chl-based growth rate
integer, parameter :: oGraz(NPHY) = (/ (i+omuChl(NPHY),i=1,NPHY) /)

! The diffusion output order must be consistent with the tracer order!
integer, parameter :: oD_VARS(NVAR) = [(oGraz(NPHY) + i, i = 1, NVAR)]
integer, parameter :: Nout          = oD_VARS(NVAR)

! Initialize Varout matrix:
real               :: Varout(Nout, nlev) = 0d0

! Initialize labels for Varout:
character(LEN=10)  :: Labelout(Nout+ow)  = 'Unknown'

! Define parameter indices:
integer, parameter :: imu00   = 1
integer, parameter :: iaI0    = imu00 +1
integer, parameter :: iKN     = iaI0  +1
integer, parameter :: iRDN    = iKN   +1
integer, parameter :: igmax   = iRDN  +1
integer, parameter :: iKp     = igmax +1
integer, parameter :: imz     = iKp   +1
integer, parameter :: iwDET   = imz   +1  ! Sinking rate of detritus
integer, parameter :: ithemax = iwDET +1  ! Maximal Chl:N ratio
integer, parameter :: iRchl   = ithemax+1 ! Basic respiration rate of Chl
integer, parameter :: NPar    = iRchl     ! Total number of parameters
real               :: params(NPar)     = 0d0  ! Define parameters
character(LEN=8)   :: ParamLabel(NPar) = 'Unknown' !Define parameter labels

!  Maxima and minima  of parameter values 
real               :: MaxValue(NPar)   = 0d0, MinValue(NPar) = 0d0
real,    parameter :: Topt(NPHY)       = [(-3. + dble(i), i = 1, NPHY)]
real,    parameter :: alphaG           = 1.1  !Kill-the-winner coefficient

CONTAINS

SUBROUTINE choose_model
IMPLICIT NONE
integer :: rc
character(len=10), parameter :: format_string = "(A3,I0)"
!End of declaration

IF(TASKID==0) WRITE(6,*) 'Geider-Simple phytoplankton model selected!'

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
Labelout(oNPP  + ow) = 'NPPe'
Labelout(oNPPr + ow) = 'NPPr'
Labelout(oPON  + ow) = 'PON'

DO i = 1, NPHY
   write(Labelout(iPHY(i)  +ow), format_string) 'PHY', i
   write(Labelout(omuNet(i)+ow), format_string) 'muN', i
   write(Labelout(omuchl(i)+ow), format_string) 'muH', i
   write(Labelout(ograz(i) +ow), format_string) 'Gra', i
   write(Labelout(oLno3(i) +ow), format_string) 'LNO', i
   write(Labelout(oSI(i)   +ow), format_string) 'SI' , i
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
ParamLabel(iRDN)   = 'RDN'
ParamLabel(igmax)  = 'gmax'
ParamLabel(iKp)    = 'Kp'
ParamLabel(imz)    = 'mz'
ParamLabel(iwDET)  = 'wDET'
ParamLabel(iRchl)  = 'Rchl'
ParamLabel(ithemax)= 'Thetam'

namelist /paramlist/ params, MaxValue, MinValue

! Check whether the namelist file exists.
inquire (file='Model.nml', iostat=rc)

if (rc /= 0) then
    write (6, '(a)') 'Error: namelist file Model.nml does not exist.'
    return
end if
!  open the namelist file and read initial paramter values
open(namlst,file='Model.nml',status='old',action='read')
read(namlst,nml=paramlist)
close(namlst)

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
